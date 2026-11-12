"""Read NASE mzTab output and render sequence coverage without network access.

Authors: Hendrik Weisser, Samuel Wein (original coverage visualization).
License: BSD 3-clause.
"""

import re
from collections import defaultdict
from html import escape
from pathlib import Path

import pandas as pd

HEADERS = {"NUH": "NUC", "OLH": "OLI", "OSH": "OSM"}
# NASequence emits D (dihydrouridine), I (inosine) and Y (pseudouridine)
# without brackets. Terminal groups are decorations, not extra nucleotides.
_NUCLEOTIDE = r"[ACGUNDIY]|\[(?![35]'-)[^\[\]]+\]"
_FIVE_PRIME = r"p|\*|\[5'-[^\[\]]+\]"
_THREE_PRIME = r"p|c|\[3'-[^\[\]]+\]"
_NUCLEOTIDE_TOKEN = re.compile(
    rf"(?:{_FIVE_PRIME})?(?P<nucleotide>{_NUCLEOTIDE})(?:{_THREE_PRIME})?"
)
_TERMINAL_ALIASES = {
    "[5'-p]": "p", "[5'-p*]": "*", "[3'-p]": "p", "[3'-c]": "c",
}


def read_mztab(path: str | Path) -> dict[str, pd.DataFrame]:
    """Read NASE's nucleotide, oligonucleotide and spectrum-match sections.

    Empty sections retain their columns. Metadata and comment rows can occur
    between sections; they do not become table rows. NASE coordinates are 1-based.
    """
    headers, rows = {}, defaultdict(list)
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            tag = fields[0]
            if tag in HEADERS:
                headers[HEADERS[tag]] = fields[1:]
            elif tag in HEADERS.values():
                if tag not in headers or len(fields) - 1 != len(headers[tag]):
                    raise ValueError(
                        f"Invalid {tag} row at line {number} in {Path(path).name}."
                    )
                rows[tag].append(fields[1:])
    tables = {}
    for tag, columns in headers.items():
        table = pd.DataFrame(rows[tag], columns=columns).replace("null", pd.NA)
        for column in [
            "start",
            "end",
            "charge",
            "retention_time",
            "exp_mass_to_charge",
            "calc_mass_to_charge",
            "sequence_coverage",
            "search_engine_score[1]",
            "search_engine_score[2]",
        ]:
            if column in table:
                table[column] = pd.to_numeric(table[column], errors="coerce")
        tables[tag] = table
    return tables


def nucleotide_token(token: str) -> str:
    """Return one native nucleotide identity without terminal decorations.

    ``split_sequence`` attaches terminal groups to its outer nucleotide tokens.
    This helper removes only those decorations, preserving modified nucleotide
    identity such as ``D``, ``Y`` or ``[m5C]`` for labels and site markers.
    """
    match = _NUCLEOTIDE_TOKEN.fullmatch(token) if isinstance(token, str) else None
    if match is None:
        raise ValueError(f"Invalid RNA nucleotide token: {token}")
    return match["nucleotide"]


def split_sequence(sequence: str) -> list[str]:
    """Tokenize native RNA identities with terminal groups on the outer bases.

    In addition to bracketed nucleotide annotations, OpenMS emits D/I/Y as
    single letters and uses a leading ``*`` for a 5′ phosphorothioate. Known
    bracketed terminal aliases canonicalize to p/*/c. Other explicit 5′/3′
    terminal groups retain their labels without increasing nucleotide count.
    """
    if not isinstance(sequence, str):
        raise ValueError("An RNA sequence must be text.")  # noqa: TRY004
    if sequence == "":
        return []
    tokens = re.findall(r"\[[^\[\]]+\]|[ACGUNDIYpc*]", sequence)
    if "".join(tokens) != sequence:
        raise ValueError(f"Invalid RNA sequence: {sequence}")
    five_prime = ""
    if re.fullmatch(_FIVE_PRIME, tokens[0]):
        terminal = tokens.pop(0)
        five_prime = _TERMINAL_ALIASES.get(terminal, terminal)
    three_prime = ""
    if tokens and re.fullmatch(_THREE_PRIME, tokens[-1]):
        terminal = tokens.pop()
        three_prime = _TERMINAL_ALIASES.get(terminal, terminal)
    if not tokens:
        raise ValueError("A terminal group must accompany a nucleotide.")
    if any(re.fullmatch(_NUCLEOTIDE, token) is None for token in tokens):
        raise ValueError("Terminal groups are only valid at the sequence termini.")
    tokens[0] = five_prime + tokens[0]
    tokens[-1] += three_prime
    return tokens


def match_table(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Format OSM rows, preserving ambiguous sequence-to-accession mappings."""
    matches = tables.get("OSM", pd.DataFrame()).copy()
    if matches.empty:
        return matches
    oligos = tables.get("OLI", pd.DataFrame())
    if not oligos.empty:
        accessions = oligos.groupby("sequence")["accession"].agg(
            lambda values: "; ".join(sorted(set(values.dropna().astype(str))))
        )
        matches.insert(0, "accession", matches["sequence"].map(accessions))
    names = {
        "accession": "Accession",
        "sequence": "Sequence",
        "retention_time": "Retention time (s)",
        "exp_mass_to_charge": "Experimental m/z",
        "calc_mass_to_charge": "Calculated m/z",
        "charge": "Charge",
        "search_engine_score[1]": "Search score",
        "spectra_ref": "Spectrum",
    }
    result = matches[[key for key in names if key in matches]].rename(columns=names)
    # Score 2 is -1 when NASE did not compute a q-value. Use the output, not
    # the current UI settings: these may describe a different run.
    q_values = matches.get("search_engine_score[2]")
    if q_values is not None and q_values.ge(0).any():
        result.insert(
            min(7, len(result.columns)),
            "q-value (%)",
            q_values.where(q_values.ge(0)) * 100,
        )
    # Keep native diagnostics (e.g. precursor intensity and scoring metadata)
    # visible and downloadable instead of discarding new optional columns.
    for column in matches.columns:
        if column.startswith("opt_"):
            result[column] = matches[column]
    return result


def coverage_html(tables: dict[str, pd.DataFrame], accession: str) -> str:
    """Render reference bases and every mapped oligo, with escaped tooltips.

    Each hit gets a separate row so overlapping and modified oligos remain
    visible. The plot uses only the result file and works entirely offline.
    """
    references = tables.get("NUC", pd.DataFrame())
    oligos = tables.get("OLI", pd.DataFrame())
    matches = tables.get("OSM", pd.DataFrame())
    reference = references.loc[references["accession"] == accession].iloc[0]
    sequence = reference.get("opt_sequence", "")
    if not isinstance(sequence, str) or not sequence:
        return "<p>No reference sequence was exported for this accession.</p>"
    bases = split_sequence(sequence)
    if not bases:
        return "<p>No reference sequence was exported for this accession.</p>"
    hits = oligos.loc[oligos["accession"] == accession].sort_values(["start", "end"])
    counts = (
        matches["sequence"].value_counts()
        if not matches.empty
        else pd.Series(dtype=int)
    )
    covered = set()
    mapped = []
    for _, hit in hits.iterrows():
        start, end = int(hit["start"]), int(hit["end"])
        if not 1 <= start <= end <= len(bases):
            raise ValueError(
                f"Oligonucleotide coordinates {start}–{end} exceed the reference sequence."
            )
        covered.update(range(start, end + 1))
        label = str(hit["sequence"])
        count = int(counts.get(label, 0))
        tooltip = escape(
            f"{label} · positions {start}–{end} · {count} spectrum matches", quote=True
        )
        mapped.append((start, end, label, count, tooltip))
    blocks = []
    for offset in range(0, len(bases), 50):
        block_bases = bases[offset : offset + 50]
        cells = "".join(
            f'<span title="Position {offset + i}">{escape(base)}</span>'
            for i, base in enumerate(block_bases, 1)
        )
        lanes = []
        for start, end, label, count, tooltip in mapped:
            left, right = max(start, offset + 1), min(end, offset + len(block_bases))
            if left <= right:
                lanes.append(
                    f'<div class="lane"><span class="hit" style="grid-column:{left - offset}/{right - offset + 1}" '
                    f'title="{tooltip}">{escape(label)} · {count}</span></div>'
                )
        blocks.append(
            f'<section style="--columns:{len(block_bases)}"><p class="positions">'
            f"Positions {offset + 1}–{offset + len(block_bases)}</p>"
            f'<div class="reference">{cells}</div>{"".join(lanes)}</section>'
        )
    return f"""<!doctype html><html><head><meta charset="utf-8"><style>
:root {{color-scheme:light dark}} body {{font:14px system-ui,sans-serif; margin:8px}}
.scroll {{overflow-x:auto; padding-bottom:12px}}
section {{min-width:700px;margin:24px 0}} .positions {{font-size:12px;opacity:.7}}
.reference,.lane {{display:grid;grid-template-columns:repeat(var(--columns),minmax(0,1fr));gap:0;margin:5px 0}}
.reference span {{text-align:center;font-family:monospace;font-weight:bold}}
.hit {{border-radius:4px;background:#29379b;color:#fff;padding:3px 5px;font:12px monospace;overflow:hidden;white-space:nowrap}}
</style></head><body><p><strong>{escape(accession)}</strong> · {len(covered) / len(bases):.1%} sequence coverage · {len(hits)} mapped oligonucleotides</p>
<p>Hover over a bar for its sequence, reference positions and spectrum count. Sequences wrap every 50 bases; bars crossing a boundary continue in the next block.</p>
<div class="scroll">{"".join(blocks)}</div></body></html>"""
