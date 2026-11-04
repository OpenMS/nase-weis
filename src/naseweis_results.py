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


def split_sequence(sequence: str) -> list[str]:
    """Tokenize bracketed modifications and attach terminal phosphate markers."""
    if not sequence:
        return []
    # OpenMS accepts bracketed terminal groups as well as p/c shorthand.
    sequence = (
        sequence.replace("[5'-p]", "p").replace("[3'-p]", "p").replace("[3'-c]", "c")
    )
    tokens = re.findall(r"\[[^\[\]]+\]|[ACGUNpc]", sequence)
    if "".join(tokens) != sequence:
        raise ValueError(f"Invalid RNA sequence: {sequence}")
    if tokens[0] == "p":
        tokens = tokens[1:]
        if not tokens:
            raise ValueError("A terminal phosphate must accompany a nucleotide.")
        tokens[0] = "p" + tokens[0]
    if tokens[-1] in ("p", "c"):
        terminal = tokens.pop()
        if not tokens:
            raise ValueError("A terminal phosphate must accompany a nucleotide.")
        tokens[-1] += terminal
    if any(token in ("p", "c") for token in tokens):
        raise ValueError("Phosphate markers are only valid at sequence termini.")
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
