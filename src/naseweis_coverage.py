"""Offline, sequence-resolved NASE coverage plots for saved acquisitions.

The original NASEWEIS view used packed oligo tracks and shared spectral-count
colors. This renderer keeps those conventions while retaining native mapping
ambiguity and exact modified-sequence identities.
"""

import math
from collections import Counter
from html import escape

import pandas as pd

from src.naseweis_results import nucleotide_token, split_sequence

# Fixed bins share one scale across conditions and across exported figures.
# Colors sample the original coverage view's reversed magma palette.
COUNT_BINS = (
    (1, "1", "#fcfdbf"),
    (2, "2", "#fdc78d"),
    (4, "3–4", "#fb8861"),
    (8, "5–8", "#e65163"),
    (16, "9–16", "#b73779"),
    (32, "17–32", "#822581"),
    (64, "33–64", "#4f127b"),
    (128, "65–128", "#1c1044"),
    (math.inf, "129+", "#000004"),
)


def _coordinates(value, field):
    try:
        numeric = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"Oligonucleotide {field} must be a positive integer.") from error
    if isinstance(value, bool) or not math.isfinite(numeric) or not numeric.is_integer() or numeric < 1:
        raise ValueError(f"Oligonucleotide {field} must be a positive integer.")
    return int(numeric)


def _reference(tables, accession, group):
    table = tables.get("NUC", pd.DataFrame())
    if table.empty:
        return None
    if "accession" not in table:
        raise ValueError(f"Coverage group {group!r} has no nucleotide accession column.")
    selected = table.loc[table["accession"] == accession]
    if selected.empty:
        return None
    sequences = selected.get("opt_sequence", pd.Series(dtype="object"))
    valid = [sequence for sequence in sequences if isinstance(sequence, str) and sequence]
    if not valid:
        return None
    tokenized = [split_sequence(sequence) for sequence in valid]
    if any(tokens != tokenized[0] for tokens in tokenized[1:]):
        raise ValueError(f"Conflicting reference sequences for accession {accession!r} in {group!r}.")
    return tokenized[0]


def _assignment_counts(tables):
    """Count each native sequence/scan pair once within one acquisition."""
    table = tables.get("OSM")
    if table is None or "sequence" not in table:
        return None
    if table.empty:
        return Counter()
    counts, seen = Counter(), set()
    has_scans = "spectra_ref" in table
    for index, row in enumerate(table.to_dict("records")):
        label = row["sequence"]
        if not isinstance(label, str) or not label:
            raise ValueError("Coverage spectrum assignments require a native RNA sequence.")
        scan = row.get("spectra_ref") if has_scans else None
        # A missing scan reference cannot be deduplicated against another row.
        scan = ("scan", scan) if isinstance(scan, str) and scan else ("row", index)
        pair = (label, scan)
        if pair not in seen:
            counts[label] += 1
            seen.add(pair)
    return counts


def build_coverage(groups, accession):
    """Build coverage from ``[(group_label, mzTab_tables), ...]``.

    Repeated group labels pool acquisitions. Every distinct exported position
    remains visible, including ambiguous mappings; duplicate OLI records do not
    duplicate evidence. Counts use exact modified sequences and distinct native
    sequence/scan pairs within each acquisition, then sum across acquisitions.
    A scan supporting multiple sequences and a sequence mapping to multiple
    positions are retained as shared evidence, not resolved into new identities.
    """
    if not isinstance(accession, str) or not accession:
        raise ValueError("Choose an RNA accession for sequence coverage.")
    entries = list(groups)
    if not entries:
        raise ValueError("No acquisition results are available for sequence coverage.")
    reference = None
    references = []
    for name, tables in entries:
        if not isinstance(name, str) or not name.strip() or not isinstance(tables, dict):
            raise ValueError("Coverage groups require a text label and native mzTab tables.")
        tokens = _reference(tables, accession, name)
        if tokens is not None:
            if reference is not None and tokens != reference:
                raise ValueError(f"Conflicting reference sequences for accession {accession!r}; these acquisitions cannot share a coverage axis.")
            reference = tokens
        references.append(tokens)
    if not reference:
        raise ValueError("No reference sequence was exported for this accession.")
    grouped, warnings = {}, []
    for (name, tables), own_reference in zip(entries, references):
        state = grouped.setdefault(name, {"name": name, "acquisitions": 0, "hits": {}, "assignments": 0})
        state["acquisitions"] += 1
        table = tables.get("OLI", pd.DataFrame())
        if table.empty:
            continue
        if "accession" not in table:
            raise ValueError(f"Coverage group {name!r} has no oligonucleotide accession column.")
        hits = table.loc[table["accession"] == accession]
        if hits.empty:
            continue
        if own_reference is None:
            raise ValueError(f"No reference sequence was exported for mapped accession {accession!r} in {name!r}.")
        if not {"sequence", "start", "end"}.issubset(hits):
            raise ValueError(f"Coverage group {name!r} has incomplete oligonucleotide coordinates.")
        counts = _assignment_counts(tables)
        if counts is None:
            warning = f"Spectrum counts are unavailable for one or more acquisitions in {name}; bars indicate exported mappings only."
            if warning not in warnings:
                warnings.append(warning)
            state["assignments"] = None
        observed, mapped_labels = set(), set()
        for row in hits.to_dict("records"):
            start, end = _coordinates(row["start"], "start"), _coordinates(row["end"], "end")
            if not start <= end <= len(reference):
                raise ValueError(f"Oligonucleotide coordinates {start}–{end} exceed the reference sequence.")
            label = row["sequence"]
            if not isinstance(label, str) or not label:
                raise ValueError("Coverage oligonucleotides require native RNA sequences.")
            tokens = split_sequence(label)
            if len(tokens) != end - start + 1:
                raise ValueError(f"Oligonucleotide {label!r} does not match its exported coordinate span {start}–{end}.")
            identity = (label, start, end)
            if identity in observed:
                continue
            observed.add(identity)
            mapped_labels.add(label)
            count = None if counts is None else counts[label]
            existing = state["hits"].get(identity)
            if existing is None:
                state["hits"][identity] = {
                    "sequence": label, "start": start, "end": end, "count": count,
                    "modifications": [
                        {"position": start + offset, "token": nucleotide_token(token)}
                        for offset, token in enumerate(tokens)
                        if nucleotide_token(token) not in {"A", "C", "G", "U", "N"}
                    ],
                }
            elif count is None or existing["count"] is None:
                existing["count"] = None
            else:
                existing["count"] += count
        if state["assignments"] is not None:
            state["assignments"] += sum(counts[label] for label in mapped_labels)
    result_groups, covered = [], set()
    for state in grouped.values():
        hits = sorted(state["hits"].values(), key=lambda hit: (hit["start"], hit["end"], hit["sequence"]))
        own_covered = {position for hit in hits for position in range(hit["start"], hit["end"] + 1)}
        covered.update(own_covered)
        result_groups.append({
            **state, "hits": hits, "covered_bases": len(own_covered),
            "coverage": len(own_covered) / len(reference),
            "oligonucleotides": len({hit["sequence"] for hit in hits}),
        })
    return {
        "accession": accession, "reference": reference, "length": len(reference),
        "covered_bases": len(covered), "coverage": len(covered) / len(reference),
        "groups": result_groups, "warnings": warnings,
    }


def _count_color(count):
    if count is None or count == 0:
        return "#e4e7ec"
    return next(color for upper, _, color in COUNT_BINS if count <= upper)


def _text(x, y, value, *, size=12, weight="normal", anchor="start", fill="#202532", title=None):
    tooltip = f"<title>{escape(str(title))}</title>" if title is not None else ""
    return (
        f'<text x="{x:g}" y="{y:g}" font-size="{size:g}" font-weight="{weight}" '
        f'text-anchor="{anchor}" fill="{fill}">{tooltip}{escape(str(value))}</text>'
    )


def _lanes(hits, left, right):
    ends, placements = [], []
    for hit in hits:
        start, end = max(left, hit["start"]), min(right, hit["end"])
        if start > end:
            continue
        lane = next((index for index, occupied in enumerate(ends) if occupied < start), len(ends))
        if lane == len(ends):
            ends.append(end)
        else:
            ends[lane] = end
        placements.append((hit, start, end, lane))
    return placements, max(1, len(ends))


def render_coverage_svg(model, *, block_size=50):
    """Export packed, shared-scale condition tracks as a standalone SVG."""
    if isinstance(block_size, bool) or not isinstance(block_size, int) or not 10 <= block_size <= 100:
        raise ValueError("Coverage blocks must contain between 10 and 100 bases.")
    reference, groups = model["reference"], model["groups"]
    cell, left = 12, 150
    width = max(700, left + min(block_size, len(reference)) * cell + 30)
    body, y = [], 32
    body.append(_text(12, y, f"{model['accession']} sequence coverage", size=19, weight="bold"))
    y += 23
    body.append(_text(12, y, f"{model['covered_bases']} / {model['length']} bases mapped ({model['coverage']:.1%}); native 1-based inclusive coordinates", size=12))
    y += 27
    for offset in range(0, len(reference), block_size):
        right = min(offset + block_size, len(reference))
        body.append(_text(12, y, f"Positions {offset + 1}–{right}", size=12, weight="bold"))
        y += 22
        body.append(_text(12, y, "Reference", size=12, weight="bold"))
        for index, token in enumerate(reference[offset:right]):
            position = offset + index + 1
            base = nucleotide_token(token)
            glyph = "*" if "[" in base else base
            x = left + (index + 0.5) * cell
            body.append(_text(x, y, glyph, size=13, anchor="middle", weight="bold", title=f"Position {position}: {token}"))
            if position % 5 == 0 or index == 0 or position == right:
                body.append(_text(x, y + 14, position, size=9, anchor="middle", fill="#626a78"))
        y += 32
        for group in groups:
            placements, count = _lanes(group["hits"], offset + 1, right)
            name = group["name"]
            short_name = name if len(name) <= 22 else name[:19] + "…"
            body.append(_text(12, y + 11, short_name, size=12, weight="bold", title=name))
            body.append(_text(12, y + 26, f"{group['coverage']:.1%} mapped · n={group['acquisitions']}", size=10, fill="#626a78"))
            if not placements:
                body.append(_text(left, y + 12, "No mapped oligonucleotides in this block", size=10, fill="#626a78"))
            for hit, start, end, lane in placements:
                x = left + (start - offset - 1) * cell
                top = y + lane * 19
                count_label = "unknown" if hit["count"] is None else str(hit["count"])
                tooltip = f"{hit['sequence']} · positions {hit['start']}–{hit['end']} · {count_label} sequence–spectrum assignments · {name}"
                body.append(
                    f'<rect class="mapped-oligo" x="{x:g}" y="{top:g}" width="{(end-start+1)*cell-1:g}" height="14" '
                    f'rx="2" fill="{_count_color(hit["count"])}" stroke="#7d8490" stroke-width="0.35" '
                    f'data-start="{hit["start"]}" data-end="{hit["end"]}" data-sequence="{escape(hit["sequence"], quote=True)}" '
                    f'data-count="{count_label}"><title>{escape(tooltip)}</title></rect>'
                )
                for modification in hit["modifications"]:
                    position = modification["position"]
                    if start <= position <= end:
                        marker_x = left + (position - offset - 0.5) * cell
                        title = f"Position {position}: {modification['token']} in {hit['sequence']} (native annotation)"
                        body.append(f'<circle class="modification" cx="{marker_x:g}" cy="{top+7:g}" r="3.3" fill="#00a6a6" stroke="#fff" stroke-width="0.8"><title>{escape(title)}</title></circle>')
            y += max(34, count * 19 + 9)
        body.append(f'<line x1="12" y1="{y}" x2="{width-20}" y2="{y}" stroke="#e0e3e9"/>')
        y += 28
    body.append(_text(12, y, "Sequence–spectrum assignments", size=12, weight="bold"))
    y += 15
    legend = list(COUNT_BINS) + [(0, "0 / unknown", "#e4e7ec")]
    for index, (_, label, color) in enumerate(legend):
        x = 12 + index * 63
        body.append(f'<rect x="{x}" y="{y}" width="59" height="14" fill="{color}" stroke="#7d8490" stroke-width="0.35"/>')
        body.append(_text(x + 29.5, y + 28, label, size=10, anchor="middle"))
    y += 51
    notes = [
        "Colors share fixed count bins across all tracks. Gray = zero or unavailable count; hover for exact counts.",
        "Teal dots mark annotated nucleotide modifications; terminal p/c groups are not modification calls.",
        "Native D/I/Y codes also mark modified nucleotides; a leading * denotes terminal phosphorothioate, not a base.",
        "Repeated/ambiguous mappings retain shared evidence at every native position; counts are not independent totals.",
        "n = acquisitions pooled in the track. Reference * = bracket-annotated nucleotide; hover for its exact identity.",
    ] + model.get("warnings", [])
    for note in notes:
        # Wrap explanatory text at words so vector exports stay self-contained.
        line = ""
        for word in note.split():
            if len(line) + len(word) > max(75, (width - 30) // 6):
                body.append(_text(12, y, line, size=10, fill="#626a78"))
                y += 15
                line = ""
            line += (" " if line else "") + word
        body.append(_text(12, y, line, size=10, fill="#626a78"))
        y += 18
    title = escape(f"{model['accession']} NASE sequence coverage")
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{y+8}" viewBox="0 0 {width} {y+8}" role="img" aria-labelledby="coverage-title" '
        'font-family="Arial, Helvetica, sans-serif">'
        f'<title id="coverage-title">{title}</title><rect width="100%" height="100%" fill="#fff"/>'
        + "".join(body) + "</svg>"
    )


def render_coverage_html(model, *, block_size=50):
    """Wrap the SVG in a self-contained page with horizontally scrollable tracks."""
    svg = render_coverage_svg(model, block_size=block_size)
    return (
        '<!doctype html><html><head><meta charset="utf-8"><title>'
        + escape(f"{model['accession']} sequence coverage")
        + '</title><style>body{margin:8px;background:#fff}.coverage{overflow-x:auto}</style>'
        + '</head><body><div class="coverage">' + svg + "</div></body></html>"
    )
