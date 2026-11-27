"""Saved-result coverage and comparative LFQ plots; viewing never reruns tools."""

import pandas as pd
import streamlit as st

from src.naseweis_insight import show_plot
from src.naseweis_comparison_plots import condition_scatter, paired_condition_rows
from src.naseweis_coverage import (
    build_coverage,
    coverage_plot,
    render_coverage_html,
    render_coverage_svg,
)
from src.naseweis_results import read_mztab


@st.cache_data(show_spinner=False)
def _saved_tables(path, mtime_ns, ctime_ns, size):
    """Cache only saved native outputs, with each output's own file revision."""
    return read_mztab(path)


def coverage_accessions(groups):
    """The RNA accessions of ``[(group, tables), ...]``, in the order they first appear."""
    accessions = []
    for _group, tables in groups:
        references = tables.get("NUC")
        if references is not None and "accession" in references:
            accessions.extend(references["accession"].dropna().astype(str))
    return list(dict.fromkeys(accessions))


def rna_summary(groups):
    """One row per RNA of ``groups``: its length, the oligonucleotides and
    spectra found for it and the share of its residues they cover (pooled over
    the groups). An RNA whose coverage cannot be built keeps empty figures."""
    rows = []
    for index, accession in enumerate(coverage_accessions(groups)):
        row = {"rna_idx": index, "accession": accession, "length": None,
               "oligonucleotides": None, "matches": None, "coverage": None}
        try:
            model = build_coverage([("all", tables) for _group, tables in groups], accession)
        except (ValueError, KeyError, IndexError, TypeError):
            rows.append(row)
            continue
        group = model["groups"][0] if model["groups"] else {}
        row.update(length=model["length"], oligonucleotides=group.get("oligonucleotides", 0),
                   matches=group.get("assignments"), coverage=round(100 * model["coverage"], 1))
        rows.append(row)
    return pd.DataFrame(rows, columns=["rna_idx", "accession", "length", "oligonucleotides",
                                       "matches", "coverage"])


def render_coverage_groups(groups, prefix, accession=None):
    """Render one accession across supplied acquisition/condition tracks.

    ``accession`` (picked elsewhere, in the RNA table) is shown without the
    accession selector; an accession these results do not hold falls back to it.
    """
    accessions = coverage_accessions(groups)
    if not accessions:
        st.info("No RNA sequence coverage to display for these results.")
        return
    if accession in accessions:
        st.markdown(f"**{accession}**")
    else:
        accession = st.selectbox("RNA accession", accessions, key=f"{prefix}-accession")
    block_size = st.select_slider(
        "Bases per line", [25, 50, 75, 100], value=50, key=f"{prefix}-block-size"
    )
    try:
        model = build_coverage(groups, accession)
        figure = coverage_plot(model, block_size=block_size)
        plot = render_coverage_html(model, block_size=block_size)
        svg = render_coverage_svg(model, block_size=block_size)
    except (ValueError, KeyError, IndexError, TypeError) as error:
        st.warning(f"Coverage is unavailable: {error}")
        return
    show_plot(figure, f"{prefix}-plot", "naseweis-coverage")
    st.download_button(
        "Download coverage (HTML)", plot, "naseweis-coverage.html", "text/html",
        key=f"{prefix}-html",
    )
    st.download_button(
        "Download coverage (SVG)", svg, "naseweis-coverage.svg", "image/svg+xml",
        key=f"{prefix}-svg",
    )


def _coverage_tracks(available, grouping):
    """Keep an unassigned acquisition separate from every named condition."""
    reserved = {sample.get("condition", "").strip() for sample, _ in available}
    groups = []
    for sample, tables in available:
        name = sample["file"]
        if grouping == "Condition":
            name = sample.get("condition", "").strip()
            if not name:
                name = f"Unassigned · {sample['file']}"
                while name in reserved:
                    name += " (unassigned acquisition)"
                reserved.add(name)
        groups.append((name, tables))
    return groups


def experiment_coverage_tables(root, samples):
    """``[(sample, tables), ...]`` of the acquisitions with saved mzTab
    evidence, including identification-only runs; warns for unreadable ones."""
    from src.naseweis_experiment import sample_results

    available = []
    for sample in samples:
        try:
            folder = sample_results(root, sample["id"]) / "mztab_results"
            path = folder / "search.mzTab"
            if not path.is_file():
                path = folder / "tab_out"
            if not path.is_file():
                continue
            stat = path.stat()
            tables = _saved_tables(
                str(path.resolve()), stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size
            )
        except (OSError, ValueError, KeyError) as error:
            st.warning(f"Could not read coverage for {sample['file']}: {error}")
        else:
            available.append((sample, tables))
    return available


def render_experiment_coverage(root, samples, prefix, accession=None):
    """Coverage of the experiment's acquisitions, as condition or acquisition tracks."""
    available = experiment_coverage_tables(root, samples)
    if not available:
        st.info("No RNA sequence coverage to display for these results.")
        return
    st.caption(
        f"Saved identification evidence from {len(available)} of {len(samples)} acquisitions. "
        "Colors show summed spectrum matches, not LFQ abundance. "
        "Condition tracks pool identification evidence across their acquisitions; "
        "use acquisition tracks to inspect unequal group sizes."
    )
    if len(available) > 1:
        assigned = all(sample.get("condition", "").strip() for sample, _ in available)
        grouping = st.radio(
            "Coverage tracks", ["Condition", "Acquisition"],
            index=0 if assigned else 1, horizontal=True, key=f"{prefix}-grouping",
        )
    else:
        grouping = "Acquisition"
    render_coverage_groups(_coverage_tracks(available, grouping), f"{prefix}-coverage", accession)


def render_lfq_scatter(contrasts, numerator, denominator, prefix, *, measure, normalization):
    """Plot the same observed means and direction used by the comparison table."""
    st.markdown("**RNA abundance comparison**")
    scale = st.radio(
        "Scatter axes", ["Logarithmic", "Linear"], horizontal=True,
        key=f"{prefix}-scatter-scale",
    )
    try:
        paired = paired_condition_rows(contrasts)
    except ValueError as error:
        st.warning(f"The abundance comparison plot is unavailable: {error}")
        return
    omitted = len(contrasts) - len(paired)
    st.caption(
        f"{len(paired)} oligonucleotides measured in both groups; {omitted} without "
        "positive measurements in both groups are omitted. No zero substitution or "
        "missing-value imputation is used. The diagonal marks equal abundance. "
        "Modified markers reflect native RNA modification annotations, not inferred modification sites."
    )
    if paired.empty:
        st.info("No oligonucleotides have positive measurements in both groups to plot.")
        return
    selected = st.multiselect(
        "Label oligonucleotides", list(paired.index), key=f"{prefix}-scatter-labels"
    )
    try:
        figure = condition_scatter(
            contrasts, numerator, denominator, measure=measure,
            normalization=normalization,
            scale="log" if scale == "Logarithmic" else "linear", labels=selected,
        )
    except ValueError as error:
        st.warning(f"The abundance comparison plot is unavailable: {error}")
        return
    show_plot(figure, f"{prefix}-scatter", "naseweis-lfq-comparison")
    st.caption("Hover for RNA identities, means and replicate counts; zoom with the mouse. The plot's menu exports SVG.")


def render_acquisition_scatter(matrix, design, labels, prefix, *, measure, normalization):
    """Offer a measured-acquisition comparison when there is only one condition."""
    from src.naseweis_comparison import compare_conditions

    if len(matrix.columns) < 2:
        return
    st.markdown("**Acquisition comparison**")
    numerator = st.selectbox(
        "Y acquisition", list(matrix.columns), index=1, format_func=labels.get,
        key=f"{prefix}-scatter-y",
    )
    denominator = st.selectbox(
        "X acquisition", [column for column in matrix.columns if column != numerator],
        format_func=labels.get, key=f"{prefix}-scatter-x",
    )
    selected = [denominator, numerator]
    comparison_design = design.loc[selected].copy()
    comparison_design["condition"] = selected
    contrasts = compare_conditions(
        matrix[selected], comparison_design, numerator, denominator
    )
    st.caption("Each point compares two acquisitions directly; these are not condition means.")
    render_lfq_scatter(
        contrasts, labels[numerator], labels[denominator], prefix,
        measure=measure, normalization=normalization,
    )
