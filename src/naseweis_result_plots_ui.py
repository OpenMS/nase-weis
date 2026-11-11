"""Saved-result coverage and comparative LFQ plots; viewing never reruns tools."""

import streamlit as st
from streamlit.components.v1 import html

from src.naseweis_comparison_plots import condition_scatter, paired_condition_rows
from src.naseweis_coverage import (
    build_coverage,
    render_coverage_html,
    render_coverage_svg,
)
from src.naseweis_results import read_mztab


@st.cache_data(show_spinner=False)
def _saved_tables(path, mtime_ns, ctime_ns, size):
    """Cache only saved native outputs, with each output's own file revision."""
    return read_mztab(path)


def render_coverage_groups(groups, prefix):
    """Render one accession across supplied acquisition/condition tracks."""
    accessions = []
    for _group, tables in groups:
        references = tables.get("NUC")
        if references is not None and "accession" in references:
            accessions.extend(references["accession"].dropna().astype(str))
    accessions = list(dict.fromkeys(accessions))
    if not accessions:
        st.info("No RNA sequence coverage to display for these results.")
        return
    accession = st.selectbox("RNA accession", accessions, key=f"{prefix}-accession")
    block_size = st.select_slider(
        "Bases per line", [25, 50, 75, 100], value=50, key=f"{prefix}-block-size"
    )
    try:
        model = build_coverage(groups, accession)
        plot = render_coverage_html(model, block_size=block_size)
        svg = render_coverage_svg(model, block_size=block_size)
    except (ValueError, KeyError, IndexError, TypeError) as error:
        st.warning(f"Coverage is unavailable: {error}")
        return
    html(plot, height=680, scrolling=True)
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


def render_experiment_coverage(root, samples, prefix):
    """Use available saved mzTab evidence, including identification-only runs."""
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
    render_coverage_groups(_coverage_tracks(available, grouping), f"{prefix}-coverage")


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
    config = {"toImageButtonOptions": {"format": "svg", "filename": "naseweis-lfq-comparison"}}
    st.plotly_chart(
        figure, use_container_width=True, key=f"{prefix}-scatter",
        config=config,
    )
    st.caption("Zoom and hover for RNA identities, means and replicate counts. Use the camera button for SVG export.")
    st.download_button(
        "Download abundance comparison (HTML)",
        figure.to_html(full_html=True, include_plotlyjs=True, config=config),
        "naseweis-lfq-comparison.html", "text/html", key=f"{prefix}-scatter-html",
    )


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
