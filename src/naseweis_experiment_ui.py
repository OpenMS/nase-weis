"""Configure RNA experiments and inspect their saved, cross-sample results."""

import hashlib
import json
import math
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from src.naseweis_tools import enabled
from src.workflow.ExperimentalDesign import (
    default_design,
    design_path,
    load_design,
    save_design,
)

NORMALIZATIONS = {"total": "Total area", "none": "None"}
MEASURES = {"intensity": "Native feature area", "model_area": "Fitted model area"}


def _validate_design_edits(state):
    for changes in state.get("edited_rows", {}).values():
        for field in ("replicate", "fraction"):
            if field not in changes:
                continue
            value = changes[field]
            try:
                number = float(value)
                valid = (
                    not isinstance(value, bool)
                    and math.isfinite(number)
                    and number.is_integer()
                    and number > 0
                )
            except (TypeError, ValueError, OverflowError):
                valid = False
            if not valid:
                raise ValueError(
                    f"Sample {field} must be a positive integer; received {value!r}."
                )


def _design_ready(workflow, selected):
    """Validate raw design input before the generic editor can coerce numbers."""
    from src.naseweis_experiment import validate_design_source

    version = f"_design-{Path(workflow.workflow_dir).name}-mzML-files"
    path = design_path(workflow.workflow_dir)
    stored_invalid = False
    try:
        if path.is_file():
            try:
                validate_design_source(path.read_bytes())
                load_design(workflow.workflow_dir, selected)
            except ValueError:
                stored_invalid = True
                raise
        uploaded = st.session_state.get(f"{version}-import")
        if uploaded is not None:
            validate_design_source(uploaded.getvalue())
        for key, state in st.session_state.items():
            if (
                key.startswith(f"{version}-")
                and isinstance(state, dict)
                and "edited_rows" in state
            ):
                _validate_design_edits(state)
    except (OSError, ValueError) as error:
        st.error(f"Could not use the sample design: {error}")
        label = (
            "Reset sample design" if stored_invalid else "Discard invalid design input"
        )
        if st.button(label, key=f"{version}-recover"):
            if stored_invalid:
                save_design(workflow.workflow_dir, default_design(selected), selected)
            for key in list(st.session_state):
                if key.startswith(f"{version}-"):
                    del st.session_state[key]
            st.session_state[version] = st.session_state.get(version, 0) + 1
            st.rerun()
        return False
    return True


def configure_experiment_inputs(workflow):
    """Migrate legacy single-file selections before rendering the multi-select."""
    key = "mzML-files"
    widget = f"{workflow.parameter_manager.param_prefix}{key}"
    value = st.session_state.get(widget, workflow.params.get(key, []))
    selected = [value] if isinstance(value, str) and value else value
    if selected is None or selected == "":
        selected = []
    if not isinstance(selected, list) or any(
        not isinstance(path, str) for path in selected
    ):
        st.error(
            "The saved acquisition selection is invalid. Select the input files again."
        )
        return
    if value != selected or workflow.params.get(key) != selected:
        st.session_state[widget] = selected
        if not workflow.parameter_manager.save_parameters():
            return
    workflow.params[key] = selected
    workflow.ui.params[key] = selected
    workflow.ui.select_input_file(
        key,
        name="MS/MS acquisitions",
        multiple=True,
        reactive=True,
        upload={"file_types": ["mzML", "raw", "d"]},
    )
    selected = st.session_state.get(widget, selected)
    workflow.params[key] = list(selected)
    workflow.ui.params[key] = list(selected)
    st.caption(
        "Every selected acquisition is included in the next experiment. "
        "Conditions and biological replicate numbers are required for LFQ and "
        "optional for identification-only experiments. Acquisitions with the "
        "same condition and replicate are technical repeats."
    )
    if not _design_ready(workflow, selected):
        return
    workflow.ui.experimental_design(key, name="Samples and conditions", suggest=True)
    if selected:
        try:
            design = load_design(workflow.workflow_dir, selected)
        except (OSError, ValueError, KeyError) as error:
            st.error(f"Could not read the sample design: {error}")
            return
        lfq_enabled = enabled(
            st.session_state.get(
                f"{workflow.parameter_manager.param_prefix}enable-lfq",
                workflow.params.get("enable-lfq", False),
            )
        )
        if lfq_enabled and design["condition"].str.strip().eq("").any():
            st.warning(
                "Condition comparisons require a condition for each selected acquisition."
            )
        if lfq_enabled and design["fraction"].ne(1).any():
            st.warning(
                "RNA experiment comparisons currently require unfractionated samples (fraction 1)."
            )


def _comparison_control(workflow, key, label, choices, default, *, visible=True):
    widget = f"{workflow.parameter_manager.param_prefix}{key}"
    value = st.session_state.get(widget, workflow.params.get(key, default))
    valid = isinstance(value, str) and value in choices
    if not valid:
        st.error(
            f"Unsupported {label.lower()} setting: {value!r}. Correct it before running the experiment."
        )
    if not visible:
        return
    if not isinstance(value, str):
        return
    options = list(choices) + ([] if valid else [value])
    if widget not in st.session_state:
        st.session_state[widget] = value
    st.selectbox(
        label,
        options,
        format_func=lambda option: choices.get(option, f"Unsupported: {option}"),
        key=widget,
        on_change=workflow.parameter_manager.save_parameters,
    )
    workflow.params[key] = st.session_state[widget]
    workflow.ui.params[key] = st.session_state[widget]
    workflow.parameter_manager.save_parameters()


def configure_comparison(workflow):
    """Configure an enabled LFQ experiment without changing its enable toggle."""
    st.markdown("**Compare acquisitions**")
    st.caption(
        "Select two or more acquisitions to compare RNA abundances. Total-area "
        "normalization divides each RNA area by the sum of quantified RNA areas "
        "in that acquisition. None retains the extracted areas. "
        "Missing measurements remain missing."
    )
    _comparison_control(
        workflow, "lfq-normalization", "Normalization", NORMALIZATIONS, "total"
    )
    _comparison_control(
        workflow,
        "lfq-measure",
        "Area measure",
        MEASURES,
        "intensity",
        visible=bool(st.session_state.get("advanced", False)),
    )
    if not st.session_state.get("advanced", False):
        measure = workflow.params.get("lfq-measure", "intensity")
        if isinstance(measure, str) and measure in MEASURES:
            st.caption(
                f"Area measure: {MEASURES[measure]}. Enable Advanced to change it."
            )


def _sample_labels(samples):
    return {
        sample["id"]: (
            f"{sample['file']} · {sample.get('condition') or 'Unassigned'} "
            f"· replicate {sample.get('replicate', 1)}"
        )
        for sample in samples
    }


def _download_frame(label, frame, filename, *, index_label):
    st.download_button(
        label,
        frame.to_csv(sep="\t", index_label=index_label, na_rep="").encode(),
        filename,
        "text/tab-separated-values",
    )


def _render_comparison(root, samples, prefix):
    from src.naseweis_comparison import compare_conditions, load_comparison

    try:
        comparison = load_comparison(root)
    except (OSError, ValueError, KeyError) as error:
        st.warning(f"The experiment abundance matrix is not available: {error}")
        return
    raw, matrix = comparison["raw"], comparison["abundances"]
    design, metadata = comparison["design"], comparison["metadata"]
    fields = ("file", "path", "condition", "replicate", "fraction")
    if list(design.index) != [sample["id"] for sample in samples] or any(
        design.loc[sample["id"], field] != sample[field]
        for sample in samples
        for field in fields
    ):
        st.error(
            "The saved abundance matrix does not match this experiment's acquisitions and sample design. Re-run the complete experiment before comparing these samples."
        )
        return
    measure = metadata.get("measure", "intensity")
    normalization = metadata.get("normalization", "none")
    st.caption(
        f"Area measure: {MEASURES.get(measure, measure)}. "
        f"Normalization: {NORMALIZATIONS.get(normalization, normalization)}. "
        "Empty cells are missing measurements. Raw areas are retained for download."
    )
    for warning in metadata.get("warnings", []):
        st.warning(str(warning))
    labels = _sample_labels(samples)
    st.markdown("**RNA abundance matrix**")
    st.dataframe(matrix.rename(columns=labels), use_container_width=True)
    with st.expander("Raw extracted areas"):
        st.dataframe(raw.rename(columns=labels), use_container_width=True)
    _download_frame(
        "Download abundance matrix (TSV)", matrix, "abundances.tsv", index_label="label"
    )
    _download_frame(
        "Download raw areas (TSV)", raw, "raw-abundances.tsv", index_label="label"
    )
    _download_frame(
        "Download sample design (TSV)",
        design,
        "sample-design.tsv",
        index_label="sample_id",
    )
    st.download_button(
        "Download comparison settings (JSON)",
        json.dumps(metadata, indent=2).encode(),
        "comparison.json",
        "application/json",
    )
    if matrix.empty:
        st.info("No quantified RNA targets are available for comparison.")
        return
    target = st.selectbox(
        "RNA oligonucleotide", list(matrix.index), key=f"{prefix}-target"
    )
    values = matrix.loc[target]
    figure = go.Figure(
        go.Bar(
            x=[labels.get(column, column) for column in matrix.columns],
            y=[None if pd.isna(value) else float(value) for value in values],
            name=str(target),
        )
    )
    figure.update_layout(
        title=str(target),
        xaxis_title="Acquisition",
        yaxis_title=(
            "Relative RNA abundance"
            if normalization == "total"
            else MEASURES.get(measure, measure)
        ),
        showlegend=False,
    )
    st.plotly_chart(figure, use_container_width=True, key=f"{prefix}-abundance")
    conditions = sorted(
        {
            str(value).strip()
            for value in design["condition"].dropna()
            if str(value).strip()
        }
    )
    if len(conditions) < 2:
        st.info(
            "Assign at least two conditions to compare biological replicate abundances."
        )
        return
    st.markdown("**Condition comparison**")
    numerator = st.selectbox(
        "Numerator condition",
        conditions,
        index=len(conditions) - 1,
        key=f"{prefix}-numerator",
    )
    denominator = st.selectbox(
        "Denominator condition",
        [value for value in conditions if value != numerator],
        key=f"{prefix}-denominator",
    )
    st.caption(
        "Technical repeats are averaged within each biological replicate first. "
        "Means, standard deviations and counts describe biological replicates. "
        "Ratios use numerator / denominator; log2 fold change uses that ratio."
    )
    try:
        contrasts = compare_conditions(matrix, design, numerator, denominator)
    except (ValueError, KeyError) as error:
        st.error(f"Could not compare these conditions: {error}")
        return
    st.dataframe(contrasts, use_container_width=True)
    exported = contrasts.copy()
    for field, value in reversed(
        [
            ("numerator_condition", numerator),
            ("denominator_condition", denominator),
            ("measure", measure),
            ("normalization", normalization),
        ]
    ):
        exported.insert(0, field, value)
    _download_frame(
        "Download condition comparison (TSV)",
        exported,
        "condition-comparison.tsv",
        index_label="label",
    )


def render_experiment_results(workflow, render_sample):
    """Render saved experiment outputs; return False for a legacy single run."""
    from src.naseweis_experiment import load_manifest, sample_results

    root = Path(workflow.results_dir)
    try:
        manifest = load_manifest(root)
    except (OSError, ValueError) as error:
        st.error(f"Could not read this experiment: {error}")
        return True
    if manifest is None:
        return False
    samples = manifest["samples"]
    prefix = (
        "experiment-" + hashlib.sha256(str(root.resolve()).encode()).hexdigest()[:12]
    )
    summary, quantification, details = st.tabs(
        ["Experiment summary", "Quantification", "Sample details"]
    )
    with summary:
        finished = sum(sample["status"] == "finished" for sample in samples)
        st.metric("Completed acquisitions", f"{finished} / {len(samples)}")
        columns = [
            "id",
            "file",
            "condition",
            "replicate",
            "fraction",
            "status",
            "error",
        ]
        table = pd.DataFrame(samples).reindex(columns=columns)
        st.dataframe(
            table.dropna(axis=1, how="all"), hide_index=True, use_container_width=True
        )
        if finished < len(samples):
            st.warning(
                "This experiment has incomplete acquisitions. Review their statuses and available sample results below."
            )
    with quantification:
        if not manifest.get("lfq_enabled", False):
            st.info(
                "This is an identification-only experiment. Label-free quantification was not enabled for this run."
            )
        elif (
            manifest.get("status") != "finished"
            or manifest.get("comparison_status") != "finished"
            or not all(sample["status"] == "finished" for sample in samples)
        ):
            st.warning(
                "Comparison is unavailable for this incomplete experiment. A complete cohort is required; available acquisition results remain accessible under Sample details."
            )
        else:
            _render_comparison(root, samples, prefix)
    with details:
        if not samples:
            st.info("This experiment has no acquisitions.")
        else:
            by_id = {sample["id"]: sample for sample in samples}
            labels = _sample_labels(samples)
            selected = st.selectbox(
                "Acquisition details",
                list(by_id),
                format_func=lambda sample_id: (
                    f"{labels[sample_id]} · {by_id[sample_id]['status']}"
                ),
                key=f"{prefix}-sample",
            )
            sample = by_id[selected]
            try:
                available = sample_results(root, selected).is_dir()
            except (OSError, ValueError) as error:
                st.error(f"Could not open this acquisition's results: {error}")
            else:
                if available:
                    render_sample(sample)
                else:
                    st.info("No sample result files are available yet.")
    return True
