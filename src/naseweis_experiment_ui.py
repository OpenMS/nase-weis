"""Configure RNA experiments and inspect their saved, cross-sample results."""

import hashlib
import json
import math
from pathlib import Path

import pandas as pd
import streamlit as st

from src.naseweis_insight import show_plot, show_table
from src.naseweis_plots import abundance_bars_plot
from src.naseweis_result_plots_ui import (
    render_acquisition_scatter,
    render_experiment_coverage,
    render_lfq_scatter,
)
from src.naseweis_tools import enabled
from src.workflow.ExperimentalDesign import (
    default_design,
    design_path,
    load_design,
    save_design,
)

EXAMPLE_MZML = ["example-data/mzML/NME1.mzML"]
EXAMPLE_FASTA = ["example-data/fasta/NME1_Target.fasta"]
_EXAMPLE_NAMES = {Path(f).name for f in EXAMPLE_MZML}

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
    """Pick the acquisitions, then describe their samples (both parts in one place)."""
    selected = select_acquisitions(workflow)
    if selected is None:
        return
    configure_sample_design(workflow, selected)


def select_acquisitions(workflow):
    """The acquisitions of the next experiment, or None while the saved selection is invalid.

    Migrates legacy single-file selections before rendering the multi-select.
    """
    key = "mzML-files"
    widget = f"{workflow.parameter_manager.param_prefix}{key}"
    value = st.session_state.get(widget, workflow.params.get(key, []))
    selected = [value] if isinstance(value, str) and value else value
    if selected is None or selected == "":
        selected = []
    if not isinstance(selected, list) or any(
        not isinstance(path, str) or not path.strip() for path in selected
    ):
        st.error(
            "The saved acquisition selection is invalid. Select the input files again."
        )
        if st.button("Reset acquisition selection", key=f"{widget}-reset"):
            params = workflow.parameter_manager.read_parameters_strict()
            params[key] = []
            workflow.parameter_manager.write_parameters(params)
            st.session_state[widget] = []
            st.session_state.pop(f"_pending-{widget}", None)
            st.rerun()
        return None
    available = set(workflow.ui._input_file_options(key))
    unavailable = [path for path in selected if path not in available]
    # The example leaves once files of the user's own are added: no warning for it
    missing = [Path(path).name for path in unavailable if Path(path).name not in _EXAMPLE_NAMES]
    if missing:
        st.warning(
            "Previously selected acquisitions are no longer available in this workspace: "
            + ", ".join(missing)
            + ". Add them again with Browse or select other inputs."
        )
    if unavailable:
        selected = [path for path in selected if path in available]
    if value != selected or workflow.params.get(key) != selected:
        st.session_state[widget] = selected
        if not workflow.parameter_manager.save_parameters():
            return None
    workflow.params[key] = selected
    workflow.ui.params[key] = selected
    workflow.ui.select_input_file(
        key,
        name="MS/MS acquisitions",
        multiple=True,
        reactive=True,
        upload={"file_types": ["mzML", "raw", "d"], "fallback": EXAMPLE_MZML},
    )
    selected = st.session_state.get(widget, workflow.ui.params.get(key, []))
    workflow.params[key] = list(selected)
    workflow.ui.params[key] = list(selected)
    return list(selected)


def configure_sample_design(workflow, selected=None):
    """Conditions and replicates of the selected acquisitions, with the LFQ checks."""
    key = "mzML-files"
    if selected is None:
        widget = f"{workflow.parameter_manager.param_prefix}{key}"
        value = st.session_state.get(widget, workflow.params.get(key, []))
        value = [value] if isinstance(value, str) and value else (value or [])
        if not isinstance(value, list) or any(
            not isinstance(path, str) or not path.strip() for path in value
        ):
            st.info("Reset the acquisition selection under **Your analysis** to describe the samples.")
            return
        selected = list(value)
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


def _download_frame(label, frame, filename, *, index_label, key=None):
    st.download_button(
        label,
        frame.to_csv(sep="\t", index_label=index_label, na_rep="").encode(),
        filename,
        "text/tab-separated-values",
        key=key,
    )


def _conditions(design):
    return sorted({str(v).strip() for v in design["condition"].dropna() if str(v).strip()})


def _choice_keys(prefix):
    return f"{prefix}-cmp-numerator", f"{prefix}-cmp-denominator", f"{prefix}-cmp-normalization"


def comparison_choices(prefix, conditions, saved_normalization):
    """(numerator, denominator, normalization) picked in Edit comparison, or
    the defaults: the last condition over the first other one, and the
    normalization the run was compared with."""
    num_key, den_key, norm_key = _choice_keys(prefix)
    numerator = st.session_state.get(num_key)
    if numerator not in conditions:
        numerator = conditions[-1] if conditions else None
    others = [c for c in conditions if c != numerator]
    denominator = st.session_state.get(den_key)
    if denominator not in others:
        denominator = others[0] if others else None
    normalization = st.session_state.get(norm_key)
    if normalization not in NORMALIZATIONS:
        normalization = saved_normalization if saved_normalization in NORMALIZATIONS else "total"
    return numerator, denominator, normalization


def reset_comparison_choices(prefix):
    for name in _choice_keys(prefix):
        st.session_state.pop(name, None)


def comparison_steps(context):
    """The steps of Edit comparison for an experiment with a finished comparison."""
    from src.naseweis_comparison import load_comparison

    prefix = context["prefix"]
    num_key, den_key, norm_key = _choice_keys(prefix)
    try:
        comparison = load_comparison(context["root"])
    except (OSError, ValueError, KeyError):
        return []
    conditions = _conditions(comparison["design"])
    saved = comparison["metadata"].get("normalization", "none")

    def conditions_step():
        numerator, denominator, _ = comparison_choices(prefix, conditions, saved)
        if len(conditions) < 2:
            st.info("Assign at least two conditions in the experimental design to compare them.")
            return
        st.caption("Ratios and log2 fold changes are numerator / denominator.")
        picked = st.selectbox("Numerator condition", conditions, index=conditions.index(numerator),
                              key=f"{num_key}-widget")
        st.session_state[num_key] = picked
        others = [c for c in conditions if c != picked]
        if denominator not in others:
            denominator = others[0]
        st.session_state[den_key] = st.selectbox(
            "Denominator condition", others, index=others.index(denominator), key=f"{den_key}-widget")

    def normalization_step():
        _, _, normalization = comparison_choices(prefix, conditions, saved)
        st.caption("Total area divides each oligonucleotide's area by the sum of the quantified "
                   "areas in its sample. None keeps the extracted areas.")
        options = list(NORMALIZATIONS)
        picked = st.segmented_control("Normalization", options, default=normalization,
                                      format_func=NORMALIZATIONS.get, key=f"{norm_key}-widget")
        st.session_state[norm_key] = picked or normalization

    return [("1. Conditions", conditions_step), ("2. Normalization", normalization_step)]


def _render_comparison(root, samples, prefix, controls=True):
    """The condition comparison of an experiment. ``controls`` draws the
    condition pickers here; without them the choices of Edit comparison
    (``comparison_choices``) are used, with the normalization picked there."""
    from src.naseweis_comparison import compare_conditions, load_comparison, normalized

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
    conditions = _conditions(design)
    if not controls:
        numerator, denominator, picked = comparison_choices(prefix, conditions, normalization)
        if picked != normalization:
            matrix, normalization = normalized(raw, picked), picked
    st.caption(
        f"Area measure: {MEASURES.get(measure, measure)}. "
        f"Normalization: {NORMALIZATIONS.get(normalization, normalization)}. "
        "Empty cells are missing measurements. Raw areas are retained for download."
    )
    for warning in metadata.get("warnings", []):
        st.warning(str(warning))
    labels = _sample_labels(samples)
    with st.expander("Abundance tables and downloads"):
        st.markdown("**RNA abundance matrix**")
        show_table(_labelled(matrix.rename(columns=labels)), f"{prefix}-abundance-table",
                   "naseweis-abundances")
        st.markdown("**Raw extracted areas**")
        show_table(_labelled(raw.rename(columns=labels)), f"{prefix}-raw-table", "naseweis-raw-areas")
        _download_frame(
            "Download abundance matrix (TSV)", matrix, "abundances.tsv", key=f"{prefix}-abundances.tsv", index_label="label"
        )
        _download_frame(
            "Download raw areas (TSV)", raw, "raw-abundances.tsv", key=f"{prefix}-raw-abundances.tsv", index_label="label"
        )
        _download_frame(
            "Download sample design (TSV)",
            design,
            "sample-design.tsv",
            key=f"{prefix}-sample-design.tsv",
            index_label="sample_id",
        )
        st.download_button(
            "Download comparison settings (JSON)",
            json.dumps(metadata, indent=2).encode(),
            "comparison.json",
            "application/json",
            key=f"{prefix}-comparison.json",
        )
    if matrix.empty:
        st.info("No quantified RNA targets are available for comparison.")
        return
    with st.expander("Individual oligonucleotide abundances"):
        target = st.selectbox(
            "RNA oligonucleotide", list(matrix.index), key=f"{prefix}-target"
        )
        values = matrix.loc[target]
        show_plot(
            abundance_bars_plot(
                [labels.get(column, column) for column in matrix.columns],
                list(values),
                str(target),
                "Relative RNA abundance" if normalization == "total" else MEASURES.get(measure, measure),
            ),
            f"{prefix}-abundance",
            "naseweis-abundance",
        )
    if len(conditions) < 2:
        st.info(
            "Assign at least two conditions to compare biological replicate abundances."
        )
        render_acquisition_scatter(
            matrix, design, labels, prefix, measure=measure, normalization=normalization
        )
        return
    st.markdown("**Condition comparison**")
    if controls:
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
    else:
        st.caption(f"**{numerator}** over **{denominator}**. Change the conditions and the "
                   "normalization with **Edit comparison** above.")
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
    render_lfq_scatter(
        contrasts, numerator, denominator, prefix,
        measure=measure, normalization=normalization,
    )
    leading = ["log2_fold_change", "ratio"]
    show_table(
        _labelled(contrasts[leading + [column for column in contrasts if column not in leading]]),
        f"{prefix}-contrast-table", "naseweis-condition-comparison",
        titles={
            "label": "RNA oligonucleotide",
            "log2_fold_change": "Log2 fold change",
            "ratio": "Ratio",
            "numerator_mean": f"{numerator} mean",
            "numerator_std": f"{numerator} SD",
            "numerator_n": f"{numerator} replicates",
            "denominator_mean": f"{denominator} mean",
            "denominator_std": f"{denominator} SD",
            "denominator_n": f"{denominator} replicates",
        },
    )
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
        key=f"{prefix}-condition-comparison.tsv",
        index_label="label",
    )


def experiment_context(workflow):
    """Read a run's experiment; None for a legacy single run.

    Raises OSError or ValueError when the manifest cannot be read.
    """
    from src.naseweis_experiment import load_manifest

    root = Path(workflow.results_dir)
    manifest = load_manifest(root)
    if manifest is None:
        return None
    partial_run = getattr(workflow, "_result_run_status", None) in {
        "failed", "cancelled", "interrupted"
    }

    def displayed_status(sample):
        # A stopped worker cannot finalize its manifest. Retain the saved file,
        # but do not present its unfinished steps as currently running.
        if partial_run and sample["status"] in {"pending", "running"}:
            return "not completed"
        return sample["status"]

    return {
        "root": root,
        "manifest": manifest,
        "samples": manifest["samples"],
        "partial_run": partial_run,
        "displayed_status": displayed_status,
        "prefix": (
            "experiment-"
            + hashlib.sha256(str(root.resolve()).encode()).hexdigest()[:12]
        ),
    }


def experiment_summary_section(context):
    samples = context["samples"]
    finished = sum(sample["status"] == "finished" for sample in samples)
    summary_tiles(context)
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
    table["status"] = [context["displayed_status"](sample) for sample in samples]
    show_table(table.dropna(axis=1, how="all"), f"{context['prefix']}-samples", "naseweis-samples",
               title="Acquisitions", titles={c: c.capitalize() for c in columns})
    if finished < len(samples):
        st.warning(
            "This experiment has incomplete acquisitions. Review their statuses and available sample results below."
        )


def experiment_figures(context):
    """What the experiment found: completed acquisitions, conditions, and the
    spectrum matches, oligonucleotides and RNAs pooled over its acquisitions."""
    from src.naseweis_result_plots_ui import coverage_accessions, experiment_coverage_tables

    samples = context["samples"]
    groups = experiment_coverage_tables(context["root"], samples)
    matches, oligos = 0, set()
    for _sample, tables in groups:
        osm = tables.get("OSM")
        if osm is not None and not osm.empty:
            matches += len(osm)
            if "sequence" in osm:
                oligos.update(osm["sequence"].dropna().astype(str))
    conditions = {str(s.get("condition") or "") for s in samples} - {""}
    return {
        "finished": sum(sample["status"] == "finished" for sample in samples),
        "acquisitions": len(samples),
        "conditions": len(conditions),
        "matches": matches,
        "oligonucleotides": len(oligos),
        "rnas": len(coverage_accessions(groups)),
    }


def summary_tiles(context):
    """The experiment's figures as tiles, each with what it counts."""
    figures = experiment_figures(context)
    tiles = [
        ("Completed acquisitions", f"{figures['finished']} / {figures['acquisitions']}",
         "Acquisitions whose search finished, of all acquisitions of the experiment."),
        ("Conditions", f"{figures['conditions']:,}", "Conditions named in the experimental design."),
        ("Spectrum matches", f"{figures['matches']:,}",
         "Spectra assigned to an oligonucleotide at the FDR cutoff, over all acquisitions."),
        ("Oligonucleotides", f"{figures['oligonucleotides']:,}",
         "Distinct oligonucleotide sequences, modifications included."),
        ("RNAs", f"{figures['rnas']:,}", "RNAs of the database the oligonucleotides map to."),
    ]
    for col, (label, value, help_text) in zip(st.columns(len(tiles)), tiles):
        col.metric(label, value, help=help_text)


def _labelled(frame):
    """``frame`` with its index (the oligonucleotides) as the column ``label``."""
    shown = frame.copy()
    shown.index = shown.index.astype(str)
    shown.index.name = "label"
    return shown.reset_index()


def experiment_coverage_section(context):
    render_experiment_coverage(context["root"], context["samples"], context["prefix"])


def comparison_ready(context):
    """True when the experiment holds a finished comparison of all its samples."""
    manifest, samples = context["manifest"], context["samples"]
    return bool(
        manifest.get("lfq_enabled", False)
        and not context["partial_run"]
        and manifest.get("status") == "finished"
        and manifest.get("comparison_status") == "finished"
        and all(sample["status"] == "finished" for sample in samples)
    )


def experiment_quantification_section(context, controls=True):
    manifest = context["manifest"]
    if not manifest.get("lfq_enabled", False):
        st.info(
            "This is an identification-only experiment. Label-free quantification was not enabled for this run."
        )
    elif not comparison_ready(context):
        st.warning(
            "Comparison is unavailable for this incomplete experiment. A complete cohort is required; available acquisition results remain accessible under Sample details."
        )
    else:
        _render_comparison(context["root"], context["samples"], context["prefix"], controls=controls)


def sample_details_section(context, render_sample):
    from src.naseweis_experiment import sample_results

    samples, root = context["samples"], context["root"]
    displayed_status = context["displayed_status"]
    if not samples:
        st.info("This experiment has no acquisitions.")
        return
    by_id = {sample["id"]: sample for sample in samples}
    labels = _sample_labels(samples)
    selected = st.selectbox(
        "Acquisition details",
        list(by_id),
        format_func=lambda sample_id: (
            f"{labels[sample_id]} · {displayed_status(by_id[sample_id])}"
        ),
        key=f"{context['prefix']}-sample",
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


def render_experiment_results(workflow, render_sample):
    """Render saved experiment outputs; return False for a legacy single run."""
    try:
        context = experiment_context(workflow)
    except (OSError, ValueError) as error:
        st.error(f"Could not read this experiment: {error}")
        return True
    if context is None:
        return False
    for title, draw in (
        ("Experiment summary", lambda: experiment_summary_section(context)),
        ("Sequence coverage", lambda: experiment_coverage_section(context)),
        ("Quantification", lambda: experiment_quantification_section(context)),
        ("Sample details", lambda: sample_details_section(context, render_sample)),
    ):
        with st.container(border=True):
            st.subheader(title)
            draw()
    return True
