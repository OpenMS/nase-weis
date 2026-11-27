"""What Explore results offers: the components of each kind of dataset, its views.

Component names are stored in saved views: never rename one.

One view per former result tab (Summary, Matches, Spectra, Coverage,
Quantification, Export; for an experiment Summary, Coverage, Quantification,
Samples), so every result section NASE-weis has shown stays in a view.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from src.layout_manager import experiment_state
from src.viewer import panels, sources, teaching
from src.workflow.guide import FEATURES, IDENTIFICATION, MS_DATA, QUANTIFICATION

RESULTS = "Results"
QUANT = "quantification"

ANALYSES = {
    QUANT: {
        "label": "quantification",
        "help": "MS1 feature intensities of the identified oligonucleotides, switched on "
                "under Quantification when the run is started.",
    }
}

ACQUISITION_COMPONENTS = [
    {"name": "summary", "label": "Summary", "group": RESULTS},
    {"name": "matches", "label": "Matches", "group": RESULTS},
    {"name": "spectrum_viewer", "label": "Spectrum viewer", "group": RESULTS},
    {"name": "coverage", "label": "Sequence coverage", "group": RESULTS},
    {"name": "rna_table", "label": "Identified RNAs", "group": IDENTIFICATION},
    {"name": "quantification", "label": "MS1 quantification", "group": RESULTS, "needs": QUANT},
    {"name": "downloads", "label": "Downloads", "group": RESULTS},
    {"name": "file_scans", "label": "Scans", "group": MS_DATA},
    {"name": "file_spectrum", "label": "Spectrum", "group": MS_DATA},
    {"name": "file_chromatogram", "label": "Chromatogram", "group": MS_DATA},
    {"name": "file_peak_map", "label": "Peak map", "group": MS_DATA},
    {"name": "match_table", "label": "Identified spectra", "group": IDENTIFICATION},
    {"name": "score_distribution", "label": "Score distribution", "group": IDENTIFICATION},
    {"name": "precursor_error", "label": "Precursor mass error", "group": IDENTIFICATION},
]

EXPERIMENT_COMPONENTS = [
    {"name": "experiment_summary", "label": "Experiment summary", "group": RESULTS},
    {"name": "experiment_coverage", "label": "Sequence coverage", "group": RESULTS},
    {"name": "rna_table", "label": "Identified RNAs", "group": IDENTIFICATION},
    {"name": "experiment_quantification", "label": "Condition comparison", "group": RESULTS,
     "needs": QUANT},
    {"name": "sample_details", "label": "Sample details", "group": RESULTS},
    {"name": "abundance_table", "label": "Abundances", "group": QUANTIFICATION, "needs": QUANT},
    {"name": "abundance_profile", "label": "Abundance profile", "group": QUANTIFICATION,
     "needs": QUANT},
    {"name": "zscore_heatmap", "label": "Abundance z-scores", "group": QUANTIFICATION,
     "needs": QUANT},
]

ACQUISITION_PRESETS = [
    {"name": "summary", "label": "Summary",
     "layout": [["summary"], ["score_distribution", "precursor_error"]]},
    {"name": "matches", "label": "Matches",
     "layout": [["match_table", "file_spectrum"], ["spectrum_viewer"], ["matches"]]},
    {"name": "spectra", "label": "Spectra",
     "layout": [["file_scans", "file_spectrum"], ["file_chromatogram", "file_peak_map"],
                ["match_table"]]},
    {"name": "coverage", "label": "Identification", "layout": [["rna_table"], ["coverage"]]},
    {"name": "quantification", "label": "Quantification", "layout": [["quantification"]]},
    {"name": "export", "label": "Export", "layout": [["downloads"]]},
    {"name": "identification", "label": "Search quality",
     "layout": [["summary"], ["score_distribution", "precursor_error"], ["match_table"]]},
]

EXPERIMENT_PRESETS = [
    {"name": "summary", "label": "Summary", "layout": [["experiment_summary"]]},
    {"name": "coverage", "label": "Identification",
     "layout": [["rna_table"], ["experiment_coverage"]]},
    {"name": "quantification", "label": "Quantification",
     "layout": [["experiment_quantification"], ["abundance_table", "abundance_profile"]]},
    {"name": "abundances", "label": "Abundances",
     "layout": [["abundance_table", "abundance_profile"], ["zscore_heatmap"]]},
    {"name": "samples", "label": "Samples", "layout": [["sample_details"]]},
]

for _kind, _presets in (("", ACQUISITION_PRESETS), ("experiment", EXPERIMENT_PRESETS)):
    for _preset in _presets:
        _preset["help"] = teaching.view(_preset["name"], _kind)

ACQUISITION_RENDERERS = {c["name"]: getattr(panels, c["name"]) for c in ACQUISITION_COMPONENTS}
EXPERIMENT_RENDERERS = {c["name"]: getattr(panels, c["name"]) for c in EXPERIMENT_COMPONENTS}


def states(kind: str):
    def state(ctx):
        return experiment_state(f"naseweis-viewer-{kind}", ctx.index, ctx.dataset)
    return state


def renderers(wf, table: dict):
    def build(ctx):
        src = sources.resolve(wf, ctx.dataset)
        ctx.state = ctx.kind.states(ctx)

        def gone() -> None:
            st.info("This run is no longer available.")

        def bind(render):
            return (lambda: render(ctx, src)) if src is not None else gone
        return {name: bind(render) for name, render in table.items()}
    return build


def available(wf):
    def find(ctx) -> set[str]:
        src = sources.resolve(wf, ctx.dataset)
        if src is None:
            return set()
        if src.sample is None and sources.kind_of(wf, ctx.dataset) == sources.EXPERIMENT:
            # Started with quantification: its section then says why a comparison
            # is missing (an incomplete cohort) instead of offering to add it.
            return {QUANT} if _experiment_quant(src)[0] else set()
        return {QUANT} if src.lfq_table() is not None or src.lfq_requested() else set()
    return find


def _experiment_quant(src) -> tuple[bool, bool]:
    """(started with quantification, holds a finished comparison) of an experiment."""
    from src.naseweis_experiment_ui import comparison_ready, experiment_context

    try:
        context = experiment_context(src.wf)
    except (OSError, ValueError):
        return False, False
    if context is None:
        return False, False
    return bool(context["manifest"].get("lfq_enabled", False)), comparison_ready(context)


def unavailable_options(wf):
    def options(missing, view, ctx):
        run_id, _ = sources.split(ctx.dataset)

        def start() -> None:
            wf.follow_up(
                run_id, changes={"enable-lfq": True}, page=wf.configure_page,
                note="Quantification added to the settings of this run.",
            )
        return [{
            "label": "Search again with quantification",
            "help": "Opens Describe experiment with this run's acquisitions and settings and "
                    "quantification switched on; the new run sits next to this one.",
            "icon": "⚖️",
            "on_click": start,
        }]
    return options


def views_path(wf, kind: str) -> Path:
    return Path(wf.workflow_dir).parent / "layouts" / f"naseweis-{kind}.json"


def comparison_button(wf):
    """**Edit comparison** beside the selector of an experiment with a finished comparison."""
    def draw(ctx) -> None:
        src = sources.resolve(wf, ctx.dataset)
        if src is None or src.sample is not None:
            return
        from src.naseweis_experiment_ui import (
            comparison_ready, comparison_steps, experiment_context, reset_comparison_choices)
        from src.layout_manager.comparison import edit_comparison

        try:
            context = experiment_context(src.wf)
        except (OSError, ValueError):
            return
        if context is None or not comparison_ready(context):
            return
        prefix = context["prefix"]
        edit_comparison(
            f"{prefix}-comparison", comparison_steps(context),
            intro="Pick the two conditions to compare and how abundances are normalized. "
                  "Every quantification component of this experiment follows.",
            on_reset=lambda: reset_comparison_choices(prefix),
            button_help="Conditions and normalization of the comparison shown below.",
        )
    return draw
