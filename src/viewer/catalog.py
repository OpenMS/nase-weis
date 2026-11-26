"""What Explore results offers: the components of each kind of dataset, its views.

Component names are stored in saved views: never rename one.

The first view of each kind, **All results**, is the page as NASE-weis always
showed it (every result section in its tabs), so nothing a run showed before
is lost; the other views arrange the same results and the OpenMS-Insight
components side by side.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from src.layout_manager import experiment_state
from src.viewer import panels, sources
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
    {"name": "all_results", "label": "All results", "group": RESULTS},
    {"name": "summary", "label": "Summary", "group": RESULTS},
    {"name": "matches", "label": "Matches", "group": RESULTS},
    {"name": "spectrum_viewer", "label": "Spectrum viewer", "group": RESULTS},
    {"name": "coverage", "label": "Sequence coverage", "group": RESULTS},
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
    {"name": "all_results", "label": "All results", "group": RESULTS},
    {"name": "experiment_summary", "label": "Experiment summary", "group": RESULTS},
    {"name": "experiment_coverage", "label": "Sequence coverage", "group": RESULTS},
    {"name": "experiment_quantification", "label": "Condition comparison", "group": RESULTS,
     "needs": QUANT},
    {"name": "sample_details", "label": "Sample details", "group": RESULTS},
    {"name": "abundance_table", "label": "Abundances", "group": QUANTIFICATION, "needs": QUANT},
    {"name": "abundance_profile", "label": "Abundance profile", "group": QUANTIFICATION,
     "needs": QUANT},
    {"name": "zscore_heatmap", "label": "Abundance z-scores", "group": QUANTIFICATION,
     "needs": QUANT},
]

ALL = [["all_results"]]

ACQUISITION_PRESETS = [
    {"name": "all", "label": "All results", "layout": ALL,
     "help": "Every result of the acquisition in tabs: summary, matches, spectra, "
             "coverage, quantification and downloads."},
    {"name": "spectra", "label": "Spectra",
     "layout": [["file_scans", "file_spectrum"], ["file_chromatogram", "file_peak_map"],
                ["match_table"]],
     "help": "The acquisition's scans, the picked spectrum, its chromatograms and peak map, "
             "with the identified spectra linked to them."},
    {"name": "identification", "label": "Identification quality",
     "layout": [["summary"], ["score_distribution", "precursor_error"], ["match_table"]],
     "help": "How many oligonucleotides were found, how their scores and precursor mass "
             "errors are distributed, and the identified spectra."},
]

EXPERIMENT_PRESETS = [
    {"name": "all", "label": "All results", "layout": ALL,
     "help": "Every result of the experiment in tabs: summary, coverage, condition "
             "comparison and the details of each sample."},
    {"name": "abundances", "label": "Abundances",
     "layout": [["abundance_table", "abundance_profile"], ["zscore_heatmap"]],
     "help": "The abundance of each oligonucleotide in every sample, the picked one by "
             "condition, and the most variable ones as z-scores."},
]

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
            from src.naseweis_experiment_ui import comparison_ready, experiment_context

            try:
                context = experiment_context(src.wf)
            except (OSError, ValueError):
                return set()
            return {QUANT} if context is not None and comparison_ready(context) else set()
        return {QUANT} if src.lfq_table() is not None else set()
    return find


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
