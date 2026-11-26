"""What Explore results offers: the components of each kind of dataset, its views.

Component names are stored in saved views: never rename one.

The first view of each kind, **All results**, is one dashboard: every result
section NASE-weis has shown and the OpenMS-Insight components on one page, so
nothing a run showed before is lost. One view per former tab (Summary, Matches,
Spectra, Coverage, Quantification, Export; for an experiment Summary, Coverage,
Quantification, Samples) shows one part of it.
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

# "All results" is one dashboard: every result of the run as components on
# one page. Components that need quantification are added by ``rows`` when the
# shown run has it, so the view is never locked for a run without it.
ACQUISITION_ALL = [
    ["summary"],
    ["match_table", "file_spectrum"],
    ["spectrum_viewer"],
    ["coverage"],
    ["file_chromatogram", "file_peak_map"],
    ["matches"],
    ["downloads"],
]
ACQUISITION_QUANT_ROWS = [["quantification"]]
EXPERIMENT_ALL = [
    ["experiment_summary"],
    ["experiment_coverage"],
    ["sample_details"],
]
EXPERIMENT_QUANT_ROWS = [["experiment_quantification"], ["abundance_table", "abundance_profile"],
                         ["zscore_heatmap"]]
# A view saved before the dashboard named this single component
ALL = [["all_results"]]

ACQUISITION_PRESETS = [
    {"name": "all", "label": "All results", "layout": ACQUISITION_ALL,
     "help": "Every result of the acquisition on one page: summary, matches with the "
             "picked spectrum, its annotated fragments, coverage, MS data, quantification "
             "and downloads."},
    {"name": "summary", "label": "Summary",
     "layout": [["summary"], ["score_distribution", "precursor_error"]],
     "help": "How many spectra, oligonucleotides and RNAs were found, and how their "
             "scores and precursor mass errors are distributed."},
    {"name": "matches", "label": "Matches",
     "layout": [["match_table", "file_spectrum"], ["spectrum_viewer"], ["matches"]],
     "help": "The spectrum matches; a click shows the matched spectrum and its annotated "
             "fragments."},
    {"name": "spectra", "label": "Spectra",
     "layout": [["file_scans", "file_spectrum"], ["file_chromatogram", "file_peak_map"],
                ["match_table"]],
     "help": "The acquisition's scans, the picked spectrum, its chromatograms and peak map, "
             "with the identified spectra linked to them."},
    {"name": "coverage", "label": "Coverage", "layout": [["coverage"]],
     "help": "Which residues of each RNA the identified oligonucleotides cover."},
    {"name": "quantification", "label": "Quantification", "layout": [["quantification"]],
     "help": "MS1 feature intensities of the identified oligonucleotides and their "
             "extracted ion chromatograms."},
    {"name": "export", "label": "Export", "layout": [["downloads"]],
     "help": "The result files, the search database, parameters and run log."},
    {"name": "identification", "label": "Identification quality",
     "layout": [["summary"], ["score_distribution", "precursor_error"], ["match_table"]],
     "help": "How many oligonucleotides were found, how their scores and precursor mass "
             "errors are distributed, and the identified spectra."},
]

EXPERIMENT_PRESETS = [
    {"name": "all", "label": "All results", "layout": EXPERIMENT_ALL,
     "help": "Every result of the experiment on one page: summary, coverage, the "
             "condition comparison with the abundances, and the details of each sample."},
    {"name": "summary", "label": "Summary", "layout": [["experiment_summary"]],
     "help": "The experiment's acquisitions, their conditions and status."},
    {"name": "coverage", "label": "Coverage", "layout": [["experiment_coverage"]],
     "help": "Which residues of each RNA the oligonucleotides cover, per condition."},
    {"name": "quantification", "label": "Quantification",
     "layout": [["experiment_quantification"], ["abundance_table", "abundance_profile"]],
     "help": "The condition comparison, the abundance table and the picked "
             "oligonucleotide by condition."},
    {"name": "abundances", "label": "Abundances",
     "layout": [["abundance_table", "abundance_profile"], ["zscore_heatmap"]],
     "help": "The abundance of each oligonucleotide in every sample, the picked one by "
             "condition, and the most variable ones as z-scores."},
    {"name": "samples", "label": "Samples", "layout": [["sample_details"]],
     "help": "The results of each acquisition of the experiment."},
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


def _dashboard(kind: str):
    if kind == sources.EXPERIMENT:
        return EXPERIMENT_ALL, EXPERIMENT_QUANT_ROWS, 2
    return ACQUISITION_ALL, ACQUISITION_QUANT_ROWS, len(ACQUISITION_ALL) - 1


def rows(wf, kind: str):
    """The ``rows`` hook: the dashboard gets the quantification of a run that has
    it, and the component ``all_results`` of a view saved before the dashboard
    becomes the dashboard's rows."""
    base, quant, at = _dashboard(kind)
    find = available(wf)

    def _ready(ctx) -> bool:
        src = sources.resolve(wf, ctx.dataset)
        return src is not None and _experiment_quant(src)[1]

    def dashboard(ctx) -> list[list[str]]:
        shown = [list(row) for row in base]
        if QUANT in find(ctx):
            # An incomplete cohort keeps only the section's note on why.
            added = quant if kind != sources.EXPERIMENT or _ready(ctx) else quant[:1]
            shown[at:at] = [list(row) for row in added]
        elif kind != sources.EXPERIMENT:
            # A run searched without quantification keeps the section, drawn as
            # the note that offers to add it, so the dashboard says why it is missing.
            shown[at:at] = [list(quant[0])]
        return shown

    def hook(ctx, shown, saved):
        if [list(row) for row in shown] == base:
            return dashboard(ctx)
        out = []
        for row in shown:
            if "all_results" in row:
                out.extend(dashboard(ctx))
                rest = [name for name in row if name != "all_results"]
                if rest:
                    out.append(rest)
            else:
                out.append(list(row))
        return out
    return hook
