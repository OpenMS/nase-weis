"""Explore results: the template's standard Viewer over NASE-weis runs.

A run of several acquisitions is shown as its **Experiment** or as one
**Single acquisition** at a time; the components of each are in
``src/viewer/catalog.py``.
"""

import streamlit as st

from src.common.common import page_setup
from src.layout_manager import DatasetKind, ViewerSpec, show_viewer
from src.viewer import catalog, panels, sources
from src.Workflow import Workflow
from src.workflow import guide
from src.workflow.run_widgets import run_import, run_label, selected_run_key, show_import_notice

page_setup()
wf = Workflow()
panels.reset()

st.title("Explore results")
guide.page_intro(
    "Choose a saved analysis to compare sample abundances and inspect each acquisition's "
    "matches, spectra and sequence coverage."
)
show_import_notice(wf.runs)

available_runs, statuses = wf.run_statuses()
datasets = {}
for run in available_runs:
    label = run_label(run)
    if statuses[run["id"]] != "finished":
        label += f" · {statuses[run['id']].capitalize()} / partial results"
    datasets[run["id"]] = label

manifests = {run_id: sources.manifest(wf, run_id) for run_id in datasets}
sample_labels = {}
for run_id, found in manifests.items():
    for sample in (found or {}).get("samples", []):
        text = f"{sample.get('file', sample['id'])} · {sample.get('condition') or 'no condition'}"
        if sample.get("replicate") not in (None, ""):
            text += f" · replicate {sample['replicate']}"
        sample_labels[sources.member_id(run_id, sample["id"])] = text


def members(run_id: str) -> list[str]:
    return [sources.member_id(run_id, s["id"]) for s in (manifests.get(run_id) or {}).get("samples", [])]


def empty() -> None:
    st.info(
        "No finished runs yet. Start the workflow on the **Run** page, "
        "or open a run someone shared as a ZIP."
    )
    run_import(wf.runs)


def partial_warning(ctx) -> None:
    run_id, _ = sources.split(ctx.experiment)
    status = statuses.get(run_id, "finished")
    if status != "finished":
        state = "failed" if status == "failed" else f"was {status}"
        st.warning(
            f"This run {state}. Available outputs may be incomplete; inspect this run's "
            "saved log on Export before interpreting them."
        )


def kind(name: str, components, presets, table, buttons=None) -> DatasetKind:
    return DatasetKind(
        name=name,
        components=components,
        presets=presets,
        default_layout=presets[0]["layout"],
        views_path=catalog.views_path(wf, name),
        renderers=catalog.renderers(wf, table),
        available=catalog.available(wf),
        analyses=catalog.ANALYSES,
        unavailable_options=catalog.unavailable_options(wf),
        states=catalog.states(name),
        before=partial_warning,
        buttons=buttons,
        dataset_label="run",
        experiment_label="Run",
    )


shared = selected_run_key(wf.runs)
# The run another page of the workflow picked (the Run page's new run, say) is
# shown once; a pick made here only follows to the other pages
if st.session_state.get(shared) in datasets and (
        st.session_state.get(shared) != st.session_state.get("_naseweis_remembered")):
    st.session_state["_naseweis_pick"] = st.session_state[shared]
if st.session_state.get("_naseweis_pick") not in datasets:
    st.session_state["_naseweis_pick"] = next(iter(datasets), None)


def remember(ctx) -> None:
    st.session_state[shared] = st.session_state["_naseweis_remembered"] = ctx.experiment


show_viewer(ViewerSpec(
    datasets=datasets,
    kinds=[
        kind(sources.ACQUISITION, catalog.ACQUISITION_COMPONENTS, catalog.ACQUISITION_PRESETS,
             catalog.ACQUISITION_RENDERERS),
        kind(sources.EXPERIMENT, catalog.EXPERIMENT_COMPONENTS, catalog.EXPERIMENT_PRESETS,
             catalog.EXPERIMENT_RENDERERS, buttons=catalog.comparison_button(wf)),
    ],
    key="naseweis-viewer",
    kind_of=lambda dataset: sources.kind_of(wf, dataset),
    members=members,
    member_mode="switch",
    member_format=lambda dataset: sample_labels.get(dataset, dataset),
    whole=lambda run_id: manifests.get(run_id) is not None,
    member_choice="Single acquisition",
    whole_choice="Experiment",
    member_select_label="Acquisition",
    whole_kind=sources.EXPERIMENT,
    label="Run",
    empty=empty,
    pick=st.session_state["_naseweis_pick"],
    on_pick=remember,
))

guide.next_steps(
    wf,
    (guide.EXPORT, "Take the result files and the processing report of a run with you."),
    (guide.DESCRIBE, "Search more acquisitions, or the same ones with other settings."),
)
