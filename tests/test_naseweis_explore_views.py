"""Regression tests for Explore results on the standard Viewer.

Every view of both kinds is opened on both demo runs: no exception, no error,
and what NASE-weis showed before (its result tabs, downloads, the spectrum and
LFQ viewers, the condition comparison) is still there.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from src.layout_manager.layout import select_view
from src.viewer import catalog
from tests.naseweis_drawn import tables
from test_gui import _demo, _run_box, app  # noqa: F401  (the app fixture)

LFQ_RUN = "20261126-072337"
ID_RUN = "20261126-072320"
SWITCH = "naseweis-viewer_show_0_widget"
ACQUISITION_TABS = ["Matches", "Spectra", "Sequence coverage map", "Quantification", "Export"]
EXPERIMENT_TABS = ["Experiment summary", "Sequence coverage map", "Quantification", "Sample details"]


def _clean(app):
    assert not app.exception, [str(e.value)[:500] for e in app.exception]
    assert not app.error, [e.value for e in app.error]


def _tabs(app):
    # Each result section is a titled panel; the coverage map is one Insight
    # SequenceCoverage, whose own frame carries its title.
    titles = [header.value for header in app.subheader]
    if any(item.component == "SequenceCoverage" for item in getattr(app, "drawn", [])) or any(
        text.value.startswith("**Sequence coverage") for text in app.markdown
    ):
        titles.append("Sequence coverage map")
    return titles


def _has(app, sections):
    shown = _tabs(app)
    return all(section in shown for section in sections)


def _insight(app):
    # One custom component is the view bar itself; components NASE-weis draws
    # through its own plots are recorded instead of drawn.
    return len(app.get("component_instance")) - 1 + len(app.drawn)


def _tour(app, kind):
    """Open every ready-made view of ``kind`` and collect what each one drew."""
    presets = catalog.EXPERIMENT_PRESETS if kind == "experiment" else catalog.ACQUISITION_PRESETS
    seen = {"tabs": set(), "buttons": set(), "selectboxes": set(), "markdown": set(),
            "downloads": 0, "drawn": False}
    for preset in presets:
        _open_view(app, kind, preset["name"])
        _clean(app)
        seen["tabs"].update(_tabs(app))
        seen["buttons"].update(b.label for b in app.button)
        seen["selectboxes"].update(box.label for box in app.selectbox)
        seen["markdown"].update(m.value for m in app.markdown)
        seen["downloads"] += len(app.get("download_button"))
        seen["drawn"] = seen["drawn"] or bool(tables(app.drawn) or app.get("component_instance"))
    return seen


def _sections(run):
    return ACQUISITION_TABS if run == LFQ_RUN else [
        t for t in ACQUISITION_TABS if t != "Quantification"]


def _views(app, kind):
    return Path(app.session_state.workspace) / "layouts" / f"naseweis-{kind}.json"


def _open_view(app, kind, name):
    components, presets = (
        (catalog.EXPERIMENT_COMPONENTS, catalog.EXPERIMENT_PRESETS)
        if kind == "experiment" else
        (catalog.ACQUISITION_COMPONENTS, catalog.ACQUISITION_PRESETS))
    assert select_view(_views(app, kind), name, components, presets, presets[0]["layout"])
    app.run(timeout=120)


def _single_acquisition(app):
    switch = next(b for b in app.button_group if b.key == SWITCH)
    switch.set_value(["member"]).run(timeout=120)


def _show(app, run):
    _demo(app, every_section=False)
    if run != LFQ_RUN:
        _run_box(app).set_value(run).run(timeout=120)
    assert _run_box(app).value == run


# --- Round 1: the page as NASE-weis showed it, on the standard Viewer --------


def test_quantified_experiment_keeps_every_result_section(app):
    _show(app, LFQ_RUN)
    _clean(app)
    seen = _tour(app, "experiment")
    assert set(EXPERIMENT_TABS) <= seen["tabs"]
    assert {"Open spectrum viewer", "Open LFQ feature viewer", "Edit comparison"} <= seen["buttons"]
    for label in ("Acquisition details", "RNA accession", "RNA oligonucleotide"):
        assert label in seen["selectboxes"]
    assert "**Condition comparison**" in seen["markdown"]
    assert seen["drawn"]


def test_experiment_switch_offers_each_acquisition(app):
    _show(app, LFQ_RUN)
    _open_view(app, "experiment", "samples")
    acquisition = next(box for box in app.selectbox if box.label == "Acquisition details")
    assert len(acquisition.options) == 2
    _single_acquisition(app)
    _clean(app)
    picker = next(box for box in app.selectbox if box.label == "Acquisition")
    assert len(picker.options) == 2
    assert all("replicate" in o or "condition" in o for o in picker.options)
    assert set(ACQUISITION_TABS) <= _tour(app, "acquisition")["tabs"]
    picker = next(box for box in app.selectbox if box.label == "Acquisition")
    picker.set_value(picker.options[1]).run(timeout=120)
    _clean(app)
    assert set(ACQUISITION_TABS) <= _tour(app, "acquisition")["tabs"]


@pytest.mark.parametrize("run", [LFQ_RUN, ID_RUN])
def test_single_acquisition_keeps_downloads_and_spectra(app, run):
    _show(app, run)
    _single_acquisition(app)
    _clean(app)
    # A run without quantification leaves that panel out;
    # its Quantification view stays in the bar, locked, with a new search.
    seen = _tour(app, "acquisition")
    assert set(_sections(run)) <= seen["tabs"]
    assert seen["downloads"]
    if run == ID_RUN:
        assert "Quantification" not in seen["tabs"]


def test_identification_only_experiment_has_no_comparison(app):
    _show(app, ID_RUN)
    _open_view(app, "experiment", "coverage")
    _clean(app)
    assert not any(m.value == "**Condition comparison**" for m in app.markdown)
    # The RNA picked in the table (or, before a pick, in the accession box) is
    # the reference the coverage map names in its own frame.
    assert (any("NME1" in str(box.options) for box in app.selectbox if "accession" in box.label)
            or any(item.settings.get("reference_name") == "NME1" for item in app.drawn
                   if item.component == "SequenceCoverage")
            or any(text.value.startswith("**Sequence coverage") and "NME1" in text.value
                   for text in app.markdown))


# --- Round 1: every ready-made view of both kinds ---------------------------


@pytest.mark.parametrize("run", [LFQ_RUN, ID_RUN])
@pytest.mark.parametrize("view, components", [("spectra", 4), ("identification", 3), ("matches", 2)])
def test_acquisition_views(app, run, view, components):
    _show(app, run)
    _single_acquisition(app)
    _open_view(app, "acquisition", view)
    _clean(app)
    assert _insight(app) >= components
    _open_view(app, "acquisition", "coverage")
    _clean(app)
    assert "Sequence coverage map" in _tabs(app)
    assert _insight(app) >= 1  # the table of identified RNAs


def test_quantified_abundances_view(app):
    _show(app, LFQ_RUN)
    _open_view(app, "experiment", "abundances")
    _clean(app)
    assert _insight(app) == 3


def test_abundances_view_without_quantification_offers_a_new_search(app):
    _show(app, ID_RUN)
    _open_view(app, "experiment", "abundances")
    _clean(app)
    # The view is locked for a run without quantification: the page stays on
    # its first view, and the view bar offers the new search instead
    assert "Experiment summary" in _tabs(app)
    assert json.loads(_views(app, "experiment").read_text())["current"] == "summary"
    calls = []
    wf = SimpleNamespace(configure_page="describe",
                         follow_up=lambda *a, **k: calls.append((a, k)))
    options = catalog.unavailable_options(wf)({catalog.QUANT}, None,
                                              SimpleNamespace(dataset=ID_RUN))
    assert [o["label"] for o in options] == ["Search again with quantification"]
    options[0]["on_click"]()
    assert calls[0][0] == (ID_RUN,)
    assert calls[0][1]["changes"] == {"enable-lfq": True}


def test_every_component_has_a_renderer():
    for components, table in ((catalog.ACQUISITION_COMPONENTS, catalog.ACQUISITION_RENDERERS),
                              (catalog.EXPERIMENT_COMPONENTS, catalog.EXPERIMENT_RENDERERS)):
        assert [c["name"] for c in components] == list(table)
        assert all(callable(render) for render in table.values())
    for presets, components in ((catalog.ACQUISITION_PRESETS, catalog.ACQUISITION_COMPONENTS),
                                (catalog.EXPERIMENT_PRESETS, catalog.EXPERIMENT_COMPONENTS)):
        names = {c["name"] for c in components}
        assert presets[0]["name"] == "summary"
        assert all(preset.get("help") for preset in presets)
        for preset in presets:
            assert all(n in names for row in preset["layout"] for n in row)


# --- Round 2: selections stay with a run, picks follow other pages ----------


def test_run_pick_survives_view_changes(app):
    _show(app, ID_RUN)
    _single_acquisition(app)
    _open_view(app, "acquisition", "spectra")
    assert _run_box(app).value == ID_RUN
    _open_view(app, "acquisition", "coverage")
    assert _run_box(app).value == ID_RUN
    _clean(app)


def test_run_picked_here_is_shared_with_other_pages(app):
    _show(app, ID_RUN)
    shared = [k for k in app.session_state.filtered_state
              if k.startswith("_runs-") and k.endswith("-selected")]
    assert shared and all(app.session_state[k] == ID_RUN for k in shared)
