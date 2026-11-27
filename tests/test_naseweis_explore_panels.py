"""Round 3: every component of Explore results drawn on its own.

Each test draws one component of the catalog for one dataset of the demo
workspace (both runs, the experiment and each acquisition), with a spectrum
and an oligonucleotide picked where the component needs one, then without a
pick, and for a run that is gone. A component must never raise or show an
error, so no view a user can build from the palette breaks.
"""
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from src.viewer import catalog

LFQ_RUN = "20261126-072337"
ID_RUN = "20261126-072320"


def _page():
    from types import SimpleNamespace

    import streamlit as st

    from src.common.common import page_setup
    from src.layout_manager import experiment_state
    from src.viewer import catalog, panels, sources
    from src.Workflow import Workflow

    page_setup()
    wf = Workflow()
    panels.reset()
    want = st.session_state["_panel_test"]
    dataset = want["dataset"]
    if dataset.startswith("first-sample-of:"):
        run_id = dataset.split(":", 1)[1]
        dataset = sources.member_id(run_id, sources.manifest(wf, run_id)["samples"][0]["id"])
    kind = sources.kind_of(wf, dataset)
    original = sources.Source.mzml
    if want.get("no_spectra"):
        def missing(self):
            if "the acquisition was moved" not in self.notes:
                self.notes.append("the acquisition was moved")
        sources.Source.mzml = missing
    state = experiment_state("panel-test", 0, dataset)
    if not st.session_state.get("_panel_picked"):
        # Picked once, as a click would; a component may clear a pick it
        # cannot show, and setting it again on every rerun would loop.
        st.session_state["_panel_picked"] = True
        for identifier, value in want.get("picks", {}).items():
            state.set_selection(identifier, value)
    ctx = SimpleNamespace(index=want.get("index", 0), state=state, dataset=dataset)
    table = catalog.EXPERIMENT_RENDERERS if kind == sources.EXPERIMENT else catalog.ACQUISITION_RENDERERS
    renderers = catalog.renderers(wf, table)(SimpleNamespace(
        dataset=dataset, index=ctx.index, kind=SimpleNamespace(states=lambda c: state)))
    st.session_state["_panel_kind"] = kind
    try:
        for name in want["names"]:
            renderers[name]()
    finally:
        sources.Source.mzml = original
    st.session_state["_panel_available"] = sorted(catalog.available(wf)(SimpleNamespace(dataset=dataset)))


@pytest.fixture
def page(tmp_path, monkeypatch):
    monkeypatch.setenv("LOCAL_WORKSPACES_DIR", str(tmp_path))
    settings = json.loads(Path("settings.json").read_text())
    settings.update(workspaces_dir=str(tmp_path), test=False)

    def draw(dataset, names, **want):
        test = AppTest.from_function(_page, default_timeout=180)
        test.session_state.settings = settings
        test.session_state["_panel_test"] = dict(want, dataset=dataset, names=list(names))
        test.run()
        assert not test.exception, [str(e.value)[:800] for e in test.exception]
        assert not test.error, [e.value for e in test.error]
        return test
    return draw


def _names(components, skip=()):
    return [c["name"] for c in components if c["name"] not in skip]


def _sample(run):
    return f"first-sample-of:{run}"


def _insight(test):
    return len(test.get("component_instance"))


def _texts(test):
    return " ".join(str(e.value) for e in list(test.info) + list(test.markdown) + list(test.caption))


ACQUISITION = _names(catalog.ACQUISITION_COMPONENTS)
EXPERIMENT = _names(catalog.EXPERIMENT_COMPONENTS)
INSIGHT_ACQUISITION = ["file_scans", "file_spectrum", "file_chromatogram", "file_peak_map",
                       "match_table", "score_distribution", "precursor_error"]


@pytest.mark.parametrize("run", [LFQ_RUN, ID_RUN])
@pytest.mark.parametrize("name", ACQUISITION)
def test_every_acquisition_component_with_a_spectrum_picked(page, run, name):
    test = page(_sample(run), [name], picks={"scan": 0})
    assert test.session_state["_panel_kind"] == "acquisition"
    if name in INSIGHT_ACQUISITION:
        assert _insight(test) >= 1, _texts(test)


@pytest.mark.parametrize("name", ["file_spectrum"])
def test_the_spectrum_asks_for_a_pick(page, name):
    test = page(_sample(LFQ_RUN), [name])
    assert "Pick a spectrum" in _texts(test)
    assert _insight(test) == 0


@pytest.mark.parametrize("name", EXPERIMENT)
def test_every_experiment_component_of_the_quantified_run(page, name):
    test = page(LFQ_RUN, [name], picks={"oligo": 0})
    assert test.session_state["_panel_kind"] == "experiment"
    assert test.session_state["_panel_available"] == ["quantification"]
    if name in ("abundance_table", "abundance_profile", "zscore_heatmap"):
        assert _insight(test) >= 1, _texts(test)


def test_the_abundance_profile_asks_for_a_pick(page):
    test = page(LFQ_RUN, ["abundance_profile"])
    assert "Pick an oligonucleotide" in _texts(test)


@pytest.mark.parametrize("name", EXPERIMENT)
def test_experiment_components_without_quantification_explain_themselves(page, name):
    test = page(ID_RUN, [name])
    assert test.session_state["_panel_available"] == []
    if name in ("abundance_table", "abundance_profile", "zscore_heatmap"):
        assert "no finished comparison" in _texts(test)


SECTIONS_ACQUISITION = list(__import__("src.viewer.panels", fromlist=["x"]).SECTIONS_SINGLE)
SECTIONS_EXPERIMENT = list(__import__("src.viewer.panels", fromlist=["x"]).SECTIONS_EXPERIMENT)


def test_result_sections_of_both_kinds(page):
    page(LFQ_RUN, SECTIONS_EXPERIMENT)
    page(_sample(ID_RUN), SECTIONS_ACQUISITION)


def test_result_sections_are_drawn_once_per_page(page):
    test = page(_sample(LFQ_RUN), ["summary", "summary", "matches"])
    assert "Shown above in this view." in _texts(test)


def test_result_sections_are_kept_for_the_first_run(page):
    test = page(_sample(LFQ_RUN), ["summary"], index=1)
    assert "shown for the first run" in _texts(test)


@pytest.mark.parametrize("dataset", ["19990101-000000", f"{LFQ_RUN}/no-such-sample"])
def test_a_run_that_is_gone(page, dataset):
    test = page(dataset, ["summary"])
    assert "no longer available" in _texts(test)
    assert test.session_state["_panel_available"] == []


@pytest.mark.parametrize("dataset, names", [(_sample(LFQ_RUN), SECTIONS_ACQUISITION),
                                            (LFQ_RUN, SECTIONS_EXPERIMENT)])
def test_every_result_section_of_a_second_run_points_to_the_first(page, dataset, names):
    test = page(dataset, names, index=1)
    assert _texts(test).count("shown for the first run") == len(names)


@pytest.mark.parametrize("dataset, names", [(_sample(LFQ_RUN), SECTIONS_ACQUISITION),
                                            (LFQ_RUN, SECTIONS_EXPERIMENT)])
def test_every_result_section_is_shown_once_per_view(page, dataset, names):
    test = page(dataset, names + names)
    assert _texts(test).count("Shown above in this view.") == len(names)


def test_spectra_components_say_when_the_acquisition_is_missing(page):
    names = ["file_scans", "file_spectrum", "file_chromatogram", "file_peak_map", "match_table"]
    test = page(_sample(LFQ_RUN), names, picks={"scan": 0}, no_spectra=True)
    # Said once per component, not repeated for each component drawn before it
    assert _texts(test).count("the acquisition was moved") == 4
    assert "moved; the acquisition was moved" not in _texts(test)
    assert _insight(test) == 1  # the matches, without their spectra


# The identification chain: RNA -> its oligonucleotides -> their matches ->
# the fragment spectrum and fragment table of the picked match.
CHAIN = ["oligo_table", "oligo_matches", "match_spectrum", "fragment_table"]


@pytest.mark.parametrize("name, note", [
    ("oligo_table", "Pick an RNA"),
    ("oligo_matches", "Pick an oligonucleotide"),
    ("match_spectrum", "Pick a spectrum match"),
    ("fragment_table", "Pick a spectrum match"),
])
def test_the_chain_asks_for_its_pick(page, name, note):
    test = page(_sample(ID_RUN), [name])
    assert note in _texts(test)
    assert _insight(test) == 0


@pytest.mark.parametrize("run", [LFQ_RUN, ID_RUN])
def test_the_chain_follows_its_picks(page, run):
    test = page(_sample(run), CHAIN, picks={"rna": 0, "oligo_id": 0, "psm": 0})
    texts = _texts(test)
    assert "Pick an RNA" not in texts
    assert "Pick an oligonucleotide" not in texts
    assert "Pick a spectrum match" not in texts
    assert _insight(test) == 4, texts


def test_the_oligonucleotide_table_of_the_experiment(page):
    test = page(LFQ_RUN, ["oligo_table"], picks={"rna": 0})
    assert _insight(test) == 1, _texts(test)


def _chain_page():
    from types import SimpleNamespace

    import streamlit as st

    from src.common.common import page_setup
    from src.viewer import panels, sources
    from src.Workflow import Workflow

    page_setup()
    wf = Workflow()
    panels.reset()
    run_id = st.session_state["_chain_run"]
    sample = sources.manifest(wf, run_id)["samples"][0]["id"]
    src = sources.resolve(wf, sources.member_id(run_id, sample))
    rnas = panels._rnas(src)
    oligos = panels.oligo_rows(panels._rna_groups(src), rnas, panels._oligo_ids(src))
    data = panels._psms(src)[0]
    found = {}
    for row in oligos.itertuples():
        psms = data.filter(data["oligo_id"] == row.oligo_id)
        if psms.height:
            found = dict(rna=int(row.rna_idx), oligo=int(row.oligo_id), sequence=row.sequence,
                         psm=int(psms["id_idx"][0]), psm_sequence=psms["sequence"][0])
            break
    picks = {"psm": found.get("psm")}
    state = SimpleNamespace(get_selection=lambda identifier: picks.get(identifier))
    picked = panels.picked_fragments(SimpleNamespace(state=state, index=0), src)
    if picked is not None:
        _match, peaks, assignments, note = picked
        found.update(peaks=len(peaks), assigned=len(assignments), note=note or "",
                     placed=int((assignments["match_status"] == "matched").sum()) if len(assignments) else 0)
    st.session_state["_chain_found"] = found


@pytest.mark.parametrize("run", [LFQ_RUN, ID_RUN])
def test_the_chain_reaches_annotated_fragments(tmp_path, monkeypatch, run):
    """A real RNA leads to its oligonucleotides, one of them to its matches,
    and the match to NASE's fragment assignments placed on its peaks."""
    monkeypatch.setenv("LOCAL_WORKSPACES_DIR", str(tmp_path))
    settings = json.loads(Path("settings.json").read_text())
    settings.update(workspaces_dir=str(tmp_path), test=False)
    test = AppTest.from_function(_chain_page, default_timeout=180)
    test.session_state.settings = settings
    test.session_state["_chain_run"] = run
    test.run()
    assert not test.exception, [str(e.value)[:800] for e in test.exception]
    found = test.session_state["_chain_found"]
    assert found, "no RNA with a matched oligonucleotide"
    assert found["psm_sequence"] == found["sequence"]
    assert found["peaks"] > 0
    assert found["assigned"] > 0, found["note"]
    assert found["placed"] > 0
