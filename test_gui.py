"""Smoke-test every NASEWEIS page and the complete example-data setup."""

import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

# Import UI dependencies before tests that temporarily mock Streamlit at collection.
from src.common.common import page_setup  # noqa: F401


@pytest.fixture
def app(tmp_path, monkeypatch):
    # The runtime image sets this variable; always isolate each test even there.
    monkeypatch.setenv("LOCAL_WORKSPACES_DIR", str(tmp_path))
    test = AppTest.from_file("app.py")
    settings = json.loads(Path("settings.json").read_text())
    settings.update(workspaces_dir=str(tmp_path), test=False)
    test.session_state.settings = settings
    return test


@pytest.mark.parametrize(
    "page",
    [
        "content/quickstart.py",
        "content/topp_workflow_file_upload.py",
        "content/naseweis.py",
        "content/topp_workflow_execution.py",
        "content/topp_workflow_results.py",
    ],
)
def test_launch(app, page):
    app.run(timeout=30)
    if page != "content/quickstart.py":
        app.switch_page(page).run(timeout=30)
    assert not app.exception


def test_load_example_and_configure(app):
    app.run(timeout=30)
    next(
        button for button in app.button if button.label == "Load example dataset"
    ).click().run(timeout=30)
    assert not app.exception
    assert app.title[0].value == "Configure your search"
    workspace = Path(app.session_state.workspace)
    params = json.loads((workspace / "naseweis/params.json").read_text())
    assert Path(params["mzML-files"]).name == "NME1.mzML"
    assert Path(params["fasta-files"]).name == "NME1_Target.fasta"
    assert any("NME1.mzML" in option for box in app.selectbox for option in box.options)
    assert not app.exception


def _finished_run(
    app,
    name="Saved search",
    *,
    sequence="ACGU",
    q_value="0.01",
    legacy=False,
    content=None,
):
    """Seed saved output; page rendering still uses the real run manager/parser."""
    from src.workflow.RunManager import RunManager

    app.run(timeout=30)
    directory = Path(app.session_state.workspace) / "naseweis"
    directory.mkdir(exist_ok=True)
    current = directory / "params.json"
    current.write_text(json.dumps({"add-decoys": False, "FDR_cutoff": 1.0}))
    runs = RunManager(directory)
    if legacy:
        results = directory / "results"
        results.mkdir(exist_ok=True)
        run_id = None
    else:
        run_id = runs.new_run(name, current)
        runs.params_file(run_id).write_text(
            json.dumps({"add-decoys": True, "FDR_cutoff": 0.05})
        )
        results = runs.results_dir(run_id)
    (results / "mztab_results").mkdir(parents=True)
    (results / "idxml_results").mkdir(parents=True)
    (results / "idxml_results" / ("id_out" if legacy else "search.idXML")).write_text(
        "<IdXML/>"
    )
    tab_file = results / "mztab_results" / ("tab_out" if legacy else "search.mzTab")
    tab_file.write_text(
        content
        if content is not None
        else (
            "NUH\taccession\topt_sequence\nNUC\trna\tACGU\n"
            f"OLH\tsequence\taccession\tstart\tend\nOLI\t{sequence}\trna\t1\t4\n"
            "OSH\tsequence\tretention_time\tcharge\texp_mass_to_charge\tsearch_engine_score[1]\tsearch_engine_score[2]\n"
            f"OSM\t{sequence}\t25.5\t-2\t500.2\t120\t{q_value}\n"
        )
    )
    if run_id:
        runs.finish(run_id, "finished")
    return results


def _results_page(app):
    app.switch_page("content/topp_workflow_results.py").run(timeout=30)
    assert not app.exception


def _downloads(app):
    return {item.proto.label for item in app.get("download_button")}


def test_configuration_survives_page_changes(app):
    app.run(timeout=30)
    next(
        button for button in app.button if button.label == "Load example dataset"
    ).click().run(timeout=30)
    # Follow the programmatic page switch explicitly in AppTest before interacting.
    app.switch_page("content/naseweis.py").run()
    app.number_input(key="naseweis-param-FDR_cutoff").set_value(0.01).run()
    # AppTest returns a fragment-only tree after a fragment widget changes.
    app.switch_page("content/topp_workflow_file_upload.py").run()
    app.switch_page("content/naseweis.py").run()
    app.checkbox(key="naseweis-param-add-decoys").uncheck().run()
    app.switch_page("content/topp_workflow_file_upload.py").run()
    app.switch_page("content/naseweis.py").run()
    assert not app.exception
    assert app.number_input(key="naseweis-param-FDR_cutoff").value == 0.01
    assert app.checkbox(key="naseweis-param-add-decoys").value is False
    saved = json.loads(
        (Path(app.session_state.workspace) / "naseweis/params.json").read_text()
    )
    assert saved["FDR_cutoff"] == 0.01 and saved["add-decoys"] is False


def test_saved_results_use_their_own_fdr_and_offer_all_exports(app):
    _finished_run(app)
    _results_page(app)
    table = app.dataframe[0].value
    assert table.loc[0, "q-value (%)"] == 1.0  # current config has FDR disabled
    assert table.loc[0, "Sequence"] == "ACGU"
    assert [str(metric.value) for metric in app.metric] == ["1", "1", "1"]
    assert _downloads(app) >= {
        "Download matches (CSV)",
        "Download mzTab",
        "Download idXML",
        "Download coverage (HTML)",
        "Download search parameters",
    }


def test_selecting_an_older_run_does_not_reuse_newer_results(app):
    earlier = _finished_run(app, "Earlier", sequence="ACGU", q_value="0.02").parent.name
    later = _finished_run(app, "Later", sequence="AAAA", q_value="-1").parent.name
    _results_page(app)
    picker = next(box for box in app.selectbox if box.label == "Run")
    for run_id, sequence, has_q in [
        (later, "AAAA", False),
        (earlier, "ACGU", True),
        (later, "AAAA", False),
    ]:
        # The widget values are run IDs; labels include a display name/date.
        picker.select(run_id).run()
        assert not app.exception
        table = app.dataframe[0].value
        assert table.loc[0, "Sequence"] == sequence
        assert ("q-value (%)" in table) is has_q
        picker = next(box for box in app.selectbox if box.label == "Run")


def test_legacy_results_are_still_readable_and_downloadable(app):
    _finished_run(app, legacy=True)
    _results_page(app)
    assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
    assert "Download legacy idXML" in _downloads(app)
    assert "Download mzTab" in _downloads(app)


def test_empty_results_are_an_explanation_not_an_exception(app):
    _finished_run(
        app,
        content="NUH\taccession\topt_sequence\nOLH\tsequence\taccession\tstart\tend\nOSH\tsequence\tcharge\n",
    )
    _results_page(app)
    assert all(str(metric.value) == "0" for metric in app.metric)
    assert any("No matches passed" in info.value for info in app.info)
    assert any("No RNA sequence coverage" in info.value for info in app.info)
    assert "Download mzTab" in _downloads(app)
    assert "Download matches (CSV)" not in _downloads(app)


def test_corrupt_results_show_an_actionable_error(app):
    _finished_run(app, content="OSH\tsequence\tcharge\nOSM\tACGU\n")
    _results_page(app)
    assert any(
        "Could not read this run's results" in error.value for error in app.error
    )
    assert not app.dataframe


def test_missing_result_file_explains_what_to_do_next(app):
    results = _finished_run(app)
    (results / "mztab_results/search.mzTab").unlink()
    _results_page(app)
    assert any("No search results yet" in info.value for info in app.info)


def test_missing_optional_export_does_not_hide_the_table(app):
    results = _finished_run(app)
    (results / "idxml_results/search.idXML").unlink()
    _results_page(app)
    assert len(app.dataframe[0].value) == 1
    assert "Download idXML" not in _downloads(app)
    assert "Download mzTab" in _downloads(app)


def test_invalid_coverage_does_not_prevent_match_export(app):
    results = _finished_run(app)
    file = results / "mztab_results/search.mzTab"
    file.write_text(file.read_text().replace("\trna\t1\t4", "\trna\t1\t99"))
    _results_page(app)
    assert any("Coverage is unavailable" in warning.value for warning in app.warning)
    assert "Download matches (CSV)" in _downloads(app)
    assert "Download coverage (HTML)" not in _downloads(app)


def test_first_visit_loads_settings_without_touching_personal_workspaces(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("LOCAL_WORKSPACES_DIR", str(tmp_path))
    fresh = AppTest.from_file("app.py").run(timeout=30)
    assert not fresh.exception
    assert fresh.title[0].value == "NASEWEIS"
    assert Path(fresh.session_state.workspace) == tmp_path / "default"


def test_upload_page_can_load_the_paired_example(app):
    app.run(timeout=30)
    app.switch_page("content/topp_workflow_file_upload.py").run()
    next(
        button for button in app.button if button.label == "Load example dataset"
    ).click().run()
    assert not app.exception
    saved = json.loads(
        (Path(app.session_state.workspace) / "naseweis/params.json").read_text()
    )
    assert Path(saved["mzML-files"]).name == "NME1.mzML"
    assert Path(saved["fasta-files"]).name == "NME1_Target.fasta"
    assert any("Continue to Configure" in message.value for message in app.success)


def test_missing_run_parameters_do_not_break_other_exports(app):
    results = _finished_run(app)
    (results.parent / "params.json").unlink()
    _results_page(app)
    assert len(app.dataframe[0].value) == 1
    assert "Download search parameters" not in _downloads(app)
    assert "Download mzTab" in _downloads(app)
