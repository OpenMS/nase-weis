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
        "content/decoy_database.py",
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


def test_nightly_configuration_preserves_exports_and_cached_digest_mode(app):
    app.run(timeout=30)
    app.switch_page("content/naseweis.py").run(timeout=30)
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    exports = next(
        widget for widget in app.multiselect if widget.label == "Result files"
    )
    exports.set_value(["id_out", "db_out", "lfq_out"]).run()
    app.switch_page("content/topp_workflow_file_upload.py").run()
    app.switch_page("content/naseweis.py").run()
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    app.selectbox(key="naseweis-param-database-source").select("Cached digest").run()
    assert not app.exception
    assert any("cached database controls digestion" in info.value for info in app.info)
    saved = json.loads(
        (Path(app.session_state.workspace) / "naseweis/params.json").read_text()
    )
    assert saved["search-exports"] == ["id_out", "db_out", "lfq_out"]
    assert saved["database-source"] == "Cached digest"


def test_all_available_search_artifacts_can_be_downloaded(app):
    from src.naseweis_tools import SEARCH_EXPORTS

    results = _finished_run(app)
    for label, relative, mime in SEARCH_EXPORTS.values():
        path = results / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("test artifact")
    (results / "database").mkdir()
    (results / "database/decoys.fasta").write_text(">DECOY_rna\nUGCA\n")
    _results_page(app)
    downloads = _downloads(app)
    for option, (label, relative, mime) in SEARCH_EXPORTS.items():
        assert (
            "Download idXML" if option == "id_out" else f"Download {label}"
        ) in downloads
    assert "Download decoys.fasta" in downloads


def test_database_page_shows_saved_counts_exports_and_parameter_snapshot(app):
    from src.workflow.RunManager import RunManager

    app.run(timeout=30)
    app.switch_page("content/decoy_database.py").run(timeout=30)
    assert not app.exception
    app.radio(key="decoy-database-step").set_value("Configure").run(timeout=30)
    assert not app.toggle(key="advanced").value
    assert not any(
        widget.label == "Relevant proteins (optional)" for widget in app.selectbox
    )
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert not app.exception
    assert any(
        widget.key == "decoy-database-TOPP-DecoyDatabase:1:decoy_string"
        for widget in app.text_input
    )
    directory = Path(app.session_state.workspace) / "decoy-database"
    runs = RunManager(directory)
    run_id = runs.new_run("Saved database", directory / "params.json")
    database = runs.results_dir(run_id) / "database"
    database.mkdir(parents=True)
    (database / "database.fasta").write_text(">rna\nACGU\n>DECOY_rna\nUGCA\n")
    runs.finish(run_id, "finished")
    # Refresh the full page so that the run picker reads its new persisted run.
    app.switch_page("content/quickstart.py").run()
    app.switch_page("content/decoy_database.py").run(timeout=30)
    assert not app.exception
    app.radio(key="decoy-database-step").set_value("Results").run(timeout=30)
    assert any(str(metric.value) == "2 sequences" for metric in app.metric)
    assert "Download database.fasta" in _downloads(app)
    assert "Download database parameters" in _downloads(app)
    runs.params_file(run_id).unlink()
    app.switch_page("content/quickstart.py").run()
    app.switch_page("content/decoy_database.py").run(timeout=30)
    app.radio(key="decoy-database-step").set_value("Results").run(timeout=30)
    assert "Download database parameters" not in _downloads(app)
    (database / "database.fasta").unlink()
    app.run(timeout=30)
    assert any("No generated databases" in info.value for info in app.info)


def test_database_step_change_refreshes_files_added_on_upload(app):
    app.run(timeout=30)
    app.switch_page("content/decoy_database.py").run(timeout=30)
    directory = Path(app.session_state.workspace) / "decoy-database"
    uploaded = directory / "input-files/database-files/new-rna.fasta"
    uploaded.parent.mkdir(parents=True, exist_ok=True)
    uploaded.write_text(">rna\nACGU\n")
    app.radio(key="decoy-database-step").set_value("Configure").run(timeout=30)
    assert not app.exception
    selector = app.multiselect(key="decoy-database-param-database-files")
    assert "new-rna.fasta" in selector.options
    selector.set_value([str(uploaded)]).run()
    app.radio(key="decoy-database-step").set_value("Run").run(timeout=30)
    assert any(button.label == "Start Workflow" for button in app.button)
    saved = json.loads((directory / "params.json").read_text())
    assert saved["database-files"] == [str(uploaded)]
    app.radio(key="decoy-database-step").set_value("Results").run(timeout=30)
    assert any("No finished runs" in info.value for info in app.info)


def test_rnome_controls_require_advanced_and_keep_saved_values(app):
    app.run(timeout=30)
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert not app.toggle(key="advanced").value
    prefix = "naseweis-TOPP-NucleicAcidSearchEngine:1:"

    def keys():
        return {
            widget.key
            for kind in [
                "selectbox",
                "checkbox",
                "number_input",
                "text_input",
                "text_area",
                "multiselect",
            ]
            for widget in app.get(kind)
        }

    advanced_keys = {
        prefix + name
        for name in [
            "scoring:method",
            "scoring:single_peak_credit",
            "precursor:correct_with_flash",
            "precursor:flash:min_cos",
            "preprocessing:deisotoping",
            "preprocessing:flash:min_cos",
            "preprocessing:filter_window_mower",
            "report:precursor_percentage",
            "modifications:resolve_ambiguities",
            "scoring:blacklist_mz",
        ]
    }
    assert not keys() & advanced_keys
    assert not any(widget.label == "Result files" for widget in app.multiselect)
    # The searchable dialog must obey the toggle too.
    app.button(key="all-params-NucleicAcidSearchEngine").click().run(timeout=30)
    assert not any("scoring:" in str(item.value) for item in app.markdown)
    app.switch_page("content/naseweis.py").run(timeout=30)
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert not app.exception
    assert advanced_keys <= keys()
    app.selectbox(key=prefix + "scoring:method").select("mvh").run(timeout=30)
    app.switch_page("content/naseweis.py").run(timeout=30)
    app.toggle(key="advanced").set_value(False).run(timeout=30)
    assert not keys() & advanced_keys
    saved = json.loads(
        (Path(app.session_state.workspace) / "naseweis/params.json").read_text()
    )
    assert saved["NucleicAcidSearchEngine"]["scoring:method"] == "mvh"
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert app.selectbox(key=prefix + "scoring:method").value == "mvh"
    assert not app.exception
