"""Out-of-order starts explain missing prerequisites without creating jobs."""

# Shared application fixtures intentionally use their exported names.
# ruff: noqa: F811

import json
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from test_naseweis_viewers_gui import acquisition

from src.workflow.ExperimentalDesign import save_design
from src.workflow.RunManager import FINISHED, RunManager
from test_gui import app, topp_schemas  # noqa: F401


def _button(app, label):
    return next(button for button in app.button if button.label == label)


def _run_page(app):
    app.switch_page("content/topp_workflow_execution.py").run(timeout=30)
    assert not app.exception


@pytest.fixture
def prerequisites(app, monkeypatch, topp_schemas):
    """Inputs exist without visiting Upload or Configure; a worker never runs."""
    app.run(timeout=30)
    directory = Path(app.session_state.workspace) / "naseweis"
    directory.mkdir(parents=True, exist_ok=True)
    input_directory = directory / "input-files/mzML-files"
    input_directory.mkdir(parents=True, exist_ok=True)
    first = acquisition(input_directory / "control.mzML")
    second = acquisition(input_directory / "treated.mzML", second=True)
    database = directory / "input-files/fasta-files/rna.fasta"
    database.parent.mkdir(parents=True, exist_ok=True)
    database.write_text(">rna\nACGUACGU\n")
    params = {
        "mzML-files": [str(first)],
        "fasta-files": str(database),
        "enable-lfq": False,
        "add-decoys": True,
        "FDR_cutoff": 0.05,
    }
    # No process may escape a negative UX test, even if validation regresses.
    worker = Mock()
    # Schema reads are real; process dispatch is stubbed, so local tool PATH
    # availability must not decide whether a valid UI start is accepted.
    monkeypatch.setattr("src.Workflow.shutil.which", lambda tool: f"/test-tools/{tool}")
    monkeypatch.setattr(
        "src.workflow.WorkflowManager.WorkflowManager._start_workflow_local", worker
    )
    return directory, [first, second], database, params, worker


def _save(directory, params):
    path = directory / "params.json"
    path.write_text(json.dumps(params))
    return path


def _lfq_design(directory, paths, *, blank_condition=False, fraction=1):
    design = pd.DataFrame({
        "file": [path.name for path in paths],
        "condition": ["control", "" if blank_condition else "treated"],
        "replicate": [1, 1],
        "fraction": [1, fraction],
    })
    save_design(directory, design, [str(path) for path in paths])


@pytest.mark.parametrize("case", [
    "nothing_selected", "database_only", "spectra_only", "deleted_acquisition",
    "deleted_database", "deleted_digest", "acquisition_is_directory",
    "unsupported_extension", "malformed_selection", "duplicate_path",
    "duplicate_basename", "lfq_one_acquisition", "lfq_blank_condition",
    "lfq_fractionated", "lfq_ms2_only", "lfq_missing_fit", "lfq_bad_normalization",
    "lfq_malformed_fit_list", "lfq_malformed_fit_dict",
    "invalid_design_number",
])
def test_out_of_order_or_incomplete_start_has_no_phantom_run(
    app, prerequisites, case
):
    directory, paths, database, params, worker = prerequisites
    if case == "nothing_selected":
        params = {"mzML-files": [], "enable-lfq": False}
    elif case == "database_only":
        params["mzML-files"] = []
    elif case == "spectra_only":
        params.pop("fasta-files")
    elif case == "deleted_acquisition":
        paths[0].unlink()
    elif case == "deleted_database":
        database.unlink()
    elif case == "deleted_digest":
        params.update({"database-source": "Cached digest", "digest-files": str(directory / "missing.oms")})
    elif case == "acquisition_is_directory":
        paths[0].unlink()
        paths[0].mkdir()
    elif case == "unsupported_extension":
        other = paths[0].with_suffix(".txt")
        other.write_bytes(paths[0].read_bytes())
        params["mzML-files"] = [str(other)]
    elif case == "malformed_selection":
        params["mzML-files"] = [str(paths[0]), None]
    elif case == "duplicate_path":
        params["mzML-files"] = [str(paths[0]), str(paths[0])]
    elif case == "duplicate_basename":
        other = directory / "other" / paths[0].name
        other.parent.mkdir()
        other.write_bytes(paths[0].read_bytes())
        params["mzML-files"] = [str(paths[0]), str(other)]
    else:
        params.update({
            "enable-lfq": True,
            "mzML-files": [str(path) for path in paths],
            "lfq-measure": "intensity",
            "lfq-normalization": "total",
        })
        if case == "lfq_one_acquisition":
            params["mzML-files"] = [str(paths[0])]
        _lfq_design(
            directory, paths, blank_condition=case == "lfq_blank_condition",
            fraction=2 if case == "lfq_fractionated" else 1,
        )
        if case == "lfq_ms2_only":
            acquisition(paths[1], ms2_only=True)
        elif case == "lfq_missing_fit":
            params.update({"lfq-measure": "model_area", "FeatureFinderMetaboIdent": {"model:type": "none"}})
        elif case in ("lfq_malformed_fit_list", "lfq_malformed_fit_dict"):
            malformed_model = ["symmetric"] if case.endswith("list") else {"type": "symmetric"}
            params.update({"lfq-measure": "model_area", "FeatureFinderMetaboIdent": {"model:type": malformed_model}})
        elif case == "lfq_bad_normalization":
            params["lfq-normalization"] = "unsupported"
        elif case == "invalid_design_number":
            design_file = directory / "experimental-design.tsv"
            lines = design_file.read_text().splitlines()
            row = lines[1].split("\t")
            row[1] = "1.5"
            lines[1] = "\t".join(row)
            design_file.write_text("\n".join(lines) + "\n")
    _save(directory, params)
    original_design = (
        (directory / "experimental-design.tsv").read_bytes()
        if (directory / "experimental-design.tsv").exists() else None
    )
    _run_page(app)
    app.text_input(key="runs-naseweis-next-name").set_value("My RNA experiment").run(timeout=30)
    _button(app, "Start Workflow").click().run(timeout=30)
    assert not app.exception
    assert app.error, f"{case} did not show actionable validation"
    assert not any("Traceback" in error.value for error in app.error)
    worker.assert_not_called()
    assert RunManager(directory).list_runs() == []
    assert app.text_input(key="runs-naseweis-next-name").value == "My RNA experiment"
    assert "_sample_design" not in json.loads((directory / "params.json").read_text())
    if original_design is not None:
        assert (directory / "experimental-design.tsv").read_bytes() == original_design


def test_failed_start_can_be_corrected_without_losing_run_name_or_previous_results(
    app, prerequisites
):
    directory, paths, _database, params, worker = prerequisites
    original_database = params.pop("fasta-files")
    _save(directory, params)
    runs = RunManager(directory)
    prior = runs.new_run("Previous successful analysis", directory / "params.json")
    marker = runs.results_dir(prior) / "keep.txt"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text("previous result evidence")
    runs.finish(prior, FINISHED)
    _run_page(app)
    app.text_input(key="runs-naseweis-next-name").set_value("Recovered experiment").run(timeout=30)
    _button(app, "Start Workflow").click().run(timeout=30)
    assert not app.exception
    assert app.error
    worker.assert_not_called()
    assert len(runs.list_runs()) == 1
    assert marker.read_text() == "previous result evidence"
    assert app.text_input(key="runs-naseweis-next-name").value == "Recovered experiment"

    # Correct the prerequisite through Configure, just as an end user would.
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert not app.exception
    app.selectbox(key="naseweis-param-fasta-files").set_value(original_database).run(timeout=30)
    _run_page(app)
    assert app.text_input(key="runs-naseweis-next-name").value == "Recovered experiment"

    def complete_without_tools():
        accepted = next(run for run in runs.list_runs() if run["name"] == "Recovered experiment")
        runs.finish(accepted["id"], FINISHED)

    worker.side_effect = complete_without_tools
    _button(app, "Start Workflow").click().run(timeout=30)
    assert not app.exception
    assert not app.error
    assert worker.call_count == 1
    assert len(runs.list_runs()) == 2
    names = [run["name"] for run in runs.list_runs()]
    assert "Recovered experiment" in names
    assert "Previous successful analysis" in names
    latest = next(run for run in runs.list_runs() if run["name"] == "Recovered experiment")
    snapshot = json.loads(runs.params_file(latest["id"]).read_text())
    assert snapshot["mzML-files"] == [str(paths[0])]
    assert snapshot["fasta-files"] == original_database
    assert marker.read_text() == "previous result evidence"


def test_repeated_invalid_start_does_not_consume_the_requested_name(app, prerequisites):
    directory, _paths, _database, params, worker = prerequisites
    params["mzML-files"] = []
    _save(directory, params)
    _run_page(app)
    app.text_input(key="runs-naseweis-next-name").set_value("Keep this name").run(timeout=30)
    for _ in range(2):
        _button(app, "Start Workflow").click().run(timeout=30)
        assert not app.exception
        assert app.error
        assert app.text_input(key="runs-naseweis-next-name").value == "Keep this name"
    assert RunManager(directory).list_runs() == []
    worker.assert_not_called()


def test_failed_start_keeps_run_name_when_visiting_configure_to_fix_inputs(app, prerequisites):
    directory, _paths, _database, params, worker = prerequisites
    params["mzML-files"] = []
    _save(directory, params)
    _run_page(app)
    app.text_input(key="runs-naseweis-next-name").set_value("Keep through correction").run(timeout=30)
    _button(app, "Start Workflow").click().run(timeout=30)
    assert app.error
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert not app.exception
    _run_page(app)
    assert app.text_input(key="runs-naseweis-next-name").value == "Keep through correction"
    assert RunManager(directory).list_runs() == []
    worker.assert_not_called()


def test_run_name_drafts_are_separate_for_each_workspace(tmp_path):
    harness = """
from pathlib import Path
import streamlit as st
from src.workflow.RunManager import RunManager
from src.workflow.run_widgets import run_name_input
space = st.selectbox('Workspace', ['first', 'second'], key='workspace_choice')
run_name_input(RunManager(Path(st.session_state.root) / space / 'naseweis'))
"""
    ui = AppTest.from_string(harness)
    ui.session_state.root = str(tmp_path)
    ui.run(timeout=30)
    assert not ui.exception
    ui.text_input(key="runs-naseweis-next-name").set_value("First experiment").run(timeout=30)
    ui.selectbox(key="workspace_choice").set_value("second").run(timeout=30)
    assert ui.text_input(key="runs-naseweis-next-name").value == ""
    ui.text_input(key="runs-naseweis-next-name").set_value("Second experiment").run(timeout=30)
    ui.selectbox(key="workspace_choice").set_value("first").run(timeout=30)
    assert ui.text_input(key="runs-naseweis-next-name").value == "First experiment"
    ui.selectbox(key="workspace_choice").set_value("second").run(timeout=30)
    assert ui.text_input(key="runs-naseweis-next-name").value == "Second experiment"
