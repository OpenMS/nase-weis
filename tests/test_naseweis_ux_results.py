"""Out-of-order result/viewer journeys recover without losing saved evidence."""

# These fixtures are imported intentionally and injected by pytest.
# ruff: noqa: F811

import json
from pathlib import Path
from unittest.mock import Mock

import psutil
import pytest
from test_naseweis_experiment_integration_gui import experiment_results  # noqa: F401
from test_naseweis_viewers_gui import (
    acquisition,
    click,
    feature_table,
    peak_coordinates,
    plots,
    raw_page,
    saved_result,
)

from src.naseweis_experiment import sample_results
from src.workflow.RunManager import RunManager
from test_gui import _downloads, _results_page, app, topp_schemas  # noqa: F401


def test_unstarted_user_can_visit_results_raw_run_and_database_results_in_any_order(app):
    app.run(timeout=30)
    for page in (
        "content/topp_workflow_results.py",
        "content/naseweis_raw_viewer.py",
        "content/topp_workflow_execution.py",
        "content/topp_workflow_results.py",
    ):
        app.switch_page(page).run(timeout=30)
        assert not app.exception
        assert not app.get("plotly_chart")
        if page.endswith("topp_workflow_results.py"):
            assert any("No finished runs yet" in item.value for item in app.info)
        if page.endswith("naseweis_raw_viewer.py"):
            assert any("Upload an acquisition" in item.value for item in app.info)
    app.switch_page("content/decoy_database.py").run(timeout=30)
    app.radio(key="decoy-database-step").set_value("Results").run(timeout=30)
    assert not app.exception
    assert any("No finished runs yet" in item.value for item in app.info)
    workspace = Path(app.session_state.workspace)
    assert not list(workspace.glob("*/runs/*/run.json"))


def test_database_start_without_sequences_is_recoverable_and_creates_no_run(app, monkeypatch):
    worker = Mock()
    monkeypatch.setattr("src.workflow.WorkflowManager.WorkflowManager._start_workflow_local", worker)
    app.run(timeout=30)
    app.switch_page("content/decoy_database.py").run(timeout=30)
    app.radio(key="decoy-database-step").set_value("Run").run(timeout=30)
    next(button for button in app.button if button.label == "Start Workflow").click().run(timeout=30)
    assert not app.exception
    assert any("FASTA" in error.value and "Select" in error.value for error in app.error)
    assert RunManager(Path(app.session_state.workspace) / "decoy-database").list_runs() == []
    worker.assert_not_called()
    app.radio(key="decoy-database-step").set_value("Upload").run(timeout=30)
    assert not app.exception
    assert not app.error


def test_raw_viewer_recovers_from_bad_file_then_clears_deleted_evidence(app):
    app.run(timeout=30)
    source = Path(app.session_state.workspace) / "naseweis/input-files/mzML-files/recovered.mzML"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"not an mzML file")
    raw_page(app, source)
    click(app, "Load acquisition")
    assert any("Could not display this acquisition" in error.value for error in app.error)
    assert not plots(app)
    acquisition(source)
    app.run(timeout=30)
    assert not app.exception
    assert not app.error
    assert not plots(app), "a repaired file must be opened explicitly rather than reusing an old loaded state"
    click(app, "Load acquisition")
    assert len(plots(app)) == 3
    source.unlink()
    app.run(timeout=30)
    assert not app.exception
    assert not plots(app)
    assert any("Upload an acquisition" in item.value for item in app.info)


def test_raw_viewer_empty_acquisition_stays_an_explanation_after_page_changes(app):
    app.run(timeout=30)
    source = acquisition(Path(app.session_state.workspace) / "naseweis/input-files/mzML-files/empty.mzML", empty=True)
    raw_page(app, source)
    click(app, "Load acquisition")
    assert any("contains no spectra" in item.value for item in app.info)
    app.switch_page("content/topp_workflow_results.py").run(timeout=30)
    app.switch_page("content/naseweis_raw_viewer.py").run(timeout=30)
    assert not app.exception
    assert not app.get("plotly_chart")
    if any(button.label == "Load acquisition" for button in app.button):
        click(app, "Load acquisition")
    assert any("contains no spectra" in item.value for item in app.info)


@pytest.mark.parametrize("missing", ["source", "snapshot", "mzTab"])
def test_saved_result_recovers_after_required_file_is_restored(app, tmp_path, missing):
    result, original, snapshot = saved_result(app, tmp_path)
    path = {"source": original, "snapshot": snapshot, "mzTab": result / "mztab_results/search.mzTab"}[missing]
    content = path.read_bytes()
    path.unlink()
    _results_page(app)
    if missing == "mzTab":
        assert any("No search results yet" in item.value for item in app.info)
        assert "Download matches (CSV)" not in _downloads(app)
    else:
        assert any("Spectrum viewing is unavailable" in item.value for item in app.info)
        assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
        assert "Download coverage (SVG)" in _downloads(app)
        assert "Download mzTab" in _downloads(app)
    path.write_bytes(content)
    _results_page(app)
    assert not app.error
    assert not any("Spectrum viewing is unavailable" in item.value for item in app.info)
    assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
    click(app, "Open spectrum viewer")
    spectrum = next(figure for figure in plots(app) if figure["layout"].get("xaxis", {}).get("title", {}).get("text") == "m/z")
    assert peak_coordinates(spectrum) == [(100.0, 50.0), (200.0, 100.0)]


def test_deleted_acquisition_preserves_lfq_map_tables_and_coverage(app, tmp_path):
    result, original, _snapshot = saved_result(app, tmp_path)
    feature_table(result)
    original.unlink()
    _results_page(app)
    assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
    assert "Download coverage (SVG)" in _downloads(app)
    assert "Download LFQ feature table (TSV)" in _downloads(app)
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 1
    assert any("feature map is available" in item.value for item in app.info)
    assert not any(button.label == "Open spectrum viewer" for button in app.button)


@pytest.mark.parametrize("outcome", ["failed", "cancelled"])
def test_partial_experiment_keeps_sample_evidence_but_hides_complete_cohort_comparison(experiment_results, outcome):
    app, root, _inputs, _snapshots = experiment_results
    manifest_path = root / "experiment.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update(status="failed" if outcome == "failed" else "running", comparison_status="failed" if outcome == "failed" else "pending")
    manifest["samples"][1].update(status="failed" if outcome == "failed" else "running", error="MS1 extraction interrupted")
    manifest_path.write_text(json.dumps(manifest))
    runs = RunManager(root.parents[2])
    runs.finish(root.parent.name, outcome)
    _results_page(app)
    assert any("incomplete experiment" in item.value for item in app.warning)
    assert "Download abundance matrix (TSV)" not in _downloads(app)
    assert "Download matches (CSV)" in _downloads(app)
    assert "Download mzTab" in _downloads(app)
    assert "Download coverage (SVG)" in _downloads(app)
    picker = next(box for box in app.selectbox if box.label == "Acquisition details")
    picker.set_value("sample_002").run(timeout=30)
    assert not app.exception
    assert "Download matches (CSV)" in _downloads(app)
    assert sample_results(root, "sample_002").is_dir()


@pytest.mark.parametrize("outcome", ["failed", "cancelled", "running"])
def test_nonfinished_run_does_not_publish_a_stale_complete_comparison(experiment_results, outcome):
    app, root, _inputs, _snapshots = experiment_results
    RunManager(root.parents[2]).finish(root.parent.name, outcome)
    before = (root / "experiment.json").read_bytes()
    _results_page(app)
    assert "Download matches (CSV)" in _downloads(app)
    assert "Download abundance matrix (TSV)" not in _downloads(app)
    assert any("incomplete experiment" in item.value for item in app.warning)
    if outcome == "running":
        assert any("was interrupted" in item.value for item in app.warning)
    assert (root / "experiment.json").read_bytes() == before


def test_active_experiment_is_hidden_until_its_worker_has_stopped(experiment_results, monkeypatch):
    app, root, _inputs, _snapshots = experiment_results
    RunManager(root.parents[2]).finish(root.parent.name, "running")
    monkeypatch.setattr("src.Workflow.Workflow.get_workflow_status", lambda _self: {"running": True})
    _results_page(app)
    assert any("No finished runs yet" in item.value for item in app.info)
    assert "Download matches (CSV)" not in _downloads(app)
    monkeypatch.setattr("src.Workflow.Workflow.get_workflow_status", lambda _self: {"running": False})
    app.run(timeout=30)
    assert not app.exception
    assert "Download matches (CSV)" in _downloads(app)
    assert "Download abundance matrix (TSV)" not in _downloads(app)
    assert any("was interrupted" in item.value for item in app.warning)


def test_orphaned_pid_after_restart_does_not_hide_partial_results(experiment_results, monkeypatch):
    app, root, _inputs, _snapshots = experiment_results
    runs = RunManager(root.parents[2])
    runs.finish(root.parent.name, "running")
    marker = runs.workflow_dir / "pids/54321"
    marker.parent.mkdir()
    marker.touch()
    before = (root / "experiment.json").read_bytes()
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", Mock(side_effect=psutil.NoSuchProcess(54321)))
    _results_page(app)
    assert "Download matches (CSV)" in _downloads(app)
    assert "Download abundance matrix (TSV)" not in _downloads(app)
    assert any("was interrupted" in item.value for item in app.warning)
    assert not marker.exists()
    assert (root / "experiment.json").read_bytes() == before


def test_results_workspace_switch_cannot_display_another_workspaces_saved_run(app, tmp_path):
    saved_result(app, tmp_path)
    _results_page(app)
    assert "Download matches (CSV)" in _downloads(app)
    original = Path(app.session_state.workspace)
    other = original.parent / "empty-results-workspace"
    other.mkdir()
    app.run(timeout=30)
    next(box for box in app.sidebar.selectbox if box.label == "choose existing workspace").set_value(other.name).run(timeout=30)
    app.run(timeout=30)
    assert not app.exception
    assert Path(app.session_state.workspace) == other
    assert not app.get("plotly_chart")
    assert not app.dataframe
    assert "Download matches (CSV)" not in _downloads(app)
    assert any("No finished runs yet" in item.value for item in app.info)
    next(box for box in app.sidebar.selectbox if box.label == "choose existing workspace").set_value(original.name).run(timeout=30)
    app.run(timeout=30)
    assert not app.exception
    assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
    assert "Download matches (CSV)" in _downloads(app)
