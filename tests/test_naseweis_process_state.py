"""Restart-orphaned PID files must not block valid UI recovery."""

import importlib
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import psutil
import pytest
from streamlit.testing.v1 import AppTest

from src.workflow.process_state import local_workflow_running


@pytest.fixture
def markers(tmp_path):
    folder = tmp_path / "pids"
    folder.mkdir()
    marker = folder / "54321"
    marker.touch()
    return folder, marker


@pytest.mark.parametrize("state,running", [
    (psutil.STATUS_SLEEPING, True),
    (psutil.STATUS_STOPPED, True),
    (psutil.STATUS_ZOMBIE, False),
    (psutil.STATUS_DEAD, False),
    ("gone", False),
    ("not-running", False),
    ("denied", True),
    ("invalid", False),
    ("overflow", False),
])
def test_pid_markers_follow_process_liveness_not_file_existence(markers, monkeypatch, state, running):
    folder, marker = markers
    process = Mock()
    process.is_running.return_value = state != "not-running"
    process.status.return_value = state
    lookup = Mock(return_value=process)
    if state == "gone":
        lookup.side_effect = psutil.NoSuchProcess(54321)
    elif state == "denied":
        process.status.side_effect = psutil.AccessDenied(54321)
    elif state == "invalid":
        lookup.side_effect = ValueError("invalid PID")
    elif state == "overflow":
        lookup.side_effect = OverflowError("PID outside platform range")
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", lookup)
    assert local_workflow_running(folder) is running
    assert marker.exists() is running
    lookup.assert_called_once_with(54321)


def test_live_worker_does_not_prevent_pruning_other_dead_markers(markers, monkeypatch):
    folder, marker = markers
    dead = folder / "54322"
    dead.touch()
    process = Mock()
    process.is_running.return_value = True
    process.status.return_value = psutil.STATUS_SLEEPING

    def lookup(pid):
        if pid == 54321:
            return process
        raise psutil.NoSuchProcess(pid)

    monkeypatch.setattr("src.workflow.process_state.psutil.Process", lookup)
    assert local_workflow_running(folder)
    assert marker.exists()
    assert not dead.exists()


def test_unrelated_files_numeric_directories_and_symlinks_are_not_pid_markers(tmp_path, monkeypatch):
    folder = tmp_path / "pids"
    folder.mkdir()
    note = folder / "conversion.json"
    note.write_text('{"state":"prepared"}')
    numeric_directory = folder / "54321"
    numeric_directory.mkdir()
    symbolic = folder / "54322"
    symbolic.touch()
    original = Path.is_symlink
    monkeypatch.setattr(Path, "is_symlink", lambda path: path == symbolic or original(path))
    lookup = Mock()
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", lookup)
    assert local_workflow_running(folder) is False
    assert note.read_text() == '{"state":"prepared"}'
    assert numeric_directory.is_dir()
    assert symbolic.exists()
    lookup.assert_not_called()


def test_missing_pid_directory_is_idle(tmp_path):
    assert local_workflow_running(tmp_path / "missing") is False


def test_unreadable_pid_directory_is_conservatively_busy(tmp_path, monkeypatch):
    folder = tmp_path / "pids"
    native = Path.iterdir

    def read(path):
        if path == folder:
            raise PermissionError("temporarily unavailable")
        return native(path)

    monkeypatch.setattr(Path, "iterdir", read)
    assert local_workflow_running(folder) is True


def test_dead_marker_cleanup_failure_does_not_invent_a_live_worker(markers, monkeypatch):
    folder, marker = markers
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", Mock(side_effect=psutil.NoSuchProcess(54321)))
    native = Path.unlink

    def unlink(path, *args, **kwargs):
        if path == marker:
            raise PermissionError("read-only marker")
        return native(path, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", unlink)
    assert local_workflow_running(folder) is False
    assert marker.exists()


@pytest.fixture
def manager(markers):
    cls = importlib.import_module("src.workflow.WorkflowManager").WorkflowManager
    instance = object.__new__(cls)
    instance.workflow_dir = markers[0].parent
    instance.executor = SimpleNamespace(pid_dir=markers[0])
    instance._queue_manager = None
    return instance


def test_manager_status_cleans_only_dead_pids_preserving_run_metadata_and_logs(manager, markers, monkeypatch):
    folder, marker = markers
    log = folder.parent / "run.log"
    record = folder.parent / "run.json"
    log.write_text("native partial output\n")
    record.write_text('{"status":"running"}')
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", Mock(side_effect=psutil.NoSuchProcess(54321)))
    assert manager.get_workflow_status()["running"] is False
    assert not marker.exists()
    assert log.read_text() == "native partial output\n"
    assert record.read_text() == '{"status":"running"}'


def test_manager_status_keeps_live_workers_running(manager, monkeypatch):
    process = Mock()
    process.is_running.return_value = True
    process.status.return_value = psutil.STATUS_SLEEPING
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", Mock(return_value=process))
    assert manager.get_workflow_status()["running"] is True


@pytest.mark.parametrize("state,expected", [("started", True), ("queued", True), ("finished", False), ("canceled", False)])
def test_queue_status_remains_authoritative(manager, markers, monkeypatch, state, expected):
    job = SimpleNamespace(
        status=SimpleNamespace(value=state), progress=0, current_step=None,
        queue_position=None, queue_length=0, enqueued_at=None, started_at=None, result=None, error=None,
    )
    queue = Mock()
    queue.is_available = True
    queue.load_job_id.return_value = "test-job"
    queue.get_job_info.return_value = job
    manager._queue_manager = queue
    lookup = Mock(side_effect=psutil.NoSuchProcess(54321))
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", lookup)
    assert manager.get_workflow_status()["running"] is expected
    assert markers[1].exists()
    lookup.assert_not_called()


HARNESS = """
from pathlib import Path
from types import SimpleNamespace
import streamlit as st
from src.workflow.StreamlitUI import StreamlitUI

ui = object.__new__(StreamlitUI)
ui.workflow_dir = Path(st.session_state.folder)
ui.executor = SimpleNamespace(pid_dir=ui.workflow_dir / 'pids')
ui.export_parameters_markdown = lambda: 'Saved test settings'
def start():
    st.session_state.start_count = st.session_state.get('start_count', 0) + 1
provider = None if st.session_state.status is None else lambda: st.session_state.status
ui.execution_section(start, get_status_function=provider)
"""


@pytest.mark.parametrize("status", [None, {}, {"running": False, "status": "idle"}])
def test_dead_pid_files_do_not_block_start_button(markers, monkeypatch, status):
    folder, marker = markers
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", Mock(side_effect=psutil.NoSuchProcess(54321)))
    app = AppTest.from_string(HARNESS)
    app.session_state.folder = str(folder.parent)
    app.session_state.status = status
    app.run(timeout=30)
    assert not app.exception
    assert any(button.label == "Start Workflow" for button in app.button)
    assert not any(button.label == "Stop Workflow" for button in app.button)
    next(button for button in app.button if button.label == "Start Workflow").click().run(timeout=30)
    assert not app.exception
    assert app.session_state.start_count == 1
    if status is None or "running" not in status:
        assert not marker.exists()


def test_ui_does_not_override_explicit_idle_status_with_a_pid_marker(markers, monkeypatch):
    folder, marker = markers
    lookup = Mock()
    monkeypatch.setattr("src.workflow.process_state.psutil.Process", lookup)
    app = AppTest.from_string(HARNESS)
    app.session_state.folder = str(folder.parent)
    app.session_state.status = {"running": False, "status": "finished"}
    app.run(timeout=30)
    assert not app.exception
    assert any(button.label == "Start Workflow" for button in app.button)
    assert not any(button.label == "Stop Workflow" for button in app.button)
    lookup.assert_not_called()
    assert marker.exists()
