"""Container uploads remain usable when no desktop file dialog is available."""

import importlib
from pathlib import Path
from unittest.mock import Mock

import pytest
from streamlit.testing.v1 import AppTest

from test_gui import app as naseweis_app  # noqa: F401

HARNESS = """
from pathlib import Path
import streamlit as st
from src.workflow.ParameterManager import ParameterManager
from src.workflow.StreamlitUI import StreamlitUI

workspace = Path(st.session_state.test_workspace)
workspace.mkdir(parents=True, exist_ok=True)
ui = StreamlitUI(workspace, None, None, ParameterManager(workspace))
ui.upload_widget("mzML-files", ["mzML", "raw", "d"], name="MS data", list_files=False)
"""


@pytest.fixture
def picker_modules():
    # Several inherited tests temporarily replace these modules at collection.
    # Import at execution so this harness always uses the real Streamlit UI.
    return (
        importlib.import_module("src.common.common"),
        importlib.import_module("src.workflow.StreamlitUI"),
    )


@pytest.mark.parametrize(
    "platform,tk,display,expected",
    [
        ("linux", False, None, False),
        ("linux", False, ":1", False),
        ("linux", True, None, False),
        ("linux", True, "", False),
        ("linux", True, ":1", True),
        ("darwin", True, None, True),
        ("win32", True, None, True),
        ("darwin", False, None, False),
        ("win32", False, None, False),
    ],
)
def test_native_dialog_capability(
    picker_modules, monkeypatch, platform, tk, display, expected
):
    common, _ui = picker_modules
    monkeypatch.setattr(common, "OS_PLATFORM", platform)
    monkeypatch.setattr(common, "TK_AVAILABLE", tk)
    if display is None:
        monkeypatch.delenv("DISPLAY", raising=False)
    else:
        monkeypatch.setenv("DISPLAY", display)
    assert common.native_dialog_available() is expected


@pytest.fixture
def picker_app(tmp_path, monkeypatch, picker_modules):
    common, _ui = picker_modules
    monkeypatch.setattr(common, "OS_PLATFORM", "linux")
    monkeypatch.setattr(common, "TK_AVAILABLE", False)
    monkeypatch.delenv("DISPLAY", raising=False)

    def create(*, local=True, mount=None):
        app = AppTest.from_string(HARNESS)
        app.session_state.test_workspace = str(tmp_path / "workspace")
        app.session_state.location = "local" if local else "online"
        app.session_state.settings = {"local_data_dir": str(mount) if mount else ""}
        app.session_state.local_dir = ""
        app.session_state.previous_dir = str(tmp_path)
        app.run(timeout=30)
        assert not app.exception
        return app

    return create


def test_headless_local_offers_browser_upload_and_manual_path_without_cloud_warning(
    picker_app,
):
    app = picker_app()
    assert len(app.get("file_uploader")) == 1
    assert app.text_input(key="path_to_folder_mzML-files").value == ""
    assert app.button(key="add_files_from_local_mzML-files")
    assert not any(
        "streamlit cloud" in warning.value.lower() or "tkinter" in warning.value.lower()
        for warning in app.warning
    )
    native = [
        button
        for button in app.button
        if button.key == "local_browse_single_mzML-files"
    ]
    assert not native or all(button.disabled for button in native)


def test_importable_tk_without_linux_display_still_offers_upload(
    picker_app, picker_modules, monkeypatch
):
    common, _ui = picker_modules
    monkeypatch.setattr(common, "TK_AVAILABLE", True)
    app = picker_app()
    assert len(app.get("file_uploader")) == 1
    assert app.text_input(key="path_to_folder_mzML-files")


def test_desktop_native_picker_adds_selected_file(
    picker_app, picker_modules, monkeypatch, tmp_path
):
    common, ui = picker_modules
    monkeypatch.setattr(common, "TK_AVAILABLE", True)
    monkeypatch.setattr(common, "OS_PLATFORM", "darwin")
    source = tmp_path / "picked.mzML"
    source.write_bytes(b"native picker acquisition")
    picker = Mock(return_value=[str(source)])
    monkeypatch.setattr(ui, "tk_file_dialog", picker)
    app = picker_app()
    assert not app.get("file_uploader")
    assert app.text_input(key="path_to_folder_mzML-files")
    app.button(key="local_browse_single_mzML-files").click().run(timeout=30)
    assert not app.exception
    picker.assert_called_once()
    copied = tmp_path / "workspace/input-files/mzML-files/picked.mzML"
    assert copied.read_bytes() == source.read_bytes()


def test_online_upload_does_not_offer_server_path_or_native_picker(picker_app):
    app = picker_app(local=False)
    assert len(app.get("file_uploader")) == 1
    assert not app.text_input
    assert not any(str(button.key).startswith("local_browse") for button in app.button)


@pytest.mark.parametrize(
    "copy", [True, False], ids=["default-copy", "external-reference"]
)
def test_headless_manual_folder_addition_supports_copy_and_reference(
    picker_app, tmp_path, copy
):
    source = tmp_path / "source"
    source.mkdir()
    (source / "one.mzML").write_bytes(b"mzML acquisition")
    (source / "two.raw").write_bytes(b"vendor acquisition")
    (source / "bundle.d").mkdir()
    (source / "bundle.d/analysis.bin").write_bytes(b"Bruker data")
    (source / "ignore.txt").write_text("unrelated")
    app = picker_app()
    assert app.checkbox(key="mzML-files-copy_files").value is True
    if not copy:
        app.checkbox(key="mzML-files-copy_files").uncheck().run(timeout=30)
    app.text_input(key="path_to_folder_mzML-files").set_value(str(source)).run(
        timeout=30
    )
    app.button(key="add_files_from_local_mzML-files").click().run(timeout=30)
    assert not app.exception
    destination = tmp_path / "workspace/input-files/mzML-files"
    if copy:
        assert (destination / "one.mzML").read_bytes() == b"mzML acquisition"
        assert (destination / "two.raw").read_bytes() == b"vendor acquisition"
        assert (destination / "bundle.d/analysis.bin").read_bytes() == b"Bruker data"
        assert not (destination / "ignore.txt").exists()
    else:
        referenced = (destination / "external_files.txt").read_text().splitlines()
        assert set(referenced) == {
            str(source / name) for name in ["one.mzML", "two.raw", "bundle.d"]
        }
        assert not (destination / "one.mzML").exists()


@pytest.mark.parametrize("invalid", ["blank", "missing", "file", "unreadable"])
def test_invalid_manual_folder_is_a_user_message(
    picker_app, tmp_path, monkeypatch, invalid
):
    directory = tmp_path / "unreadable"
    directory.mkdir()
    paths = {
        "blank": "",
        "missing": str(tmp_path / "missing"),
        "file": str(tmp_path / "not-directory.mzML"),
        "unreadable": str(directory),
    }
    (tmp_path / "not-directory.mzML").touch()
    if invalid == "unreadable":
        original = Path.iterdir

        def iterdir(path):
            if path == directory:
                raise PermissionError("Permission denied for test folder")
            return original(path)

        monkeypatch.setattr(Path, "iterdir", iterdir)
    app = picker_app()
    app.text_input(key="path_to_folder_mzML-files").set_value(paths[invalid]).run(
        timeout=30
    )
    add = app.button(key="add_files_from_local_mzML-files")
    if invalid == "blank" and add.disabled:
        assert not app.exception
        return
    add.click().run(timeout=30)
    assert not app.exception
    assert app.error or app.warning
    assert not any("Successfully" in success.value for success in app.success)
    destination = tmp_path / "workspace/input-files/mzML-files"
    assert not destination.exists() or not any(destination.iterdir())


@pytest.mark.parametrize("local", [True, False], ids=["local", "online"])
def test_mounted_root_is_browsable_refreshes_and_replaces_manual_paths(
    picker_app, monkeypatch, tmp_path, local
):
    mount = tmp_path / "mounted-data"
    mount.mkdir()
    (mount / "mounted.mzML").write_bytes(b"mounted acquisition")
    (mount / "batch").mkdir()
    (mount / "batch/nested.mzML").write_bytes(b"nested mounted acquisition")
    (mount / "bundle.d").mkdir()
    (mount / "bundle.d/analysis.bin").write_bytes(b"Bruker data")
    monkeypatch.setattr(
        "src.workflow.StreamlitUI.os.path.ismount", lambda path: Path(path) == mount
    )
    app = picker_app(local=local, mount=mount)
    assert len(app.get("file_uploader")) == 1
    assert not app.text_input
    assert not any(item.key == "mzML-files-copy_files" for item in app.checkbox)
    assert not any(button.key == "local_browse_mzML-files" for button in app.button)
    assert app.button(key="mounted_refresh_mzML-files")
    assert any("mounted.mzML" in item.label for item in app.checkbox)
    app.button(key="mounted_dir_mzML-files_batch").click().run(timeout=30)
    assert not app.exception
    assert any("nested.mzML" in item.label for item in app.checkbox)
    app.button(key="mounted_parent_mzML-files").click().run(timeout=30)
    assert not app.exception
    assert any("mounted.mzML" in item.label for item in app.checkbox)
    (mount / "new-file.raw").write_bytes(b"new mounted acquisition")
    assert not any("new-file.raw" in item.label for item in app.checkbox)
    app.button(key="mounted_refresh_mzML-files").click().run(timeout=30)
    assert not app.exception
    assert any("new-file.raw" in item.label for item in app.checkbox)
    selected = next(item for item in app.checkbox if "mounted.mzML" in item.label)
    selected.check().run(timeout=30)
    selected_folder = next(item for item in app.checkbox if "bundle.d" in item.label)
    selected_folder.check().run(timeout=30)
    app.button(key="mounted_add_mzML-files").click().run(timeout=30)
    assert not app.exception
    references = tmp_path / "workspace/input-files/mzML-files/external_files.txt"
    assert set(references.read_text().splitlines()) == {
        str(mount / "mounted.mzML"),
        str(mount / "bundle.d"),
    }
    assert not (references.parent / "mounted.mzML").exists()
    assert not (references.parent / "bundle.d").exists()


def test_unmounted_configured_directory_does_not_expose_a_browser(
    picker_app, monkeypatch, tmp_path
):
    directory = tmp_path / "not-mounted"
    directory.mkdir()
    monkeypatch.setattr("src.workflow.StreamlitUI.os.path.ismount", lambda _path: False)
    app = picker_app(mount=directory)
    assert not any(button.key == "mounted_refresh_mzML-files" for button in app.button)
    assert app.text_input(key="path_to_folder_mzML-files")


def test_upload_page_explains_the_host_data_folder(naseweis_app, monkeypatch, tmp_path):  # noqa: F811
    inbox = tmp_path / "data"
    inbox.mkdir()
    monkeypatch.setenv("LOCAL_DATA_HOST_DIR", str(inbox))
    naseweis_app.run(timeout=30)
    naseweis_app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    assert not naseweis_app.exception
    instructions = next(
        info.value for info in naseweis_app.info if str(inbox) in info.value
    )
    assert "Refresh folder" in instructions
    assert "saved runs" in instructions
    assert f"](<{inbox.as_uri()}/>)" in instructions
    assert any(code.value == str(inbox) for code in naseweis_app.code)
    assert any("If the link does not open your file browser" in caption.value for caption in naseweis_app.caption)


def test_upload_page_keeps_relative_host_folder_copyable_without_a_container_uri(
    naseweis_app, monkeypatch  # noqa: F811
):
    monkeypatch.setenv("LOCAL_DATA_HOST_DIR", "./data")
    naseweis_app.run(timeout=30)
    naseweis_app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    assert not naseweis_app.exception
    instructions = next(info.value for info in naseweis_app.info if "Your data folder:" in info.value)
    assert "file:" not in instructions
    assert any(code.value == "./data" for code in naseweis_app.code)
    assert any("relative to the folder containing docker-compose.yml" in caption.value for caption in naseweis_app.caption)
