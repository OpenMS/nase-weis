"""Smoke-test every NASEWEIS page and the complete example-data setup."""

import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

# Import UI dependencies before tests that temporarily mock Streamlit at collection.
from src.common.common import page_setup  # noqa: F401


@pytest.fixture
def app(tmp_path):
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
