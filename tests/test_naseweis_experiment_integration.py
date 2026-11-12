"""Experiment integration preserves sample design and native viewer evidence."""

import json
from pathlib import Path

import pandas as pd
import pytest
from test_naseweis_lfq_workflow import lfq_workflow  # noqa: F401
from test_naseweis_workflow import workflow  # noqa: F401


@pytest.mark.parametrize("selected", [None, "", "   ", [], [None], ["/data/one.mzML", None], {"bad": "selection"}])
def test_invalid_selection_cannot_create_a_run_or_rewrite_the_design(workflow, monkeypatch, selected):  # noqa: F811
    import streamlit as st

    from src.workflow.ExperimentalDesign import design_path
    from src.workflow.RunManager import RunManager

    workflow.params["mzML-files"] = selected
    workflow.parameter_manager.write_parameters(workflow.params)
    before = workflow.parameter_manager.params_file.read_bytes()
    workflow.runs = RunManager(workflow.workflow_dir)
    monkeypatch.setattr(st, "session_state", {})
    with pytest.raises(ValueError, match="Select.*Configure"):
        workflow._new_run()
    assert workflow.runs.list_runs() == []
    assert not design_path(workflow.workflow_dir).exists()
    assert workflow.parameter_manager.params_file.read_bytes() == before
    workflow.executor.run_topp.assert_not_called()


def test_start_snapshots_the_design_in_selected_file_order(workflow, monkeypatch):  # noqa: F811
    import streamlit as st

    from src.workflow.ExperimentalDesign import save_design
    from src.workflow.RunManager import RunManager

    first = Path(workflow.params["mzML-files"])
    second = first.with_name("second.mzML")
    second.touch()
    selected = [str(second), str(first)]
    workflow.params["mzML-files"] = selected
    workflow.parameter_manager.write_parameters(workflow.params)
    design = pd.DataFrame({
        "file": [first.name, second.name],
        "condition": ["control", "treated"],
        "replicate": [2, 1],
        "fraction": [1, 1],
    })
    save_design(workflow.workflow_dir, design)
    workflow.runs = RunManager(workflow.workflow_dir)
    monkeypatch.setattr(st, "session_state", {})
    workflow._new_run()
    snapshot = workflow.runs.params_file(workflow.run_id)
    saved = json.loads(snapshot.read_text())
    assert saved["mzML-files"] == selected
    assert saved["_sample_design"] == [
        {"file": str(second), "condition": "treated", "replicate": 1, "fraction": 1},
        {"file": str(first), "condition": "control", "replicate": 2, "fraction": 1},
    ]
    before = snapshot.read_bytes()
    design["condition"] = "changed after starting"
    save_design(workflow.workflow_dir, design)
    assert snapshot.read_bytes() == before


@pytest.mark.parametrize("settings", [{}, None, [], "corrupt"])
def test_fitted_preset_is_explicit_and_does_not_enable_lfq(workflow, monkeypatch, settings):  # noqa: F811
    import streamlit as st

    workflow.params["enable-lfq"] = False
    workflow.params["FeatureFinderMetaboIdent"] = settings
    workflow.parameter_manager.write_parameters(workflow.params)
    monkeypatch.setattr(st, "session_state", {})
    workflow._apply_rna_comparison_preset()
    saved = workflow.parameter_manager.read_parameters_strict()
    assert saved["enable-lfq"] is False
    assert saved["lfq-measure"] == "model_area"
    assert saved["lfq-normalization"] == "total"
    assert saved["FeatureFinderMetaboIdent"] == {
        "model:type": "symmetric", "extract:mz_window": 10.0,
        "extract:rt_window": 0.0, "detect:peak_width": 60.0,
        "extract:n_isotopes": 5, "extract:isotope_pmin": 0.0,
        "EMGScoring:init_mom": False,
    }


def test_legacy_fitted_measure_does_not_fall_back_without_a_model(lfq_workflow):  # noqa: F811
    lfq_workflow.params["lfq-measure"] = "model_area"
    assert lfq_workflow.execution() is False
    lfq_workflow.executor.run_topp.assert_not_called()
    log = lfq_workflow.workflow_dir / "logs/minimal.log"
    assert "require symmetric or asymmetric" in log.read_text()


def test_raw_viewer_discovers_all_selected_external_acquisitions(tmp_path):
    from src.naseweis_viewers import input_acquisitions

    first, second = tmp_path / "one.mzML", tmp_path / "two.mzML"
    first.touch()
    second.touch()
    workflow_dir = tmp_path / "workspace/naseweis"
    workflow_dir.mkdir(parents=True)
    params = {"mzML-files": [str(first), str(second)]}
    assert input_acquisitions(workflow_dir, params) == [str(first), str(second)]
    assert input_acquisitions(workflow_dir, params, online=True) == []
