"""Only execution can capture trusted provenance for its selected acquisition."""

import os
from copy import deepcopy
from pathlib import Path
from unittest.mock import Mock

import pytest
from test_naseweis_workflow import workflow  # noqa: F401

from src.naseweis_provenance import acquisition_fingerprint


def test_execution_replaces_imported_record_with_selected_acquisition(workflow):  # noqa: F811
    selected = Path(workflow.params["mzML-files"])
    selected.write_bytes(b"selected acquisition, no XML parse required")
    workflow.params.update(
        {
            "_acquisition_fingerprint": {
                "path": "/wrong/file.mzML",
                "sha256": "forged",
            },
            "NucleicAcidSearchEngine": {"scoring:method": "mvh"},
        }
    )
    original = workflow.params
    before = deepcopy(original)
    assert workflow.execution() is True
    assert original == before
    saved = workflow.parameter_manager.read_parameters_strict()
    assert saved["_acquisition_fingerprint"] == acquisition_fingerprint(selected)
    assert saved["NucleicAcidSearchEngine"]["scoring:method"] == "mvh"
    assert workflow.executor.run_topp.call_args.args[1]["in"] == [str(selected)]


@pytest.mark.parametrize("suffix", [".raw", ".d"])
def test_vendor_execution_clears_stale_record_without_hashing(
    workflow,  # noqa: F811
    monkeypatch,
    suffix,
):
    selected = workflow.workflow_dir / f"vendor{suffix}"
    if suffix == ".d":
        selected.mkdir()
    else:
        selected.touch()
    workflow.params.update(
        {"mzML-files": str(selected), "_acquisition_fingerprint": {"sha256": "stale"}}
    )
    original = workflow.params
    before = deepcopy(original)
    hasher = Mock(
        side_effect=AssertionError("Vendor acquisitions must not be fingerprinted.")
    )
    monkeypatch.setattr("src.naseweis_provenance.acquisition_fingerprint", hasher)
    assert workflow.execution() is True
    assert (
        "_acquisition_fingerprint"
        not in workflow.parameter_manager.read_parameters_strict()
    )
    assert "_acquisition_fingerprint" not in workflow.params
    assert original == before
    hasher.assert_not_called()


def test_failed_fingerprint_does_not_save_snapshot_or_run_tools(workflow, monkeypatch):  # noqa: F811
    before = workflow.parameter_manager.params_file.read_bytes()
    monkeypatch.setattr(
        "src.naseweis_provenance.acquisition_fingerprint",
        Mock(side_effect=ValueError("The acquisition changed; rerun the search.")),
    )
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert workflow.parameter_manager.params_file.read_bytes() == before
    assert "changed; rerun" in (workflow.workflow_dir / "logs/minimal.log").read_text()


@pytest.mark.parametrize("lfq_enabled", [True, False])
def test_changed_post_search_input_fails_before_quantification_or_success(
    workflow,  # noqa: F811
    monkeypatch,
    lfq_enabled,
):
    selected = Path(workflow.params["mzML-files"])
    workflow.params.update(
        {"enable-lfq": lfq_enabled, "add-decoys": False, "estimate-fdr": False}
    )
    monkeypatch.setattr("src.Workflow.prepare_lfq", Mock(return_value=None))
    monkeypatch.setattr("src.Workflow.validate_lfq_spectra", Mock(return_value=None))
    quantification = Mock(return_value=True)
    monkeypatch.setattr("src.Workflow.run_lfq", quantification)

    def search(tool, _io):
        assert tool == "NucleicAcidSearchEngine"
        selected.write_bytes(b"changed during search")
        return True

    workflow.executor.run_topp.side_effect = search
    assert workflow.execution() is False
    assert workflow.executor.run_topp.call_count == 1
    quantification.assert_not_called()
    assert (
        "acquisition changed during the search"
        in (workflow.workflow_dir / "logs/minimal.log").read_text()
    )


def test_successful_search_rechecks_input_before_lfq(workflow, monkeypatch):  # noqa: F811
    workflow.params.update(
        {"enable-lfq": True, "add-decoys": False, "estimate-fdr": False}
    )
    monkeypatch.setattr("src.Workflow.prepare_lfq", Mock(return_value=None))
    monkeypatch.setattr("src.Workflow.validate_lfq_spectra", Mock(return_value=None))
    quantification = Mock(return_value=True)
    monkeypatch.setattr("src.Workflow.run_lfq", quantification)
    assert workflow.execution() is True
    quantification.assert_called_once()
    saved = workflow.parameter_manager.read_parameters_strict()
    assert saved["_acquisition_fingerprint"] == acquisition_fingerprint(
        workflow.params["mzML-files"]
    )


def test_changed_input_during_successful_lfq_cannot_complete(workflow, monkeypatch):  # noqa: F811
    workflow.params.update(
        {"enable-lfq": True, "add-decoys": False, "estimate-fdr": False}
    )
    monkeypatch.setattr("src.Workflow.prepare_lfq", Mock(return_value=None))
    monkeypatch.setattr("src.Workflow.validate_lfq_spectra", Mock(return_value=None))

    def replace_acquisition(_workflow, _targets):
        Path(workflow.params["mzML-files"]).write_bytes(
            b"replaced during quantification"
        )
        return True

    quantification = Mock(side_effect=replace_acquisition)
    monkeypatch.setattr("src.Workflow.run_lfq", quantification)
    assert workflow.execution() is False
    quantification.assert_called_once()
    assert (
        "acquisition changed during the search"
        in (workflow.workflow_dir / "logs/minimal.log").read_text()
    )


@pytest.mark.parametrize("stage", ["search", "lfq"])
@pytest.mark.parametrize("replacement", ["overwrite", "replace"])
def test_same_size_replacement_with_restored_mtime_cannot_complete(
    workflow,  # noqa: F811
    monkeypatch,
    stage,
    replacement,
):
    selected = Path(workflow.params["mzML-files"])
    selected.write_bytes(b"original")
    workflow.params.update(
        {"enable-lfq": True, "add-decoys": False, "estimate-fdr": False}
    )
    monkeypatch.setattr("src.Workflow.prepare_lfq", Mock(return_value=None))
    monkeypatch.setattr("src.Workflow.validate_lfq_spectra", Mock(return_value=None))

    def replace_acquisition():
        record = workflow.params["_acquisition_fingerprint"]
        if replacement == "overwrite":
            selected.write_bytes(b"replaced")
        else:
            other = selected.with_name("replacement.mzML")
            other.write_bytes(b"replaced")
            other.replace(selected)
        os.utime(selected, ns=(record["mtime_ns"], record["mtime_ns"]))
        assert selected.stat().st_size == record["size"]
        assert selected.stat().st_mtime_ns == record["mtime_ns"]

    def search(_tool, _io):
        if stage == "search":
            replace_acquisition()
        return True

    def quantify(_workflow, _targets):
        if stage == "lfq":
            replace_acquisition()
        return True

    quantification = Mock(side_effect=quantify)
    workflow.executor.run_topp.side_effect = search
    monkeypatch.setattr("src.Workflow.run_lfq", quantification)
    assert workflow.execution() is False
    assert quantification.call_count == int(stage == "lfq")
    assert (
        "acquisition changed during the search"
        in (workflow.workflow_dir / "logs/minimal.log").read_text()
    )
