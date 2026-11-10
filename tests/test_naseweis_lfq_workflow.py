"""Optional LFQ dispatch, parameter ownership and fail-fast contracts."""

from copy import deepcopy
from pathlib import Path
from unittest.mock import Mock

import pytest
from test_naseweis_workflow import workflow  # noqa: F401

from src.naseweis_tools import tool_schema

SCHEMAS = Path(__file__).parent / "data/naseweis/nightly"
LFQ_TOOLS = ["FeatureFinderMetaboIdent", "TextExporter"]
FEATURES = "quantification/features.featureXML"
TABLE = "quantification/features.tsv"
TARGETS = "exports/lfq-targets.tsv"


@pytest.fixture
def lfq_workflow(workflow, monkeypatch):  # noqa: F811
    workflow.params.update(
        {"add-decoys": False, "estimate-fdr": False, "enable-lfq": True}
    )
    for tool in LFQ_TOOLS:
        (workflow.parameter_manager.ini_dir / f"{tool}.ini").write_bytes(
            (SCHEMAS / f"{tool}.ini").read_bytes()
        )
    # Input content validation has separate real mzML tests. Here the empty
    # fixture keeps only workflow dispatch and failure boundaries under test.
    monkeypatch.setattr("src.Workflow.validate_lfq_spectra", Mock(return_value=None))
    workflow.executor.run_topp.side_effect = successful_tool
    return workflow


def write_targets(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "CompoundName\tSumFormula\tMass\tCharge\tRetentionTime\tRetentionTimeRange\tIsoDistribution\n"
        "rna\tC10H15N5O10P2\t0\t-2\t10\t0\t0\n"
    )


def successful_tool(tool, io):
    # The native target export is the boundary between search and extraction.
    # Keep its real header so malformed/missing-target safeguards stay active.
    if tool == "NucleicAcidSearchEngine" and "lfq_out" in io:
        write_targets(Path(io["lfq_out"][0]))
    for key in ["out", "id_out"]:
        if key in io:
            destination = Path(io[key][0])
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text("test output\n")
    return True


def calls(subject):
    return subject.executor.run_topp.call_args_list


def log(subject):
    return (subject.workflow_dir / "logs/minimal.log").read_text()


@pytest.mark.parametrize("setting", [None, False, "false"])
def test_lfq_is_disabled_by_default_and_requires_no_extra_tools(lfq_workflow, setting):
    if setting is None:
        lfq_workflow.params.pop("enable-lfq")
    else:
        lfq_workflow.params["enable-lfq"] = setting
    assert lfq_workflow.execution() is True
    assert [call.args[0] for call in calls(lfq_workflow)] == ["NucleicAcidSearchEngine"]
    assert "lfq_out" not in calls(lfq_workflow)[0].args[1]
    checked = {
        call.args[0]
        for call in lfq_workflow.parameter_manager.create_ini.call_args_list
    }
    assert not checked.intersection(LFQ_TOOLS)


def test_target_export_alone_does_not_enable_quantification(lfq_workflow):
    lfq_workflow.params.update({"enable-lfq": False, "search-exports": ["lfq_out"]})
    assert lfq_workflow.execution() is True
    assert [call.args[0] for call in calls(lfq_workflow)] == ["NucleicAcidSearchEngine"]
    assert calls(lfq_workflow)[0].args[1]["lfq_out"] == [
        str(lfq_workflow.results_dir / TARGETS)
    ]


@pytest.mark.parametrize("exports", [[], ["id_out"]])
def test_lfq_uses_original_spectra_and_adds_target_export_without_changing_selection(
    lfq_workflow, exports
):
    lfq_workflow.params["search-exports"] = exports
    original = lfq_workflow.params
    before = deepcopy(original)
    assert lfq_workflow.execution() is True
    assert original == before
    assert [call.args[0] for call in calls(lfq_workflow)] == [
        "NucleicAcidSearchEngine",
        *LFQ_TOOLS,
    ]
    nase, features, table = [call.args[1] for call in calls(lfq_workflow)]
    assert nase["lfq_out"] == [str(lfq_workflow.results_dir / TARGETS)]
    assert features == {
        "in": [before["mzML-files"]],
        "id": [str(lfq_workflow.results_dir / TARGETS)],
        "out": [str(lfq_workflow.results_dir / FEATURES)],
    }
    assert table == {
        "in": [str(lfq_workflow.results_dir / FEATURES)],
        "out": [str(lfq_workflow.results_dir / TABLE)],
    }
    saved = lfq_workflow.parameter_manager.read_parameters_strict()
    assert saved["search-exports"] == exports
    assert saved["enable-lfq"] is True


@pytest.mark.parametrize("suffix", [".raw", ".d"])
def test_lfq_pipeline_passes_original_vendor_acquisition(lfq_workflow, suffix):
    acquisition = lfq_workflow.workflow_dir / f"acquisition{suffix}"
    if suffix == ".d":
        acquisition.mkdir()
    else:
        acquisition.touch()
    lfq_workflow.params["mzML-files"] = str(acquisition)
    assert lfq_workflow.execution() is True
    nase, features, _ = [call.args[1] for call in calls(lfq_workflow)]
    assert nase["in"] == [str(acquisition)]
    assert features["in"] == [str(acquisition)]


@pytest.mark.parametrize("missing", LFQ_TOOLS)
def test_missing_lfq_tool_is_reported_before_search(lfq_workflow, monkeypatch, missing):
    monkeypatch.setattr(
        "src.Workflow.shutil.which",
        lambda tool: None if tool == missing else f"/tools/{tool}",
    )
    before = lfq_workflow.parameter_manager.params_file.read_bytes()
    assert lfq_workflow.execution() is False
    lfq_workflow.executor.run_topp.assert_not_called()
    assert missing in log(lfq_workflow)
    assert lfq_workflow.parameter_manager.params_file.read_bytes() == before


@pytest.mark.parametrize("missing", LFQ_TOOLS)
def test_missing_lfq_schema_is_reported_before_search(lfq_workflow, missing):
    lfq_workflow.parameter_manager.create_ini.side_effect = lambda tool: tool != missing
    before = lfq_workflow.parameter_manager.params_file.read_bytes()
    assert lfq_workflow.execution() is False
    lfq_workflow.executor.run_topp.assert_not_called()
    assert missing in log(lfq_workflow)
    assert lfq_workflow.parameter_manager.params_file.read_bytes() == before


@pytest.mark.parametrize(
    "tool,options,option",
    [
        ("FeatureFinderMetaboIdent", {"unknown_lfq_setting": 2}, "unknown_lfq_setting"),
        ("TextExporter", {"unknown_lfq_setting": 2}, "unknown_lfq_setting"),
        ("FeatureFinderMetaboIdent", {"extract:n_isotopes": 1}, "extract:n_isotopes"),
        ("FeatureFinderMetaboIdent", {"extract:mz_window": -1.0}, "extract:mz_window"),
        ("FeatureFinderMetaboIdent", {"model:type": "invalid"}, "model:type"),
        ("TextExporter", {"id:add_metavalues": 101}, "id:add_metavalues"),
    ],
)
def test_invalid_lfq_settings_fail_before_search_and_do_not_rewrite_inputs(
    lfq_workflow, tool, options, option
):
    lfq_workflow.params[tool] = options
    original = deepcopy(lfq_workflow.params)
    before = lfq_workflow.parameter_manager.params_file.read_bytes()
    assert lfq_workflow.execution() is False
    lfq_workflow.executor.run_topp.assert_not_called()
    assert option in log(lfq_workflow)
    assert lfq_workflow.params == original
    assert lfq_workflow.parameter_manager.params_file.read_bytes() == before


def test_invalid_native_lfq_default_is_rejected_before_search(lfq_workflow):
    lfq_workflow.params["_defaults"] = {
        "FeatureFinderMetaboIdent": {"extract:n_isotopes": 1}
    }
    before = lfq_workflow.parameter_manager.params_file.read_bytes()
    assert lfq_workflow.execution() is False
    lfq_workflow.executor.run_topp.assert_not_called()
    assert "extract:n_isotopes" in log(lfq_workflow)
    assert lfq_workflow.parameter_manager.params_file.read_bytes() == before


def test_unusable_spectra_stop_lfq_before_any_tool_runs(lfq_workflow, monkeypatch):
    monkeypatch.setattr(
        "src.Workflow.validate_lfq_spectra", lambda path: "LFQ needs MS1 spectra"
    )
    before = lfq_workflow.parameter_manager.params_file.read_bytes()
    assert lfq_workflow.execution() is False
    lfq_workflow.executor.run_topp.assert_not_called()
    assert "MS1 spectra" in log(lfq_workflow)
    assert lfq_workflow.parameter_manager.params_file.read_bytes() == before


@pytest.mark.parametrize("failed_stage", [0, 1, 2])
def test_lfq_pipeline_stops_at_first_failed_tool(lfq_workflow, failed_stage):
    stage = 0

    def execute(tool, io):
        nonlocal stage
        current = stage
        stage += 1
        return successful_tool(tool, io) if current < failed_stage else False

    lfq_workflow.executor.run_topp.side_effect = execute
    assert lfq_workflow.execution() is False
    assert [call.args[0] for call in calls(lfq_workflow)] == [
        "NucleicAcidSearchEngine",
        *LFQ_TOOLS,
    ][: failed_stage + 1]


def test_lfq_file_options_in_both_layers_cannot_redirect_workspace_outputs(
    lfq_workflow,
):
    from src.workflow.CommandExecutor import CommandExecutor

    keys = {
        tool: {
            key
            for key, entry in tool_schema(SCHEMAS / f"{tool}.ini").items()
            if entry["type"] in {"input-file", "output-file"}
        }
        for tool in LFQ_TOOLS
    }
    lfq_workflow.params["_defaults"] = {
        tool: {key: f"/imported-default/{key}.wrong" for key in files}
        for tool, files in keys.items()
    }
    for tool, files in keys.items():
        lfq_workflow.params[tool] = {
            key: f"/imported-override/{key}.wrong" for key in files
        }
    lfq_workflow.params["_defaults"]["FeatureFinderMetaboIdent"][
        "extract:mz_window"
    ] = 12.0
    lfq_workflow.params["FeatureFinderMetaboIdent"]["detect:peak_width"] = 45.0
    lfq_workflow.params["_defaults"]["TextExporter"]["id:annotations"] = "default"
    lfq_workflow.params["TextExporter"].update(
        {
            "quoting": "double",
            "out_type": "csv",
            "no_ids": False,
            "feature:minimal": True,
            "feature:add_metavalues": 101,
        }
    )
    original = lfq_workflow.params
    before = deepcopy(original)
    executor = CommandExecutor(
        lfq_workflow.workflow_dir, lfq_workflow.logger, lfq_workflow.parameter_manager
    )
    executor._max_threads = 1
    write_targets(lfq_workflow.results_dir / TARGETS)
    executor.run_command = Mock(return_value=True)
    lfq_workflow.executor = executor
    assert lfq_workflow.execution() is True
    commands = {
        call.args[0][0]: call.args[0] for call in executor.run_command.call_args_list
    }
    for tool in LFQ_TOOLS:
        command = commands[tool]
        assert not any(argument.startswith("/imported-") for argument in command)
        active = (
            {"in", "out", "id"} if tool == "FeatureFinderMetaboIdent" else {"in", "out"}
        )
        for key in keys[tool]:
            assert command.count(f"-{key}") == int(key in active)
    features = commands["FeatureFinderMetaboIdent"]
    table = commands["TextExporter"]
    assert features[features.index("-extract:mz_window") + 1] == "12.0"
    assert features[features.index("-detect:peak_width") + 1] == "45.0"
    assert table[table.index("-id:annotations") + 1] == "default"
    assert table[table.index("-quoting") + 1] == "none"
    assert table[table.index("-out_type") + 1] == "tsv"
    assert table[table.index("-feature:add_metavalues") + 1] == "0"
    assert "-no_ids" in table
    assert "-feature:minimal" not in table
    saved = lfq_workflow.parameter_manager.read_parameters_strict()
    assert original == before
    for tool in LFQ_TOOLS:
        assert not keys[tool].intersection(saved[tool])
        assert not keys[tool].intersection(saved["_defaults"][tool])


def store_scan_levels(path, levels):
    import pyopenms as oms

    experiment = oms.MSExperiment()
    for level in levels:
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(level)
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)


def test_real_ms2_only_mzml_is_rejected_for_lfq(tmp_path):
    from src.Workflow import validate_lfq_spectra

    acquisition = tmp_path / "ms2-only.mzML"
    store_scan_levels(acquisition, [2, 2])
    assert "MS1" in validate_lfq_spectra(str(acquisition))


def test_real_ms1_and_ms2_metadata_passes_lfq_input_check(tmp_path):
    from src.Workflow import validate_lfq_spectra

    acquisition = tmp_path / "ms1-and-ms2.mzML"
    store_scan_levels(acquisition, [1, 2])
    assert validate_lfq_spectra(str(acquisition)) is None


def test_malformed_mzml_gives_actionable_lfq_input_error(tmp_path):
    from src.Workflow import validate_lfq_spectra

    acquisition = tmp_path / "malformed.mzML"
    acquisition.write_text("<mzML><not-valid")
    error = validate_lfq_spectra(str(acquisition))
    assert error and "Cannot inspect" in error


@pytest.mark.parametrize("suffix", [".raw", ".d"])
def test_vendor_lfq_input_defers_to_native_reader(tmp_path, suffix):
    from src.Workflow import validate_lfq_spectra

    acquisition = tmp_path / f"vendor{suffix}"
    if suffix == ".d":
        acquisition.mkdir()
    else:
        acquisition.write_bytes(b"vendor data handled by native reader")
    assert validate_lfq_spectra(str(acquisition)) is None
