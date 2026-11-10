"""Contracts for optional, per-acquisition RNA feature extraction and export."""

import json
import shutil
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pandas as pd
import pytest

from src.naseweis_lfq import (
    FEATURE_FINDER_DEFAULTS,
    FEATURE_FINDER_INPUTS,
    FEATURE_FINDER_OUTPUTS,
    LFQ_FEATURES,
    LFQ_TABLE,
    LFQ_TARGET_COLUMNS,
    LFQ_TOOLS,
    TEXT_EXPORTER_DEFAULTS,
    TEXT_EXPORTER_INPUTS,
    TEXT_EXPORTER_OUTPUTS,
    prepare_lfq,
    read_lfq_table,
    run_lfq,
)
from src.naseweis_tools import tool_schema


@pytest.fixture
def manager(tmp_path):
    directory = tmp_path / "ini"
    directory.mkdir()
    for tool in LFQ_TOOLS:
        shutil.copyfile(
            Path("tests/data/naseweis/nightly") / f"{tool}.ini",
            directory / f"{tool}.ini",
        )
    return SimpleNamespace(ini_dir=directory, create_ini=Mock(return_value=True))


@pytest.mark.parametrize(
    "tool,inputs,outputs",
    [
        (LFQ_TOOLS[0], FEATURE_FINDER_INPUTS, FEATURE_FINDER_OUTPUTS),
        (LFQ_TOOLS[1], TEXT_EXPORTER_INPUTS, TEXT_EXPORTER_OUTPUTS),
    ],
)
def test_pinned_native_file_options_have_workflow_owners(manager, tool, inputs, outputs):
    schema = tool_schema(manager.ini_dir / f"{tool}.ini")
    assert {key for key, entry in schema.items() if entry["type"] == "input-file"} == inputs
    assert {key for key, entry in schema.items() if entry["type"] == "output-file"} == outputs


@pytest.mark.parametrize("extension", ["mzML", "raw", "d", "MZML"])
def test_herbert_defaults_and_actual_native_flags_are_recorded(manager, extension):
    params = {"other": {"setting": "preserved"}}
    assert prepare_lfq(params, manager, f"spectra.{extension}") is None
    finder = params[LFQ_TOOLS[0]]
    assert finder == FEATURE_FINDER_DEFAULTS | {"EMGScoring:init_mom": "false"}
    assert params[LFQ_TOOLS[1]] == TEXT_EXPORTER_DEFAULTS
    assert "EMGScoring:init_mom" not in params["_flag_params"][LFQ_TOOLS[0]]
    assert "model:each_trace" in params["_flag_params"][LFQ_TOOLS[0]]
    assert "no_ids" in params["_flag_params"][LFQ_TOOLS[1]]
    assert params["other"] == {"setting": "preserved"}


def test_imported_paths_and_export_layout_cannot_override_the_workflow(manager):
    params = {
        "_defaults": {
            LFQ_TOOLS[0]: {"in": "/outside/spectra", "lib_out": "/outside/library", "log": "/outside/default-finder.log", "extract:mz_window": 25},
            LFQ_TOOLS[1]: {"out": "/outside/export", "log": "/outside/default-exporter.log", "feature:add_metavalues": -1, "out_type": "csv"},
            "OtherTool": {"in": "keep-this-setting"},
        },
        LFQ_TOOLS[0]: {"id": "/outside/targets", "out": "/outside/features", "log": "/outside/override-finder.log", "extract:mz_window": "12.5", "model:each_trace": "true"},
        LFQ_TOOLS[1]: {"in": "/outside/features", "consensus:elements": "/outside/elements", "log": "/outside/override-exporter.log", "quoting": "double", "no_ids": False, "feature:minimal": True, "debug": 2},
    }
    assert prepare_lfq(params, manager, "spectra.mzML") is None
    assert params[LFQ_TOOLS[0]]["extract:mz_window"] == 12.5
    assert params[LFQ_TOOLS[0]]["model:each_trace"] is True
    assert params[LFQ_TOOLS[1]] == TEXT_EXPORTER_DEFAULTS | {"debug": 2}
    for tool, owned in [
        (LFQ_TOOLS[0], FEATURE_FINDER_INPUTS | FEATURE_FINDER_OUTPUTS | {"log"}),
        (LFQ_TOOLS[1], TEXT_EXPORTER_INPUTS | TEXT_EXPORTER_OUTPUTS | {"log"}),
    ]:
        assert not owned & params[tool].keys()
        assert not owned & params["_defaults"][tool].keys()
    assert params["_defaults"][LFQ_TOOLS[1]] == TEXT_EXPORTER_DEFAULTS
    assert params["_defaults"]["OtherTool"] == {"in": "keep-this-setting"}


@pytest.mark.parametrize("tool", LFQ_TOOLS)
@pytest.mark.parametrize("layer", ["defaults", "overrides"])
def test_imported_native_log_paths_are_removed_from_each_layer(manager, tool, layer):
    settings = {"log": "/outside/imported.log"}
    params = {"_defaults": {tool: settings}} if layer == "defaults" else {tool: settings}
    params.setdefault("_defaults", {}).setdefault("FeatureFinderMetaboIdent", {})[
        "extract:mz_window"
    ] = 15.0
    params.setdefault("FeatureFinderMetaboIdent", {})["detect:peak_width"] = 45.0

    assert prepare_lfq(params, manager, "spectra.mzML") is None

    assert "log" not in params[tool]
    assert "log" not in params["_defaults"][tool]
    assert params["FeatureFinderMetaboIdent"]["extract:mz_window"] == 15.0
    assert params["FeatureFinderMetaboIdent"]["detect:peak_width"] == 45.0


def test_new_nightly_file_paths_are_also_sanitized(manager):
    path = manager.ini_dir / f"{LFQ_TOOLS[0]}.ini"
    path.write_text(path.read_text().replace('<NODE name="extract"', '<ITEM name="future_out" value="" type="output-file"/><NODE name="extract"'))
    params = {LFQ_TOOLS[0]: {"future_out": "bad"}, "_defaults": {LFQ_TOOLS[0]: {"future_out": "bad"}}}
    assert prepare_lfq(params, manager, "spectra.mzML") is None
    assert "future_out" not in params[LFQ_TOOLS[0]]
    assert "future_out" not in params["_defaults"][LFQ_TOOLS[0]]


@pytest.mark.parametrize(
    "changes,reason",
    [
        ({"extract:mz_window": -1}, "within 0.0:"),
        ({"extract:n_isotopes": 1}, "within 2:"),
        ({"extract:n_isotopes": 2.5}, "integer"),
        ({"extract:isotope_pmin": 1.1}, "within 0.0:1.0"),
        ({"extract:isotope_pmin": float("nan")}, "finite number"),
        ({"detect:peak_width": float("inf")}, "finite number"),
        ({"detect:signal_to_noise": "bad"}, "finite number"),
        ({"detect:signal_to_noise": True}, "finite number"),
        ({"detect:signal_to_noise": None}, "finite number"),
        ({"model:type": "unknown"}, "one of"),
        ({"model:type": None}, "must be text"),
        ({"model:each_trace": "yes"}, "true or false"),
        ({"extract:rt_window": 0, "detect:peak_width": 0}, "positive"),
        ({"obsolete:option": 1}, "Unsupported parameters"),
    ],
)
def test_invalid_scientific_settings_do_not_mutate_snapshot(manager, changes, reason):
    params = {LFQ_TOOLS[0]: changes, "unrelated": {"keep": 1}}
    before = json.dumps(params)
    assert reason in prepare_lfq(params, manager, "spectra.mzML")
    assert json.dumps(params) == before


@pytest.mark.parametrize("params", [
    {"_defaults": []}, {"_flag_params": []},
    {LFQ_TOOLS[0]: []}, {"_defaults": {LFQ_TOOLS[1]: []}},
])
def test_invalid_parameter_structure_is_actionable(manager, params):
    assert "Invalid stored" in prepare_lfq(params, manager, "spectra.mzML")


@pytest.mark.parametrize("spectra", [None, [], "spectra.mgf", "spectra"])
def test_unsupported_acquisition_type_fails_before_reading_tool_schemas(manager, spectra):
    assert "acquisition" in prepare_lfq({}, manager, spectra)
    manager.create_ini.assert_not_called()


@pytest.mark.parametrize("tool", LFQ_TOOLS)
def test_missing_tool_is_an_actionable_failure(manager, tool):
    manager.create_ini.side_effect = lambda name: name != tool
    params = {}
    assert f"Cannot read {tool} parameters" in prepare_lfq(params, manager, "spectra.mzML")
    assert params == {}


@pytest.mark.parametrize("change,reason", [("broken", "schema"), ("missing", "schema"), ("no_io", "required LFQ")])
def test_unreadable_or_incompatible_schema_cannot_start_a_run(manager, change, reason):
    path = manager.ini_dir / f"{LFQ_TOOLS[0]}.ini"
    if change == "broken":
        path.write_text("<broken")
    elif change == "missing":
        path.unlink()
    else:
        path.write_text(path.read_text().replace('name="id"', 'name="missing_id"'))
    assert reason in prepare_lfq({}, manager, "spectra.mzML")


@pytest.mark.parametrize("value,valid", [([1, "2.5"], True), ("1,2", False), ([-1], False)])
def test_new_native_numeric_lists_use_the_schema(manager, value, valid):
    path = manager.ini_dir / f"{LFQ_TOOLS[0]}.ini"
    path.write_text(path.read_text().replace('<NODE name="extract"', '<ITEMLIST name="future_values" type="double" restrictions="0.0:"/><NODE name="extract"'))
    params = {LFQ_TOOLS[0]: {"future_values": value}}
    error = prepare_lfq(params, manager, "spectra.mzML")
    assert (error is None) == valid
    if valid:
        assert params[LFQ_TOOLS[0]]["future_values"] == [1.0, 2.5]


@pytest.fixture
def workflow(tmp_path):
    targets = tmp_path / "targets.tsv"
    targets.write_text("\t".join(LFQ_TARGET_COLUMNS) + "\nACGU\tC10H10O2\t0\t-2\t60\t0\t0\n")
    return SimpleNamespace(
        results_dir=tmp_path / "run-results",
        params={"mzML-files": "/inputs/original.mzML"},
        logger=Mock(), executor=Mock(run_topp=Mock(return_value=True)),
        targets=targets,
    )


def test_feature_extraction_and_export_receive_only_run_owned_paths(workflow):
    assert run_lfq(workflow, workflow.targets) is True
    calls = workflow.executor.run_topp.call_args_list
    assert [call.args[0] for call in calls] == list(LFQ_TOOLS)
    assert calls[0].args[1] == {
        "in": [workflow.params["mzML-files"]], "id": [str(workflow.targets)],
        "out": [str(workflow.results_dir / LFQ_FEATURES)],
    }
    assert calls[1].args[1] == {"in": [str(workflow.results_dir / LFQ_FEATURES)], "out": [str(workflow.results_dir / LFQ_TABLE)]}


@pytest.mark.parametrize("outcomes,count", [([False], 1), ([True, False], 2)])
def test_native_failure_stops_the_lfq_chain(workflow, outcomes, count):
    workflow.executor.run_topp.side_effect = outcomes
    assert run_lfq(workflow, workflow.targets) is False
    assert workflow.executor.run_topp.call_count == count


@pytest.mark.parametrize("contents", [None, "", "not\ta\ttarget\n", "\t".join(LFQ_TARGET_COLUMNS) + "\na\tb\n", '"broken'])
def test_missing_or_malformed_native_targets_are_reported(workflow, contents):
    if contents is None:
        workflow.targets.unlink()
    else:
        workflow.targets.write_text(contents)
    assert run_lfq(workflow, workflow.targets) is False
    workflow.executor.run_topp.assert_not_called()
    assert "ERROR" in workflow.logger.log.call_args.args[0]


def test_zero_identified_targets_make_a_native_empty_feature_map(workflow):
    import pyopenms as poms

    workflow.targets.write_text("\t".join(LFQ_TARGET_COLUMNS) + "\n\n")
    assert run_lfq(workflow, workflow.targets) is True
    workflow.executor.run_topp.assert_called_once_with(
        "TextExporter", {"in": [str(workflow.results_dir / LFQ_FEATURES)], "out": [str(workflow.results_dir / LFQ_TABLE)]}
    )
    features = poms.FeatureMap()
    poms.FeatureXMLFile().load(str(workflow.results_dir / LFQ_FEATURES), features)
    assert features.size() == 0


def test_native_measurements_are_numeric_but_labels_and_meta_are_preserved(tmp_path):
    path = tmp_path / "features.tsv"
    path.write_text("# arbitrary exporter comment\n#rt\tmz\tintensity\tcharge\texpected_rt\tlabel\tidentifier\tsum_formula\tnative_meta\n60\t500.25\t1.2e6\t-2\t60\t[m5C]ACGU\t00012\tC10H12O3\t[1,2,3]\n\n61\t501.25\t0\t-3\tNA\tACGU\t00013\tC20H22O3\tunknown\n")
    table = read_lfq_table(path)
    assert table["intensity"].tolist() == [1200000.0, 0]
    assert table["charge"].tolist() == [-2, -3]
    assert pd.isna(table.loc[1, "expected_rt"])
    assert table["label"].tolist() == ["[m5C]ACGU", "ACGU"]
    assert table["identifier"].tolist() == ["00012", "00013"]
    assert table["native_meta"].tolist() == ["[1,2,3]", "unknown"]


@pytest.mark.parametrize("contents", ["#rt\tmz\tintensity\tlabel\n", "# comments\n#rt\tmz\tintensity\n"])
def test_empty_native_exports_remain_empty_tables(tmp_path, contents):
    path = tmp_path / "empty.tsv"
    path.write_text(contents)
    table = read_lfq_table(path)
    assert table.empty
    assert {"rt", "mz", "intensity"} <= set(table)


def test_typed_native_export_records_are_supported(tmp_path):
    path = tmp_path / "features.tsv"
    path.write_text("#RUN\tidentifier\n#FEATURE\trt\tmz\tintensity\tlabel\nRUN\trun1\nPROTEIN\tignored\nFEATURE\t60\t500\t10\tACGU\nPEPTIDE\tignored\n")
    table = read_lfq_table(path)
    assert list(table) == ["rt", "mz", "intensity", "label"]
    assert table.loc[0, "label"] == "ACGU"


@pytest.mark.parametrize("contents,reason", [
    ("", "Missing"),
    ("# comments only\n", "Missing"),
    ("1\t2\t3\n", "Missing"),
    ("#rt\tbad\n", "Invalid feature header"),
    ("#rt\tmz\tintensity\tmz\n", "duplicate"),
    ("#rt\tmz\tintensity\n#rt\tmz\tintensity\n", "duplicate"),
    ("#rt\tmz\tintensity\n1\t2\n", "columns"),
    ("#rt\tmz\tintensity\n1\t2\tbad\n", "numeric"),
    ('#rt\tmz\tintensity\n1\t2\t"broken', "Invalid TextExporter"),
])
def test_corrupt_exports_raise_actionable_errors(tmp_path, contents, reason):
    path = tmp_path / "invalid.tsv"
    path.write_text(contents)
    with pytest.raises(ValueError, match=reason):
        read_lfq_table(path)
