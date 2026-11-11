"""Batch isolation, immutable design provenance and complete-cohort comparison."""

import json
from copy import deepcopy
from types import SimpleNamespace
from unittest.mock import Mock

import pandas as pd
import pytest

from src.naseweis_experiment import (
    load_manifest,
    run_experiment,
    sample_results,
    sample_snapshot,
    validate_design_source,
)


@pytest.fixture
def experiment(tmp_path):
    from src.workflow.Logger import Logger
    from src.workflow.ParameterManager import ParameterManager

    directory = tmp_path / "workspace/naseweis"
    results = directory / "runs/cohort/results"
    results.mkdir(parents=True)
    inputs = directory / "input-files/mzML-files"
    inputs.mkdir(parents=True)
    files = [inputs / name for name in ["sample A.mzML", "sample B.mzML"]]
    for path in files:
        path.touch()
    params = {
        "mzML-files": [str(path) for path in files],
        "max_threads": 2,
        "NucleicAcidSearchEngine": {"scoring:method": "mvh"},
        "_sample_design": [
            {"file": str(path), "condition": name, "replicate": 1, "fraction": 1}
            for path, name in zip(files, ["Control", "Treatment"])
        ],
    }
    manager = ParameterManager(directory)
    manager.params_file = results.parent / "params.json"
    manager.write_parameters(params)
    return SimpleNamespace(
        workflow_dir=directory,
        results_dir=results,
        params=params,
        parameter_manager=manager,
        logger=Logger(directory),
        executor=Mock(),
    )


def feature_output(child):
    folder = child.results_dir / "quantification"
    folder.mkdir()
    area = 10 if "sample A" in child.params["mzML-files"] else 20
    (folder / "features.tsv").write_text(
        "#rt\tmz\tintensity\tcharge\tlabel\tsum_formula\tmodel_area\n"
        f"10\t500\t{area}\t-2\tACGU\tC38H49N15O26P3\t{area * 2}\n"
        "11\t600\t10\t-3\tUUUU\tC36H46N8O29P3\t40\n"
    )
    return True


def test_each_sample_uses_an_isolated_snapshot_dispatcher_results_and_logs(experiment):
    before, saved = (
        deepcopy(experiment.params),
        experiment.parameter_manager.params_file.read_bytes(),
    )
    children = []

    def execute(child):
        assert child is not experiment
        number = len(children) + 1
        identity = f"sample_{number:03d}"
        assert child.workflow_dir == experiment.workflow_dir
        assert child.results_dir == sample_results(experiment.results_dir, identity)
        assert child.file_manager.results_dir == child.results_dir
        assert child.parameter_manager.params_file == sample_snapshot(
            experiment.results_dir, identity
        )
        assert child.executor.pid_dir == experiment.workflow_dir / "pids"
        assert child.executor.pid_dir.is_dir()
        assert child.executor.parameter_manager is child.parameter_manager
        assert child.params["mzML-files"] == before["mzML-files"][number - 1]
        assert child.params["_sample_design"] == [before["_sample_design"][number - 1]]
        child.logger.log(f"Only {identity} ran this command.")
        child.params["NucleicAcidSearchEngine"]["scoring:method"] = "hyperscore"
        child.parameter_manager.write_parameters(child.params)
        children.append(child)
        return True

    assert run_experiment(experiment, execute) is True
    assert len(children) == 2
    assert experiment.params == before
    assert experiment.parameter_manager.params_file.read_bytes() == saved
    experiment.executor.run_topp.assert_not_called()
    manifest = load_manifest(experiment.results_dir)
    assert manifest["lfq_enabled"] is False
    assert manifest["comparison_status"] == "not_requested"
    assert manifest["status"] == "finished"
    assert [sample["status"] for sample in manifest["samples"]] == ["finished"] * 2
    assert not (experiment.results_dir / "comparison").exists()
    for number in [1, 2]:
        identity = f"sample_{number:03d}"
        local = (
            sample_snapshot(experiment.results_dir, identity).parent / "logs/all.log"
        )
        assert f"Only {identity}" in local.read_text()
        assert f"Only sample_{3 - number:03d}" not in local.read_text()
    assert "Only sample_001" in (experiment.workflow_dir / "logs/all.log").read_text()
    assert "Only sample_002" in (experiment.workflow_dir / "logs/all.log").read_text()


def test_identification_batch_defaults_to_lfq_off_and_allows_one_file_without_design(
    experiment,
):
    experiment.params.pop("_sample_design")
    experiment.params["mzML-files"] = experiment.params["mzML-files"][:1]
    child = Mock(return_value=True)
    assert run_experiment(experiment, child) is True
    child.assert_called_once()
    sample = load_manifest(experiment.results_dir)["samples"][0]
    assert (sample["condition"], sample["replicate"], sample["fraction"]) == ("", 1, 1)


@pytest.mark.parametrize(
    "measure,normalization",
    [("intensity", "total"), ("intensity", "none"), ("model_area", "total")],
)
def test_completed_lfq_cohort_generates_the_real_comparison(
    experiment, measure, normalization
):
    experiment.params.update(
        {"enable-lfq": True, "lfq-measure": measure, "lfq-normalization": normalization}
    )
    if measure == "model_area":
        experiment.params["FeatureFinderMetaboIdent"] = {"model:type": "symmetric"}
    assert run_experiment(experiment, feature_output) is True
    manifest = load_manifest(experiment.results_dir)
    assert manifest["comparison_status"] == "finished"
    raw = pd.read_csv(
        experiment.results_dir / "comparison/raw-abundances.tsv", sep="\t", index_col=0
    )
    expected = [20, 40] if measure == "model_area" else [10, 20]
    assert raw.loc["ACGU"].tolist() == expected
    relative = pd.read_csv(
        experiment.results_dir / "comparison/abundances.tsv", sep="\t", index_col=0
    )
    if normalization == "total":
        assert relative.sum().tolist() == pytest.approx([1, 1])
    else:
        pd.testing.assert_frame_equal(relative, raw)


def test_design_is_read_only_from_the_snapshot_and_order_follows_selection(experiment):
    # The mutable framework design must never be consulted by the worker.
    (experiment.workflow_dir / "experimental-design.tsv").write_text(
        "garbage mutable design"
    )
    experiment.params["_sample_design"].reverse()
    for row in experiment.params["_sample_design"]:
        row["file"] = row["file"].split("/")[-1]
    assert run_experiment(experiment, lambda child: True) is True
    assert [
        sample["condition"]
        for sample in load_manifest(experiment.results_dir)["samples"]
    ] == ["Control", "Treatment"]


def test_path_alias_is_supported_but_conflicting_design_identity_is_rejected(
    experiment,
):
    first = experiment.params["_sample_design"][0]
    first["path"] = first["file"]
    assert run_experiment(experiment, lambda child: True) is True
    # Use another output root, preserving already produced results.
    experiment.results_dir = experiment.results_dir.parent / "another-results"
    first["path"] = experiment.params["mzML-files"][1]
    callback = Mock()
    assert run_experiment(experiment, callback) is False
    callback.assert_not_called()


@pytest.mark.parametrize(
    "selection", [None, [], "a.mzML", [None], [""], ["missing.mzML"], [["nested.mzML"]]]
)
def test_invalid_acquisition_selection_fails_before_dispatch(experiment, selection):
    experiment.params["mzML-files"] = selection
    callback = Mock()
    assert run_experiment(experiment, callback) is False
    callback.assert_not_called()
    assert load_manifest(experiment.results_dir) is None


@pytest.mark.parametrize(
    "kind", ["unsupported", "directory", "duplicate", "symlink", "basename", "casefold"]
)
def test_unsafe_or_ambiguous_input_paths_are_rejected(experiment, tmp_path, kind):
    first = experiment.params["mzML-files"][0]
    if kind == "unsupported":
        path = tmp_path / "input.txt"
        path.touch()
        experiment.params["mzML-files"][1] = str(path)
    elif kind == "directory":
        path = tmp_path / "folder.mzML"
        path.mkdir()
        experiment.params["mzML-files"][1] = str(path)
    elif kind == "duplicate":
        experiment.params["mzML-files"][1] = first
    elif kind == "symlink":
        path = tmp_path / "linked.mzML"
        path.symlink_to(first)
        experiment.params["mzML-files"][1] = str(path)
    else:
        from pathlib import Path

        path = tmp_path / (
            Path(first).name if kind == "basename" else Path(first).name.upper()
        )
        path.touch()
        experiment.params["mzML-files"][1] = str(path)
    callback = Mock()
    assert run_experiment(experiment, callback) is False
    callback.assert_not_called()


def test_valid_raw_file_and_bruker_directory_are_dispatched_without_viewer_conversion(
    experiment, tmp_path
):
    raw, bruker = tmp_path / "thermo.RAW", tmp_path / "bruker.d"
    raw.touch()
    bruker.mkdir()
    experiment.params["mzML-files"] = [str(raw), str(bruker)]
    experiment.params.pop("_sample_design")
    selected = []
    assert run_experiment(
        experiment, lambda child: selected.append(child.params["mzML-files"]) or True
    )
    assert selected == [str(raw), str(bruker)]


@pytest.mark.parametrize(
    "edit",
    [
        "one-file",
        "missing-design",
        "wrong-count",
        "not-list",
        "bad-row",
        "missing-file",
        "wrong-file",
        "duplicate-row",
        "blank-condition",
        "bad-condition",
        "rep-zero",
        "rep-float",
        "rep-bool",
        "fraction-two",
        "fraction-missing",
        "bad-measure",
        "bad-normalization",
        "model-none",
        "bad-defaults",
        "bad-model-settings",
    ],
)
def test_lfq_preflights_the_complete_design_and_options_before_any_sample(
    experiment, edit
):
    params = experiment.params
    params["enable-lfq"] = True
    if edit == "one-file":
        params["mzML-files"] = params["mzML-files"][:1]
    elif edit == "missing-design":
        params.pop("_sample_design")
    elif edit == "wrong-count":
        params["_sample_design"].pop()
    elif edit == "not-list":
        params["_sample_design"] = {}
    elif edit == "bad-row":
        params["_sample_design"][0] = []
    elif edit == "missing-file":
        params["_sample_design"][0].pop("file")
    elif edit == "wrong-file":
        params["_sample_design"][0]["file"] = "not-selected.mzML"
    elif edit == "duplicate-row":
        params["_sample_design"][1] = deepcopy(params["_sample_design"][0])
    elif edit in {"blank-condition", "bad-condition"}:
        params["_sample_design"][0]["condition"] = (
            " " if edit == "blank-condition" else None
        )
    elif edit.startswith("rep-"):
        params["_sample_design"][0]["replicate"] = {
            "rep-zero": 0,
            "rep-float": 1.5,
            "rep-bool": True,
        }[edit]
    elif edit == "fraction-two":
        params["_sample_design"][0]["fraction"] = 2
    elif edit == "fraction-missing":
        params["_sample_design"][0].pop("fraction")
    elif edit == "bad-measure":
        params["lfq-measure"] = "unknown"
    elif edit == "bad-normalization":
        params["lfq-normalization"] = "unknown"
    else:
        params["lfq-measure"] = "model_area"
        if edit == "bad-defaults":
            params["_defaults"] = []
        elif edit == "bad-model-settings":
            params["FeatureFinderMetaboIdent"] = None
    original = deepcopy(params)
    callback = Mock()
    assert run_experiment(experiment, callback) is False
    callback.assert_not_called()
    assert experiment.params == original
    assert load_manifest(experiment.results_dir) is None


@pytest.mark.parametrize("failure", [False, RuntimeError("native tool failed")])
def test_failed_samples_are_recorded_later_samples_still_run_and_no_comparison_is_claimed(
    experiment, monkeypatch, failure
):
    import src.naseweis_comparison as comparison

    experiment.params["enable-lfq"] = True
    compare = Mock()
    monkeypatch.setattr(comparison, "build_comparison", compare)
    seen = []

    def execute(child):
        seen.append(child.params["mzML-files"])
        (child.results_dir / "partial.txt").write_text("retained partial output")
        if len(seen) == 1:
            if isinstance(failure, Exception):
                raise failure
            return failure
        return True

    assert run_experiment(experiment, execute) is False
    assert seen == experiment.params["mzML-files"]
    compare.assert_not_called()
    manifest = load_manifest(experiment.results_dir)
    assert manifest["status"] == manifest["comparison_status"] == "failed"
    assert [sample["status"] for sample in manifest["samples"]] == [
        "failed",
        "finished",
    ]
    assert manifest["samples"][0]["error"]
    assert (
        sample_results(experiment.results_dir, "sample_001") / "partial.txt"
    ).is_file()
    assert (
        "ERROR"
        in (
            sample_snapshot(experiment.results_dir, "sample_001").parent
            / "logs/all.log"
        ).read_text()
    )


def test_comparison_failure_preserves_finished_searches_without_claiming_a_matrix(
    experiment,
):
    experiment.params["enable-lfq"] = True
    assert (
        run_experiment(experiment, lambda child: True) is False
    )  # Missing native feature tables.
    manifest = load_manifest(experiment.results_dir)
    assert manifest["status"] == manifest["comparison_status"] == "failed"
    assert all(sample["status"] == "finished" for sample in manifest["samples"])
    assert manifest["error"]


@pytest.mark.parametrize("first_succeeds", [True, False])
def test_native_generated_database_is_recycled_across_samples_even_after_partial_failure(
    experiment, first_succeeds
):
    experiment.params["_shared_database"] = "/imported/never-used.fasta"
    experiment._shared_database = "/stale-object/database.fasta"
    experiment._generated_database = "/stale-object/database.fasta"
    children, generated = [], []

    def execute(child):
        children.append(child)
        assert child._experiment_sample_id == f"sample_{len(children):03d}"
        if len(children) == 1:
            assert not hasattr(child, "_shared_database")
            assert not hasattr(child, "_generated_database")
            path = child.results_dir / "database/target-decoy.fasta"
            path.parent.mkdir()
            path.write_text(">target\nACGU\n>DECOY_target\nUGCA\n")
            generated.append(path)
            child._generated_database = str(path)
            return first_succeeds
        assert child._shared_database == str(generated[0])
        child._generated_database = child._shared_database
        return True

    original = deepcopy(experiment.params)
    assert run_experiment(experiment, execute) is first_succeeds
    assert len(children) == 2 and len(generated) == 1
    assert experiment.params == original
    assert experiment._shared_database == "/stale-object/database.fasta"
    assert experiment._generated_database == "/stale-object/database.fasta"


@pytest.mark.parametrize(
    "kind", ["outside", "linked-outside", "missing", "invalid-type"]
)
def test_generated_database_must_exist_inside_the_controlled_experiment_before_reuse(
    experiment, tmp_path, kind
):
    seen = []

    def execute(child):
        assert not hasattr(child, "_shared_database")
        seen.append(child)
        if len(seen) == 1:
            external = tmp_path / "external.fasta"
            external.write_text(">outside\nACGU\n")
            if kind == "linked-outside":
                linked = child.results_dir / "database.fasta"
                linked.symlink_to(external)
                child._generated_database = str(linked)
            elif kind == "missing":
                child._generated_database = str(child.results_dir / "missing.fasta")
            elif kind == "invalid-type":
                child._generated_database = []
            else:
                child._generated_database = str(external)
        return True

    assert run_experiment(experiment, execute) is True
    assert len(seen) == 2


@pytest.mark.parametrize(
    "identity",
    [
        "../outside",
        "sample_001/../../outside",
        "/sample_001",
        "sample_000",
        "SAMPLE_001",
        "sample_1",
        None,
        [],
    ],
)
def test_artifact_paths_require_a_controlled_sample_id(tmp_path, identity):
    for helper in [sample_results, sample_snapshot]:
        with pytest.raises(ValueError, match="sample ID"):
            helper(tmp_path, identity)


def test_controlled_artifacts_cannot_follow_symlinked_directories_or_files(tmp_path):
    root, outside = tmp_path / "results", tmp_path / "outside"
    outside.mkdir()
    (root / "samples").mkdir(parents=True)
    (root / "samples/sample_001").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic links"):
        sample_results(root, "sample_001")
    (root / "samples/sample_001").unlink()
    (root / "samples/sample_001").mkdir()
    (root / "samples/sample_001/params.json").symlink_to(outside / "params.json")
    with pytest.raises(ValueError, match="symbolic links"):
        sample_snapshot(root, "sample_001")


def test_manifest_cannot_be_a_symlink_to_another_file(experiment, tmp_path):
    outside = tmp_path / "other-manifest.json"
    outside.write_text("{}")
    (experiment.results_dir / "experiment.json").symlink_to(outside)
    with pytest.raises(ValueError, match="symbolic link"):
        load_manifest(experiment.results_dir)


def test_manifest_reader_handles_absent_corrupt_and_inconsistent_records(experiment):
    path = experiment.results_dir / "experiment.json"
    assert load_manifest(experiment.results_dir) is None
    path.write_text("broken JSON")
    with pytest.raises(ValueError, match="manifest"):
        load_manifest(experiment.results_dir)
    assert run_experiment(experiment, lambda child: True) is True
    original = json.loads(path.read_text())
    changes = [
        lambda value: [],
        lambda value: value | {"schema": True},
        lambda value: value | {"status": []},
        lambda value: value | {"comparison_status": {}},
        lambda value: value | {"samples": []},
        lambda value: value | {"samples": [None]},
        lambda value: value | {"samples": [value["samples"][0], value["samples"][0]]},
        lambda value: (
            value | {"samples": [value["samples"][0] | {"results": "../../elsewhere"}]}
        ),
        lambda value: (
            value | {"samples": [value["samples"][0] | {"file": "../outside.mzML"}]}
        ),
        lambda value: (
            value
            | {"samples": [value["samples"][0] | {"path": "/different-name.mzML"}]}
        ),
        lambda value: value | {"samples": [value["samples"][0] | {"status": "failed"}]},
        lambda value: value | {"lfq_enabled": True},
    ]
    for edit in changes:
        path.write_text(json.dumps(edit(deepcopy(original))))
        with pytest.raises(ValueError):
            load_manifest(experiment.results_dir)


@pytest.mark.parametrize(
    "content",
    [
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\tMSstats_Condition\tMSstats_BioReplicate\n1\t1\ta.mzML\t1\t1\tcontrol\t2\n",
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n1\t1\ta.mzML\t1\t1\n\nSample\tMSstats_Condition\tMSstats_BioReplicate\n1\tcontrol\t2\n",
        "file\tcondition\treplicate\tfraction\na.mzML\tcontrol\t2.0\t1e0\n",
        "source name\tcomment[data file]\tfactor value[condition]\tcharacteristics[biological replicate]\tcomment[fraction identifier]\tcharacteristics[age]\ncontrol1\ta.mzML\tcontrol\t2\t1\tnot a number\n",
        "source name\tcomment[data file]\tfactor value[condition]\ncontrol1\ta.mzML\tcontrol\n",
        b"\xef\xbb\xbffile\tReplicate\tFraction\r\na.mzML\t1\t1\r\n",
    ],
)
def test_strict_design_source_accepts_openms_one_two_table_and_sdrf_formats(content):
    assert validate_design_source(content) is None


@pytest.mark.parametrize(
    "column",
    [
        "Fraction",
        "MSstats_BioReplicate",
        "replicate",
        "fraction",
        "biological replicate",
        "biological_replicate",
        "characteristics[biological replicate]",
        "comment[biological replicate]",
        "comment[fraction identifier]",
    ],
)
@pytest.mark.parametrize(
    "value",
    [
        "0",
        "-1",
        "1.5",
        "nan",
        "inf",
        "",
        "True",
        "False",
        "1.00000000000000000000000000001",
    ],
)
def test_raw_scientific_numeric_fields_are_rejected_before_normalization(column, value):
    with pytest.raises(ValueError, match="positive integer"):
        validate_design_source(f"file\t{column}\na.mzML\t{value}\n")


def test_strict_design_validation_checks_both_tables_and_present_missing_cells():
    content = (
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n"
        "1\t1\ta.mzML\t1\t1\n\n"
        "Sample\tMSstats_Condition\tMSstats_BioReplicate\n"
        "1\tcontrol\t1.5\n"
    )
    with pytest.raises(ValueError, match="MSstats_BioReplicate on line 5"):
        validate_design_source(content)
    with pytest.raises(ValueError, match="replicate on line 2"):
        validate_design_source("file\treplicate\na.mzML\n")
    with pytest.raises(ValueError, match="positive integer"):
        validate_design_source("replicate\tfraction\n\t\n")


@pytest.mark.parametrize(
    "content", [None, 123, b"\xff", 'file\treplicate\na.mzML\t"unterminated']
)
def test_invalid_design_text_has_an_actionable_error(content):
    with pytest.raises(ValueError, match="experimental design"):
        validate_design_source(content)
