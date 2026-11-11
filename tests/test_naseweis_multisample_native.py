"""Pinned native cohort checks with known twofold MS1 abundance differences."""

import hashlib
import json
import shutil
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_naseweis_lfq_native import (  # noqa: F401
    lfq_runtime,
    lfq_subject,
    rna_ms1,
)
from test_naseweis_nightly import ROOT, nightly_runtime  # noqa: F401


def checksum(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def double_ms1(source, destination):
    """Scale MS1 signal only; MS2 peak arrays and native IDs stay unchanged."""
    import pyopenms as oms

    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(source), experiment)
    spectra = experiment.getSpectra()
    for spectrum in spectra:
        if spectrum.getMSLevel() == 1:
            masses, intensities = spectrum.get_peaks()
            spectrum.set_peaks((masses, intensities * 2))
    experiment.setSpectra(spectra)
    oms.MzMLFile().store(str(destination), experiment)

    scaled = oms.MSExperiment()
    oms.MzMLFile().load(str(destination), scaled)
    assert scaled.size() == experiment.size()
    for before, after in zip(experiment, scaled):
        assert before.getNativeID() == after.getNativeID()
        assert before.getMSLevel() == after.getMSLevel()
        assert before.getRT() == after.getRT()
        np.testing.assert_array_equal(before.get_peaks()[0], after.get_peaks()[0])
        np.testing.assert_array_equal(before.get_peaks()[1], after.get_peaks()[1])


def cohort_design(files):
    return [
        {"file": str(path), "condition": condition, "replicate": 1, "fraction": 1}
        for path, condition in zip(files, ["control", "treatment"])
    ]


def test_native_nme1_cohort_shares_database_and_recovers_twofold_ms1_signal(
    tmp_path,
    lfq_runtime,  # noqa: F811
    monkeypatch,
):
    import pyopenms as oms

    from src.naseweis_comparison import compare_conditions, load_comparison
    from src.naseweis_experiment import load_manifest, sample_results, sample_snapshot
    from src.naseweis_lfq import read_lfq_table
    from src.naseweis_results import read_mztab
    from src.workflow.CommandExecutor import CommandExecutor
    from src.workflow.RunManager import RunManager
    from src.workflow.tasks import execute_workflow

    source = ROOT / "example-data/mzML/NME1.mzML"
    source_hash = checksum(source)
    control, treatment = tmp_path / "control.mzML", tmp_path / "treatment.mzML"
    shutil.copyfile(source, control)
    double_ms1(source, treatment)
    # Verify the intentional difference independently of the scaling helper.
    original, scaled = oms.MSExperiment(), oms.MSExperiment()
    oms.MzMLFile().load(str(control), original)
    oms.MzMLFile().load(str(treatment), scaled)
    for before, after in zip(original, scaled):
        expected = before.get_peaks()[1] * (2 if before.getMSLevel() == 1 else 1)
        np.testing.assert_array_equal(after.get_peaks()[1], expected)

    files = [control, treatment]
    params = {
        "mzML-files": [str(path) for path in files],
        "_sample_design": cohort_design(files),
        "fasta-files": str(ROOT / "example-data/fasta/NME1_Target.fasta"),
        "add-decoys": True,
        "estimate-fdr": True,
        "FDR_cutoff": 0.05,
        "enable-lfq": True,
        "lfq-measure": "intensity",
        "lfq-normalization": "total",
        "max_threads": 2,
    }
    workflow = tmp_path / "workspace/naseweis"
    workflow.mkdir(parents=True)
    configuration = workflow / "params.json"
    configuration.write_text(json.dumps(params))
    runs = RunManager(workflow)
    run_id = runs.new_run("NME1 twofold cohort", configuration)
    saved = runs.params_file(run_id).read_bytes()
    commands = []
    execute_command = CommandExecutor.run_command

    def record_command(executor, command, *args, **kwargs):
        commands.append(command[:])
        return execute_command(executor, command, *args, **kwargs)

    monkeypatch.setattr(CommandExecutor, "run_command", record_command)
    try:
        execute_workflow(str(workflow), "Workflow", "src.Workflow", run_id)
    except Exception:
        print((workflow / "logs/all.log").read_text())
        raise
    results = runs.results_dir(run_id)
    manifest = load_manifest(results)
    assert runs.get(run_id)["status"] == "finished"
    assert manifest["status"] == manifest["comparison_status"] == "finished"
    assert runs.params_file(run_id).read_bytes() == saved
    assert json.loads(configuration.read_text()) == params
    assert checksum(source) == source_hash

    tools = [Path(command[0]).name for command in commands]
    assert tools.count("DecoyDatabase") == 1
    for tool in ["NucleicAcidSearchEngine", "FeatureFinderMetaboIdent", "TextExporter"]:
        assert tools.count(tool) == 2
    databases, identifications = [], []
    for sample, acquisition in zip(manifest["samples"], files):
        assert sample["status"] == "finished"
        snapshot = json.loads(sample_snapshot(results, sample["id"]).read_text())
        assert snapshot["mzML-files"] == str(acquisition)
        assert snapshot["_sample_design"] == [params["_sample_design"][len(databases)]]
        assert snapshot["_acquisition_fingerprint"]["sha256"] == checksum(acquisition)
        assert snapshot["_acquisition_fingerprint"]["path"] == str(
            acquisition.resolve()
        )
        databases.append(snapshot["_database_used"])
        directory = sample_results(results, sample["id"])
        assert (directory / "idxml_results/search.idXML").is_file()
        assert (directory / "quantification/features.featureXML").is_file()
        table = read_lfq_table(directory / "quantification/features.tsv")
        assert not table.empty and (table["charge"] < 0).all()
        identifications.append(
            read_mztab(directory / "mztab_results/search.mzTab")["OSM"]
        )
    assert databases[0] == databases[1]
    assert Path(databases[0]).is_relative_to(results)
    native_searches = [
        command
        for command in commands
        if Path(command[0]).name == "NucleicAcidSearchEngine"
    ]
    assert [
        command[command.index("-database") + 1] for command in native_searches
    ] == databases
    columns = ["sequence", "spectra_ref", "charge", "search_engine_score[1]"]
    pd.testing.assert_frame_equal(
        identifications[0][columns].sort_values(columns).reset_index(drop=True),
        identifications[1][columns].sort_values(columns).reset_index(drop=True),
    )

    comparison = load_comparison(results)
    raw, relative = comparison["raw"], comparison["abundances"]
    common = raw.dropna()
    assert not common.empty
    abundant = common.sum(axis=1).idxmax()
    # Peak picking can move a boundary slightly when signal is scaled. The
    # strongest common RNA avoids threshold effects; tolerate 0.5% in raw area
    # and 1% after total normalization across the complete extracted cohort.
    assert raw.loc[abundant, "sample_002"] / raw.loc[
        abundant, "sample_001"
    ] == pytest.approx(2, rel=0.005)
    assert relative.loc[abundant, "sample_002"] / relative.loc[
        abundant, "sample_001"
    ] == pytest.approx(1, rel=0.01)
    summary = compare_conditions(relative, comparison["design"], "treatment", "control")
    assert summary.loc[abundant, "ratio"] == pytest.approx(1, rel=0.01)
    assert (
        summary.loc[abundant, "numerator_n"]
        == summary.loc[abundant, "denominator_n"]
        == 1
    )


def test_native_fitted_rna_area_cohort_recovers_known_gaussian_ratio(
    tmp_path,
    rna_ms1,  # noqa: F811
):
    from src.naseweis_comparison import compare_conditions, load_comparison
    from src.naseweis_experiment import load_manifest, run_experiment, sample_results
    from src.naseweis_lfq import prepare_lfq, read_lfq_table, run_lfq

    first, targets, expected = rna_ms1
    second = tmp_path / "synthetic-rna-twofold.mzML"
    double_ms1(first, second)
    parent = tmp_path / "execution"
    parent.mkdir()
    subject = lfq_subject(parent, first)
    subject.params.update(
        {
            "mzML-files": [str(first), str(second)],
            "_sample_design": cohort_design([first, second]),
            "lfq-measure": "model_area",
            "lfq-normalization": "total",
        }
    )
    subject.params["FeatureFinderMetaboIdent"]["model:type"] = "symmetric"
    subject.parameter_manager.write_parameters(subject.params)
    before, saved = (
        deepcopy(subject.params),
        subject.parameter_manager.params_file.read_bytes(),
    )

    def extract_native(child):
        error = prepare_lfq(
            child.params, child.parameter_manager, child.params["mzML-files"]
        )
        assert error is None
        child.parameter_manager.write_parameters(child.params)
        return run_lfq(child, targets)

    assert run_experiment(subject, extract_native) is True
    assert subject.params == before
    assert subject.parameter_manager.params_file.read_bytes() == saved
    manifest = load_manifest(subject.results_dir)
    assert manifest["status"] == manifest["comparison_status"] == "finished"
    for sample in manifest["samples"]:
        table = read_lfq_table(
            sample_results(subject.results_dir, sample["id"])
            / "quantification/features.tsv"
        )
        assert len(table) == 1
        assert table.loc[0, "label"] == "ACGUAG"
        assert table.loc[0, "charge"] == -2
        assert table.loc[0, "sum_formula"] == expected["formula"]
        assert float(table.loc[0, "model_area"]) > 0
    comparison = load_comparison(subject.results_dir)
    raw, relative = comparison["raw"], comparison["abundances"]
    assert comparison["metadata"]["measure"] == "model_area"
    assert raw.loc["ACGUAG", "sample_001"] == pytest.approx(expected["area"], rel=0.05)
    # The two independent nonlinear fits may have slightly different numerical
    # convergence; a 1% ratio tolerance is much tighter than the area gate.
    assert raw.loc["ACGUAG", "sample_002"] / raw.loc[
        "ACGUAG", "sample_001"
    ] == pytest.approx(2, rel=0.01)
    assert relative.loc["ACGUAG"].tolist() == pytest.approx([1, 1])
    summary = compare_conditions(raw, comparison["design"], "treatment", "control")
    assert summary.loc["ACGUAG", "ratio"] == pytest.approx(2, rel=0.01)
    assert summary.loc["ACGUAG", "log2_fold_change"] == pytest.approx(1, abs=0.015)
