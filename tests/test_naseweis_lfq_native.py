"""Pinned-tool LFQ integration, including real RNA data and a known peak area.

The small MS1 fixture is generated here from an RNA empirical formula and a
Gaussian elution profile. It is original synthetic test data under the
repository's BSD-3-Clause license; no experimental data are copied into it.
"""

import json
import math
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from test_naseweis_nightly import ROOT, nightly_runtime, run_search  # noqa: F401


@pytest.fixture
def lfq_runtime(nightly_runtime):  # noqa: F811
    missing = [
        tool
        for tool in ["FeatureFinderMetaboIdent", "TextExporter"]
        if not shutil.which(tool)
    ]
    if missing:
        pytest.skip(
            "The pinned nightly LFQ tools are not installed: " + ", ".join(missing)
        )
    revision = json.loads((ROOT / "openms-runtime.json").read_text())["revision"]
    for tool in ["FeatureFinderMetaboIdent", "TextExporter"]:
        result = subprocess.run(
            [tool, "--help"], capture_output=True, text=True, timeout=30, check=False
        )
        if f"Revision: {revision}" not in result.stdout + result.stderr:
            pytest.skip(f"{tool} does not match the pinned OpenMS nightly runtime.")


@pytest.fixture
def rna_ms1(tmp_path, lfq_runtime):
    """One negative-mode RNA isotope envelope with a known, unscaled area."""
    import pyopenms as oms

    sequence = oms.NASequence.fromString("ACGUAG")
    formula = sequence.getFormula()
    # Neutral formula, two removed protons, and C13 spacing for charge -2.
    monoisotopic_mz = (formula.getMonoWeight() - 2 * 1.007276466621) / 2
    weights = [
        peak.getIntensity()
        for peak in formula.getIsotopeDistribution(
            oms.CoarseIsotopePatternGenerator(4)
        ).getContainer()
    ]
    experiment = oms.MSExperiment()
    for rt in range(121):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(1)
        spectrum.setRT(float(rt))
        spectrum.setNativeID(f"scan={rt + 1}")
        spectrum.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
        settings = oms.InstrumentSettings()
        settings.setPolarity(oms.IonSource.Polarity.NEGATIVE)
        spectrum.setInstrumentSettings(settings)
        height = 200_000 * math.exp(-0.5 * ((rt - 60) / 8) ** 2)
        # A fixed low baseline and unrelated peaks make peak finding exercise
        # actual extracted chromatograms rather than a prebuilt feature map.
        peaks = sorted(
            [
                (monoisotopic_mz + index * 1.0033548378 / 2, height * weight + 1)
                for index, weight in enumerate(weights)
            ]
            + [(900 + index * 2.3, 10.0) for index in range(10)]
        )
        spectrum.set_peaks(
            (
                np.array([peak[0] for peak in peaks]),
                np.array([peak[1] for peak in peaks]),
            )
        )
        experiment.addSpectrum(spectrum)
    spectra = tmp_path / "synthetic-rna.mzML"
    oms.MzMLFile().store(str(spectra), experiment)
    assert spectra.stat().st_size < 300_000
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "CompoundName\tSumFormula\tMass\tCharge\tRetentionTime\t"
        "RetentionTimeRange\tIsoDistribution\n"
        f"ACGUAG\t{formula.toString()}\t0\t-2\t60\t0\t0\n"
    )
    return spectra, targets, {
        "formula": formula.toString(),
        "mz": monoisotopic_mz,
        "area": 200_000 * math.sqrt(2 * math.pi) * 8,
    }


def lfq_subject(tmp_path, spectra):
    """Use the real parameter manager, dispatcher and LFQ helper without NASE."""
    from src.naseweis_lfq import prepare_lfq
    from src.workflow.CommandExecutor import CommandExecutor
    from src.workflow.FileManager import FileManager
    from src.workflow.Logger import Logger
    from src.workflow.ParameterManager import ParameterManager

    directory = tmp_path / "naseweis"
    directory.mkdir()
    (directory / "pids").mkdir()
    results = directory / "results"
    results.mkdir()
    manager = ParameterManager(directory)
    params = {"enable-lfq": True, "mzML-files": str(spectra), "max_threads": 1}
    assert prepare_lfq(params, manager, str(spectra)) is None
    manager.write_parameters(params)
    logger = Logger(directory)
    return SimpleNamespace(
        workflow_dir=directory,
        results_dir=results,
        params=params,
        logger=logger,
        parameter_manager=manager,
        file_manager=FileManager(directory, results),
        executor=CommandExecutor(directory, logger, manager),
    )


def test_native_rna_extraction_and_export_preserve_identity_and_peak_area(
    tmp_path, rna_ms1
):
    import pyopenms as oms

    from src.naseweis_lfq import run_lfq

    spectra, targets, expected = rna_ms1
    subject = lfq_subject(tmp_path, spectra)
    assert run_lfq(subject, targets) is True
    features = oms.FeatureMap()
    oms.FeatureXMLFile().load(
        str(subject.results_dir / "quantification/features.featureXML"), features
    )
    assert features.size() == 1
    feature = features[0]
    assert feature.getCharge() == -2
    assert feature.getRT() == pytest.approx(60, abs=1)
    assert feature.getMZ() == pytest.approx(expected["mz"], rel=2e-6)
    assert feature.getIntensity() == pytest.approx(expected["area"], rel=0.05)
    table = pd.read_csv(
        subject.results_dir / "quantification/features.tsv", sep="\t"
    )
    assert len(table) == 1
    assert table.loc[0, "label"] == "ACGUAG"
    assert table.loc[0, "sum_formula"] == expected["formula"]
    assert table.loc[0, "charge"] == -2
    assert table.loc[0, "expected_rt"] == 60
    assert table.loc[0, "intensity"] == pytest.approx(feature.getIntensity(), rel=1e-5)
    commands = (subject.workflow_dir / "logs/commands-and-run-times.log").read_text()
    assert "FeatureFinderMetaboIdent" in commands
    assert "TextExporter" in commands


def test_native_empty_target_list_is_a_successful_empty_result(tmp_path, rna_ms1):
    import pyopenms as oms

    from src.naseweis_lfq import run_lfq

    spectra, targets, _ = rna_ms1
    targets.write_text(targets.read_text().splitlines()[0] + "\n")
    subject = lfq_subject(tmp_path, spectra)
    assert run_lfq(subject, targets) is True
    features = oms.FeatureMap()
    oms.FeatureXMLFile().load(
        str(subject.results_dir / "quantification/features.featureXML"), features
    )
    assert features.size() == 0
    assert pd.read_csv(
        subject.results_dir / "quantification/features.tsv", sep="\t"
    ).empty


def test_native_workflow_quantifies_the_original_nme1_ms1_data(tmp_path, lfq_runtime):
    """Exercise NASE target generation and the full worker execution path."""
    import pyopenms as oms

    directory, snapshot = run_search(
        tmp_path,
        "nme1-lfq",
        {"enable-lfq": True, "search-exports": []},
    )
    assert snapshot["enable-lfq"] is True
    assert snapshot["FeatureFinderMetaboIdent"]["extract:mz_window"] == 50.0
    assert snapshot["FeatureFinderMetaboIdent"]["model:type"] == "none"
    assert snapshot["TextExporter"]["feature:add_metavalues"] == 0
    assert not pd.read_csv(directory / "exports/lfq-targets.tsv", sep="\t").empty
    features = oms.FeatureMap()
    oms.FeatureXMLFile().load(
        str(directory / "quantification/features.featureXML"), features
    )
    table = pd.read_csv(directory / "quantification/features.tsv", sep="\t")
    assert features.size() == len(table) > 0
    assert {"label", "sum_formula", "charge", "expected_rt", "intensity"} <= set(table)
    assert np.isfinite(table["intensity"]).all()
    assert (table["intensity"] > 0).any()
    assert (table["charge"] < 0).all()
