"""Pinned NME1 integration from NASE/LFQ outputs to the existing viewers."""

import json
import re

import numpy as np
import pandas as pd
import pytest
from test_naseweis_lfq_native import lfq_runtime  # noqa: F401
from test_naseweis_nightly import ROOT, nightly_runtime, run_search  # noqa: F401


def test_real_nme1_run_preserves_native_scan_identity_signal_and_rna_annotations(
    tmp_path,
    lfq_runtime,  # noqa: F811
    monkeypatch,
):
    import pyopenms as oms

    from src.naseweis_annotations import (
        annotate_peaks,
        fragment_assignment_table,
        match_annotations,
        native_spectrum_id,
        read_annotations,
    )
    from src.naseweis_lfq import read_lfq_table
    from src.naseweis_msdata import ms1_overview, scan_index, spectrum_peaks, xic
    from src.naseweis_plots import (
        chromatogram_figure,
        feature_map_figure,
        peakmap_figure,
        spectrum_figure,
    )
    from src.naseweis_results import read_mztab

    directory, snapshot = run_search(tmp_path, "nme1-viewers", {"enable-lfq": True})
    acquisition = ROOT / "example-data/mzML/NME1.mzML"
    assert snapshot["mzML-files"] == str(acquisition)
    assert snapshot["enable-lfq"] is True
    manifest = json.loads(
        (ROOT / "tests/data/naseweis/nightly/manifest.json").read_text()
    )
    assert (
        snapshot["_acquisition_fingerprint"]["sha256"]
        == manifest["inputs"]["example-data/mzML/NME1.mzML"]
    )
    assert (directory / "idxml_results/search.idXML").is_file()

    # The independent reference deliberately uses the native experiment API,
    # rather than the viewer's scan index or overview to compute expectations.
    original = oms.MSExperiment()
    oms.MzMLFile().load(str(acquisition), original)
    expected_ids = [spectrum.getNativeID() for spectrum in original]
    ms1 = [spectrum for spectrum in original if spectrum.getMSLevel() == 1]
    assert len(ms1) == 512

    def no_full_dataframe(*_args, **_kwargs):
        raise AssertionError(
            "Viewer readers must not materialize a full experiment DataFrame."
        )

    monkeypatch.setattr(oms.MSExperiment, "get_df", no_full_dataframe)
    scans = scan_index(acquisition)
    assert scans["native_id"].tolist() == expected_ids
    assert scans["index"].tolist() == list(range(len(original)))
    by_id = scans.set_index("native_id", verify_integrity=True)
    matches = read_mztab(directory / "mztab_results/search.mzTab")["OSM"]
    native_ids = matches["spectra_ref"].map(native_spectrum_id)
    assert not matches.empty
    assert set(native_ids) <= set(by_id.index)
    for native_id in native_ids.unique():
        index = int(by_id.loc[native_id, "index"])
        assert original[index].getNativeID() == native_id
        assert by_id.loc[native_id, "ms_level"] == 2
    # This fixture distinguishes actual scan numbers from zero-based indices.
    assert any(
        int(re.search(r"scan=(\d+)$", native_id).group(1))
        not in {
            int(by_id.loc[native_id, "index"]),
            int(by_id.loc[native_id, "index"]) + 1,
        }
        for native_id in native_ids.unique()
    )

    annotations = read_annotations(directory / "idxml_results/search.idXML")
    rna_specific = annotations["ion"].str.contains(
        r"^(?:a\d+-B|a-B|[dw]\d)", regex=True
    )
    assert rna_specific.any()
    candidate = annotations.loc[rna_specific].iloc[0]
    selected_matches = matches.loc[
        (native_ids == candidate["native_id"])
        & (matches["sequence"] == candidate["sequence"])
        & (matches["charge"] == candidate["charge"])
    ]
    assert not selected_matches.empty
    match = selected_matches.iloc[0]
    selected = match_annotations(
        annotations, match["spectra_ref"], match["sequence"], int(match["charge"])
    )
    assert not selected.empty
    assert (selected["charge"] < 0).all()
    assert selected["native_id"].eq(candidate["native_id"]).all()
    assert selected["sequence"].eq(match["sequence"]).all()
    assert match_annotations(
        annotations, match["spectra_ref"], match["sequence"], -int(match["charge"])
    ).empty

    index = int(by_id.loc[candidate["native_id"], "index"])
    peaks = spectrum_peaks(acquisition, index)
    expected_mz, expected_intensity = original[index].get_peaks()
    np.testing.assert_array_equal(peaks["mz"], expected_mz)
    np.testing.assert_array_equal(peaks["intensity"], expected_intensity)
    labelled, unmatched = annotate_peaks(peaks, selected)
    pd.testing.assert_frame_equal(labelled[["mz", "intensity"]], peaks)
    assert labelled["ion"].ne("").any()
    for assignment in selected.itertuples(index=False):
        positions = np.flatnonzero(np.abs(expected_mz - assignment.mz) <= 1e-5)
        if len(positions) == 1:
            assert assignment.ion in labelled.iloc[positions[0]]["ion"].split("; ")
        else:
            assert assignment.ion in unmatched["ion"].values
    spectrum = spectrum_figure(labelled, match["sequence"])
    assert spectrum.layout.annotations
    assert any(
        re.match(r"^(?:a\d+-B|a-B|[dw]\d)", annotation.text)
        for annotation in spectrum.layout.annotations
    )
    evidence = fragment_assignment_table(peaks, selected)
    placed = evidence.loc[evidence["match_status"].eq("matched")].iloc[0]
    raw_position = int(placed["peak_index"])
    assert placed["observed_mz"] == expected_mz[raw_position]
    assert placed["observed_intensity"] == expected_intensity[raw_position]
    highlighted = spectrum_figure(
        labelled, match["sequence"], highlighted_mz=[placed["observed_mz"]]
    )
    selected_trace = next(trace for trace in highlighted.data if trace.name == "Selected ion")
    assert list(selected_trace.x) == [expected_mz[raw_position]]
    np.testing.assert_allclose(
        selected_trace.y,
        [float(expected_intensity[raw_position]) / float(expected_intensity.max()) * 100],
    )
    pd.testing.assert_frame_equal(labelled[["mz", "intensity"]], peaks)

    overview, traces = ms1_overview(acquisition)
    expected_rt = [spectrum.getRT() for spectrum in ms1]
    expected_tic = [
        float(np.sum(spectrum.get_peaks()[1], dtype=np.float64)) for spectrum in ms1
    ]
    expected_bpc = [float(np.max(spectrum.get_peaks()[1])) for spectrum in ms1]
    for name, values in [("TIC", expected_tic), ("BPC", expected_bpc)]:
        actual = traces.loc[traces["trace"] == name]
        assert len(actual) == 512
        np.testing.assert_array_equal(actual["rt"], expected_rt)
        np.testing.assert_allclose(actual["intensity"], values, rtol=0, atol=1e-8)
    assert sum(spectrum.size() for spectrum in ms1) > 40000
    assert 0 < len(overview) <= 40000
    assert overview["intensity"].sum() == pytest.approx(sum(expected_tic), rel=1e-12)
    peakmap = peakmap_figure(overview, "NME1 MS1 signal")
    assert len(peakmap.data[0].x) == len(overview)
    chromatograms = chromatogram_figure(traces, "NME1 TIC/BPC")
    assert {trace.name for trace in chromatograms.data} == {"TIC", "BPC"}

    features = read_lfq_table(directory / "quantification/features.tsv")
    assert not features.empty
    feature = features.loc[features["intensity"] > 0].iloc[0]
    tolerance = snapshot["FeatureFinderMetaboIdent"]["extract:mz_window"] / 2
    width = float(feature["mz"]) * tolerance / 1e6
    trace = xic(acquisition, float(feature["mz"]), tolerance, "ppm")
    expected_xic = []
    for spectrum in ms1:
        masses, intensities = spectrum.get_peaks()
        mask = (masses >= feature["mz"] - width) & (masses <= feature["mz"] + width)
        expected_xic.append(float(np.sum(intensities[mask], dtype=np.float64)))
    assert len(trace) == 512
    np.testing.assert_array_equal(trace["rt"], expected_rt)
    np.testing.assert_allclose(trace["intensity"], expected_xic, rtol=0, atol=1e-8)
    assert (trace["intensity"] >= 0).all() and (trace["intensity"] > 0).any()
    assert len(chromatogram_figure(trace, str(feature["label"])).data) == 1
    feature_map = feature_map_figure(features, "NME1 RNA features")
    hover = feature_map.data[0].customdata
    assert list(hover[:, 0]) == features["label"].tolist()
    np.testing.assert_array_equal(hover[:, 1], features["charge"])
    np.testing.assert_array_equal(hover[:, 2], features["intensity"])
    assert (features["charge"] < 0).all()
