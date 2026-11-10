"""Real native mzML coverage for bounded acquisition views and exact traces."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms
import pytest

from src import naseweis_msdata as msdata


def _spectrum(level, rt, masses=(), intensities=(), native_id="scan=1", precursor=None):
    spectrum = pyopenms.MSSpectrum()
    spectrum.setMSLevel(level)
    spectrum.setRT(rt)
    spectrum.setNativeID(native_id)
    spectrum.set_peaks(
        (np.array(masses, dtype=float), np.array(intensities, dtype=np.float32))
    )
    if precursor is not None:
        selected = pyopenms.Precursor()
        selected.setMZ(precursor[0])
        selected.setCharge(precursor[1])
        spectrum.setPrecursors([selected])
    return spectrum


def _write(path, spectra, indexed):
    experiment = pyopenms.MSExperiment()
    experiment.setSpectra(spectra)
    writer = pyopenms.MzMLFile()
    options = writer.getOptions()
    options.setWriteIndex(indexed)
    writer.setOptions(options)
    writer.store(str(path), experiment)
    return path


@pytest.fixture(params=[True, False], ids=["indexed", "streamed"])
def acquisition(tmp_path, request):
    return _write(
        tmp_path / "acquisition.mzML",
        [
            _spectrum(1, 10, [99.9, 100, 100.0005, 100.01], [1, 2, 3, 0]),
            _spectrum(2, 10.5, [200, 300], [9, 8], precursor=(500.2, -2)),
            _spectrum(1, 11, [101, 102], [0, 4], native_id="scan=3"),
            _spectrum(1, 12, native_id="scan=4"),
            _spectrum(2, 12.5, [500], [100000], native_id="scan=5"),
        ],
        request.param,
    )


def test_revision_is_canonical_and_changes_with_source(tmp_path):
    source = tmp_path / "data.mzML"
    source.write_text("first")
    revision = msdata.file_revision(source.parent / "." / source.name)
    info = source.stat()
    assert revision == (
        str(source.resolve()),
        info.st_mtime_ns,
        5,
        info.st_ctime_ns,
        info.st_ino,
        info.st_dev,
    )
    source.write_text("second version")
    assert msdata.file_revision(source) != revision


@pytest.mark.parametrize("replacement", ["overwrite", "replace"])
def test_revision_detects_same_size_replacement_with_restored_mtime(
    tmp_path, replacement
):
    source = tmp_path / "data.mzML"
    source.write_bytes(b"original")
    before = msdata.file_revision(source)
    if replacement == "overwrite":
        source.write_bytes(b"replaced")
    else:
        other = tmp_path / "other.mzML"
        other.write_bytes(b"replaced")
        other.replace(source)
    os.utime(source, ns=(before[1], before[1]))
    after = msdata.file_revision(source)
    assert after[:3] == before[:3]
    # Windows may expose creation-time ctime for in-place writes; callers must
    # verify content there rather than cache a fingerprint by stat revision.
    if os.name != "nt" or replacement == "replace":
        assert after != before


@pytest.mark.parametrize("suffix", [".raw", ".d", ".mzXML", ".txt"])
def test_vendor_formats_have_actionable_conversion_message(tmp_path, suffix):
    with pytest.raises(ValueError, match="Convert.*mzML.*NASE searches"):
        msdata.file_revision(tmp_path / f"acquisition{suffix}")


def test_revision_rejects_missing_file_directory_and_invalid_path(tmp_path):
    with pytest.raises(OSError, match="Cannot access mzML"):
        msdata.file_revision(tmp_path / "missing.mzML")
    directory = tmp_path / "directory.mzML"
    directory.mkdir()
    with pytest.raises(OSError, match="not a file"):
        msdata.file_revision(directory)
    with pytest.raises(ValueError, match="existing mzML"):
        msdata.file_revision(None)


def test_scan_index_preserves_acquisition_indices_duplicate_ids_and_precursors(
    acquisition,
):
    scans = msdata.scan_index(acquisition)
    assert list(scans.columns) == msdata.SCAN_COLUMNS
    assert scans["index"].tolist() == list(range(5))
    assert scans["native_id"].tolist() == [
        "scan=1",
        "scan=1",
        "scan=3",
        "scan=4",
        "scan=5",
    ]
    assert scans["ms_level"].tolist() == [1, 2, 1, 1, 2]
    assert scans["rt"].tolist() == [10, 10.5, 11, 12, 12.5]
    assert scans.loc[1, "precursor_mz"] == pytest.approx(500.2)
    assert scans.loc[1, "charge"] == -2
    assert np.isnan(scans.loc[0, "precursor_mz"])
    assert scans.loc[0, "charge"] == 0


def test_spectrum_is_full_resolution_and_selected_by_index_not_native_id(acquisition):
    first = msdata.spectrum_peaks(acquisition, 0)
    np.testing.assert_array_equal(first["mz"], [99.9, 100, 100.0005, 100.01])
    np.testing.assert_array_equal(first["intensity"], [1, 2, 3, 0])
    second = msdata.spectrum_peaks(acquisition, np.int64(1))
    np.testing.assert_array_equal(second["mz"], [200, 300])
    np.testing.assert_array_equal(second["intensity"], [9, 8])
    assert list(msdata.spectrum_peaks(acquisition, 3).columns) == ["mz", "intensity"]
    assert msdata.spectrum_peaks(acquisition, 3).empty


@pytest.mark.parametrize("index", [-1, 99, True, 1.5, "0"])
def test_invalid_scan_indices_are_actionable(acquisition, index):
    with pytest.raises(ValueError, match="scan index|Scan index"):
        msdata.spectrum_peaks(acquisition, index)


@pytest.mark.parametrize("limit", [1, 2, 4, 40000])
def test_overview_is_bounded_sums_peaks_and_keeps_exact_traces(acquisition, limit):
    points, traces = msdata.ms1_overview(acquisition, limit)
    assert list(points.columns) == ["rt", "mz", "intensity"]
    assert 0 < len(points) <= limit
    assert points["intensity"].sum() == 10
    assert points["rt"].between(10, 12).all()
    assert points["mz"].between(99.9, 102).all()
    assert list(traces.columns) == ["rt", "intensity", "trace"]
    assert traces[traces["trace"] == "TIC"]["intensity"].tolist() == [6, 4, 0]
    assert traces[traces["trace"] == "BPC"]["intensity"].tolist() == [3, 4, 0]
    assert traces[traces["trace"] == "TIC"]["rt"].tolist() == [10, 11, 12]


@pytest.mark.parametrize("limit", [0, -1, True, 1.5, "10"])
def test_invalid_overview_limits(acquisition, limit):
    with pytest.raises(ValueError, match="positive integer"):
        msdata.ms1_overview(acquisition, limit)


def test_xic_uses_full_resolution_peaks_and_retains_zero_scans(acquisition):
    trace = msdata.xic(acquisition, 100, 0.001, "Da")
    pd.testing.assert_frame_equal(
        trace,
        pd.DataFrame(
            {
                "rt": [10.0, 11.0, 12.0],
                "intensity": [5.0, 0.0, 0.0],
                "trace": ["XIC"] * 3,
            }
        ),
    )
    assert msdata.xic(acquisition, 100, 0, "Da")["intensity"].tolist() == [2, 0, 0]
    assert msdata.xic(acquisition, 100, 10, "ppm", (10, 11))["rt"].tolist() == [10, 11]
    assert msdata.xic(acquisition, 100, 10, "ppm", (11.5, 12))[
        "intensity"
    ].tolist() == [0]
    assert msdata.xic(acquisition, 100, 10, "ppm", (20, 30)).empty


@pytest.mark.parametrize("indexed", [True, False])
@pytest.mark.parametrize(
    "unit,tolerance,width", [("ppm", 10, 0.005), ("Da", 0.005, 0.005)]
)
def test_xic_includes_both_mass_boundaries_and_duplicate_peaks(
    tmp_path, indexed, unit, tolerance, width
):
    path = _write(
        tmp_path / "boundaries.mzML",
        [
            _spectrum(
                1,
                1,
                [
                    500 - width - 1e-6,
                    500 - width,
                    500,
                    500,
                    500 + width,
                    500 + width + 1e-6,
                ],
                [100, 1, 2, 3, 4, 100],
            )
        ],
        indexed,
    )
    assert msdata.xic(path, 500, tolerance, unit)["intensity"].tolist() == [10]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"mz": 0},
        {"mz": -1},
        {"mz": float("nan")},
        {"mz": True},
        {"mz": "100"},
        {"mz": 10**1000},
        {"tolerance": -1},
        {"tolerance": float("inf")},
        {"tolerance": True},
        {"unit": "mDa"},
        {"unit": None},
        {"rt_range": (2, 1)},
        {"rt_range": (1,)},
        {"rt_range": (1, 2, 3)},
        {"rt_range": (None, 1)},
        {"rt_range": (0, float("inf"))},
        {"rt_range": "12"},
        {"rt_range": (True, 1)},
        {"rt_range": (0, 10**1000)},
        {"mz": 1e308, "tolerance": 1e308},
    ],
)
def test_invalid_xic_settings_fail_before_reading_a_file(tmp_path, kwargs):
    parameters = {"mz": 100, "tolerance": 10} | kwargs
    with pytest.raises(ValueError, match="XIC"):
        msdata.xic(tmp_path / "does-not-exist.mzML", **parameters)


@pytest.mark.parametrize("indexed", [True, False])
def test_no_ms1_has_actionable_error_but_ms2_remains_viewable(tmp_path, indexed):
    path = _write(tmp_path / "ms2.mzML", [_spectrum(2, 1, [100], [5])], indexed)
    assert msdata.scan_index(path)["ms_level"].tolist() == [2]
    assert msdata.spectrum_peaks(path, 0)["intensity"].tolist() == [5]
    for operation in [msdata.ms1_overview, lambda source: msdata.xic(source, 100, 1)]:
        with pytest.raises(ValueError, match="no MS1.*Choose an mzML"):
            operation(path)


@pytest.mark.parametrize("indexed", [True, False])
def test_large_overview_retains_weak_peaks_with_a_small_display_budget(
    tmp_path, indexed
):
    masses = np.linspace(100, 1000, 1000)
    spectra = [
        _spectrum(1, rt, masses, [1] * 999 + [1000], native_id=f"scan={rt}")
        for rt in range(20)
    ]
    path = _write(tmp_path / "dense.mzML", spectra, indexed)
    points, traces = msdata.ms1_overview(path, max_points=17)
    assert len(points) <= 17
    assert points["intensity"].sum() == 20 * 1999
    assert traces[traces["trace"] == "TIC"]["intensity"].tolist() == [1999] * 20
    assert traces[traces["trace"] == "BPC"]["intensity"].tolist() == [1000] * 20


def test_scan_metadata_does_not_materialize_peaks(acquisition, monkeypatch):
    original = msdata._scan_row

    def metadata_row(index, spectrum):
        assert spectrum.size() == 0
        return original(index, spectrum)

    monkeypatch.setattr(msdata, "_scan_row", metadata_row)
    assert len(msdata.scan_index(acquisition)) == 5


@pytest.mark.parametrize("indexed", [True, False])
def test_empty_ms1_and_single_coordinate_overviews(tmp_path, indexed):
    empty = _write(tmp_path / "empty-peaks.mzML", [_spectrum(1, 2)], indexed)
    points, traces = msdata.ms1_overview(empty)
    assert points.empty
    assert traces["intensity"].tolist() == [0, 0]
    assert msdata.xic(empty, 100, 1)["intensity"].tolist() == [0]
    single = _write(
        tmp_path / "single.mzML", [_spectrum(1, 2, [100, 100], [2, 3])], indexed
    )
    points, traces = msdata.ms1_overview(single, 10)
    assert points.to_dict("records") == [{"rt": 2.0, "mz": 100.0, "intensity": 5.0}]
    assert traces["intensity"].tolist() == [5, 3]


def test_malformed_mzml_errors_explain_how_to_recover(tmp_path):
    path = tmp_path / "broken.mzML"
    path.write_text("<mzML>broken")
    for operation in [
        msdata.scan_index,
        lambda source: msdata.spectrum_peaks(source, 0),
        msdata.ms1_overview,
        lambda source: msdata.xic(source, 100, 10),
    ]:
        with pytest.raises(RuntimeError, match="complete, valid mzML"):
            operation(path)


def test_indexed_access_does_not_use_the_streaming_or_eager_parser(
    tmp_path, monkeypatch
):
    path = _write(tmp_path / "indexed.mzML", [_spectrum(1, 1, [100], [2])], True)

    def forbidden_parser():
        pytest.fail("An indexed acquisition should use native random access.")

    monkeypatch.setattr(msdata.pyopenms, "MzMLFile", forbidden_parser)
    assert msdata.scan_index(path)["index"].tolist() == [0]
    assert msdata.spectrum_peaks(path, 0)["mz"].tolist() == [100]
    assert msdata.ms1_overview(path)[0]["intensity"].tolist() == [2]
    assert msdata.xic(path, 100, 1)["intensity"].tolist() == [2]


def test_unindexed_access_streams_without_loading_an_experiment(tmp_path, monkeypatch):
    path = _write(tmp_path / "unindexed.mzML", [_spectrum(1, 1, [100], [2])], False)
    native_parser = pyopenms.MzMLFile
    passes = []

    class StreamingParser:
        def __init__(self):
            self.parser = native_parser()

        def getOptions(self):
            return self.parser.getOptions()

        def setOptions(self, options):
            self.parser.setOptions(options)

        def transform(self, *args):
            passes.append(1)
            self.parser.transform(*args)

        def load(self, *_args):
            pytest.fail("Unindexed acquisitions must be streamed, not eagerly loaded.")

    monkeypatch.setattr(msdata.pyopenms, "MzMLFile", StreamingParser)
    assert msdata.scan_index(path)["index"].tolist() == [0]
    assert msdata.spectrum_peaks(path, 0)["mz"].tolist() == [100]
    assert msdata.ms1_overview(path)[0]["intensity"].tolist() == [2]
    assert msdata.xic(path, 100, 1)["intensity"].tolist() == [2]
    assert len(passes) == 5


def test_revision_change_during_read_requires_refresh(acquisition, monkeypatch):
    original = msdata._visit_ms1

    def mutate(path, reader, consume):
        original(path, reader, consume)
        with Path(path).open("ab") as source:
            source.write(b"\n")

    monkeypatch.setattr(msdata, "_visit_ms1", mutate)
    with pytest.raises(RuntimeError, match="changed.*Refresh"):
        msdata.xic(acquisition, 100, 1)
