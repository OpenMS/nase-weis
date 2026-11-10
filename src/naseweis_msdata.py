"""Read acquisition data without loading an entire mzML experiment into memory.

Only the overview is binned for display. Spectrum peaks, extracted ion traces,
TICs, and BPCs always use the original full-resolution data.
"""

import math
from numbers import Integral, Real
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms

SCAN_COLUMNS = ["index", "native_id", "ms_level", "rt", "precursor_mz", "charge"]
PEAK_COLUMNS = ["mz", "intensity"]
POINT_COLUMNS = ["rt", "mz", "intensity"]
TRACE_COLUMNS = ["rt", "intensity", "trace"]


def file_revision(path):
    """Return path, mtime, size, ctime, inode, and device for cache invalidation."""
    try:
        source = Path(path)
    except TypeError as exc:
        raise ValueError(
            "Choose an existing mzML file to view acquisition data."
        ) from exc
    if source.suffix.lower() != ".mzml":
        raise ValueError(
            "The acquisition viewer requires mzML. Convert this acquisition to "
            "mzML to view it; vendor files can still be used for NASE searches."
        )
    try:
        source = source.resolve(strict=True)
        if not source.is_file():
            raise OSError(f"Acquisition is not a file: {source}")
        info = source.stat()
    except OSError as exc:
        raise OSError(f"Cannot access mzML acquisition {source}: {exc}") from exc
    return (
        str(source),
        info.st_mtime_ns,
        info.st_size,
        info.st_ctime_ns,
        info.st_ino,
        info.st_dev,
    )


def _ensure_revision(revision):
    if file_revision(revision[0]) != revision:
        raise RuntimeError(
            "The mzML file changed while reading it. Refresh the viewer."
        )


def _read_error(path, exc):
    return RuntimeError(
        f"Could not read mzML acquisition {path}. Check that it is a complete, "
        f"valid mzML file: {exc}"
    )


def _indexed_reader(path, metadata=True):
    # Unindexed files are common. Avoid asking the native indexed reader to open
    # those files, since it emits a warning before returning False.
    with open(path, "rb") as source:
        source.seek(0, 2)
        source.seek(max(0, source.tell() - 65536))
        if b"indexListOffset" not in source.read():
            return None
    reader = pyopenms.OnDiscMSExperiment()
    try:
        if reader.openFile(path, not metadata):
            return reader
    except RuntimeError:
        # An index can be stale or unsupported. The streaming parser remains
        # authoritative and supplies an actionable error for malformed XML.
        pass
    return None


class _Consumer:
    """Minimal native streaming consumer; no experiment or peak list is kept."""

    def __init__(self, consume):
        self.consume = consume

    def consumeSpectrum(self, spectrum):
        """Process one spectrum while its native storage is alive."""
        self.consume(spectrum)

    def consumeChromatogram(self, _chromatogram):
        """Acquisition views derive traces from MS1 spectra."""

    def setExperimentalSettings(self, _settings):
        """Experiment-level settings are not needed for the data views."""

    def setExpectedSize(self, _spectra, _chromatograms):
        """Do not allocate for the total peak count."""


def _stream(path, consume, *, ms1_only=False, fill_data=True):
    parser = pyopenms.MzMLFile()
    options = parser.getOptions()
    options.setFillData(fill_data)
    if ms1_only:
        options.setMSLevels([1])
    parser.setOptions(options)
    try:
        parser.transform(path, _Consumer(consume), True, True)
    except (RuntimeError, OSError) as exc:
        raise _read_error(path, exc) from exc


def _visit_ms1(path, reader, consume):
    if reader is None:
        _stream(path, consume, ms1_only=True)
        return
    metadata = reader.getMetaData()
    try:
        for index in range(reader.getNrSpectra()):
            if metadata[index].getMSLevel() == 1:
                consume(reader.getSpectrum(index))
    except (RuntimeError, OSError) as exc:
        raise _read_error(path, exc) from exc


def _scan_row(index, spectrum):
    precursors = spectrum.getPrecursors()
    precursor = precursors[0] if precursors else None
    return (
        index,
        spectrum.getNativeID(),
        spectrum.getMSLevel(),
        spectrum.getRT(),
        precursor.getMZ() if precursor is not None else np.nan,
        precursor.getCharge() if precursor is not None else 0,
    )


def scan_index(path):
    """Read scan metadata, preserving acquisition order and duplicate native IDs."""
    revision = file_revision(path)
    source = revision[0]
    reader = _indexed_reader(source)
    rows = []
    if reader is not None:
        metadata = reader.getMetaData()
        rows = [_scan_row(i, metadata[i]) for i in range(reader.getNrSpectra())]
    else:
        _stream(
            source,
            lambda spectrum: rows.append(_scan_row(len(rows), spectrum)),
            fill_data=False,
        )
    _ensure_revision(revision)
    return pd.DataFrame(rows, columns=SCAN_COLUMNS)


def spectrum_peaks(path, index):
    """Return all original peaks of a scan selected by acquisition index."""
    if isinstance(index, bool) or not isinstance(index, Integral) or index < 0:
        raise ValueError("Select a non-negative integer scan index.")
    revision = file_revision(path)
    source = revision[0]
    reader = _indexed_reader(source, metadata=False)
    peaks = None
    if reader is not None:
        if index < reader.getNrSpectra():
            try:
                peaks = tuple(
                    array.copy() for array in reader.getSpectrum(int(index)).get_peaks()
                )
            except (RuntimeError, OSError) as exc:
                raise _read_error(source, exc) from exc
    else:
        scan_number = 0

        def select(spectrum):
            nonlocal scan_number, peaks
            if scan_number == index:
                peaks = tuple(array.copy() for array in spectrum.get_peaks())
            scan_number += 1

        _stream(source, select)
    if peaks is None:
        raise ValueError(
            f"Scan index {index} is outside this acquisition. Refresh the scan list."
        )
    _ensure_revision(revision)
    return pd.DataFrame(dict(zip(PEAK_COLUMNS, peaks)))


def _ms1_values(spectrum):
    rt = spectrum.getRT()
    mz, intensity = spectrum.get_peaks()
    if (
        not math.isfinite(rt)
        or not np.isfinite(mz).all()
        or not np.isfinite(intensity).all()
    ):
        raise ValueError(
            "The mzML contains non-finite MS1 values. Check or reconvert the acquisition."
        )
    return rt, mz, intensity


def _require_ms1(count):
    if not count:
        raise ValueError(
            "This acquisition contains no MS1 scans. Choose an mzML containing MS1 data."
        )


def ms1_overview(path, max_points=40000):
    """Return a bounded, summed MS1 overview and exact per-scan TIC/BPC traces.

    The overview uses at most ``max_points`` occupied rectangular RT/m/z bins.
    Coordinates are bin centers, and intensities are sums of the original peaks.
    Two streaming passes keep memory bounded without sampling away weak peaks.
    """
    if (
        isinstance(max_points, bool)
        or not isinstance(max_points, Integral)
        or max_points < 1
    ):
        raise ValueError("The overview point limit must be a positive integer.")
    revision = file_revision(path)
    source = revision[0]
    reader = _indexed_reader(source)
    traces = []
    ms1_count = peak_count = 0
    rt_min = mz_min = math.inf
    rt_max = mz_max = -math.inf

    def summarize(spectrum):
        nonlocal ms1_count, peak_count, rt_min, rt_max, mz_min, mz_max
        rt, mz, intensity = _ms1_values(spectrum)
        ms1_count += 1
        peak_count += len(mz)
        rt_min, rt_max = min(rt_min, rt), max(rt_max, rt)
        traces.extend(
            [
                (rt, float(np.sum(intensity, dtype=np.float64)), "TIC"),
                (rt, float(np.max(intensity)) if len(intensity) else 0.0, "BPC"),
            ]
        )
        if len(mz):
            mz_min, mz_max = (
                min(mz_min, float(np.min(mz))),
                max(mz_max, float(np.max(mz))),
            )

    _visit_ms1(source, reader, summarize)
    _require_ms1(ms1_count)
    trace_frame = pd.DataFrame(traces, columns=TRACE_COLUMNS)
    if not peak_count:
        _ensure_revision(revision)
        return pd.DataFrame(columns=POINT_COLUMNS), trace_frame

    budget = min(int(max_points), peak_count)
    if rt_min == rt_max:
        rt_bins, mz_bins = 1, budget
    elif mz_min == mz_max:
        rt_bins, mz_bins = min(ms1_count, budget), 1
    else:
        rt_bins = min(ms1_count, max(1, math.isqrt(budget)))
        mz_bins = max(1, budget // rt_bins)
    sums = np.zeros((rt_bins, mz_bins), dtype=np.float64)
    occupied = np.zeros((rt_bins, mz_bins), dtype=bool)

    def add(spectrum):
        rt, mz, intensity = _ms1_values(spectrum)
        rt_bin = (
            0
            if rt_min == rt_max
            else min(int((rt - rt_min) / (rt_max - rt_min) * rt_bins), rt_bins - 1)
        )
        if mz_min == mz_max:
            mz_bin = np.zeros(len(mz), dtype=np.intp)
        else:
            mz_bin = np.clip(
                ((mz - mz_min) / (mz_max - mz_min) * mz_bins).astype(np.intp),
                0,
                mz_bins - 1,
            )
        np.add.at(sums[rt_bin], mz_bin, intensity)
        occupied[rt_bin, mz_bin] = True

    _visit_ms1(source, reader, add)
    _ensure_revision(revision)
    rt_cells, mz_cells = np.nonzero(occupied)
    rt_centers = rt_min + (rt_cells + 0.5) / rt_bins * (rt_max - rt_min)
    mz_centers = mz_min + (mz_cells + 0.5) / mz_bins * (mz_max - mz_min)
    return pd.DataFrame(
        {"rt": rt_centers, "mz": mz_centers, "intensity": sums[rt_cells, mz_cells]},
        columns=POINT_COLUMNS,
    ), trace_frame


def _finite_number(value, name, *, positive=False, nonnegative=True):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number.")  # noqa: TRY004 - viewer validation uses ValueError
    try:
        number = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number.")
    if (nonnegative and number < 0) or (positive and number == 0):
        raise ValueError(
            f"{name} must be {'positive' if positive else 'non-negative'}."
        )
    return number


def xic(path, mz, tolerance, unit="ppm", rt_range=None):
    """Sum full-resolution MS1 peaks inclusively inside a ppm or Da window.

    Every MS1 scan in the optional inclusive RT window contributes a row,
    including scans with no peaks inside the mass window (intensity zero).
    """
    center = _finite_number(mz, "XIC m/z", positive=True)
    tolerance = _finite_number(tolerance, "XIC tolerance")
    if not isinstance(unit, str) or unit.lower() not in {"ppm", "da"}:
        raise ValueError("XIC tolerance unit must be ppm or Da.")
    width = center * (tolerance * 1e-6) if unit.lower() == "ppm" else tolerance
    if not math.isfinite(width) or not math.isfinite(center + width):
        raise ValueError("XIC mass window must have finite bounds.")
    if rt_range is not None:
        try:
            start, end = rt_range
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "XIC RT range must contain a start and end in seconds."
            ) from exc
        start = _finite_number(start, "XIC RT range", nonnegative=False)
        end = _finite_number(end, "XIC RT range", nonnegative=False)
        if start > end:
            raise ValueError(
                "XIC RT range must have finite, increasing bounds in seconds."
            )
    revision = file_revision(path)
    source = revision[0]
    reader = _indexed_reader(source)
    rows = []
    ms1_count = 0

    def extract(spectrum):
        nonlocal ms1_count
        rt, masses, intensities = _ms1_values(spectrum)
        ms1_count += 1
        if rt_range is None or start <= rt <= end:
            selected = (masses >= center - width) & (masses <= center + width)
            rows.append(
                (rt, float(np.sum(intensities[selected], dtype=np.float64)), "XIC")
            )

    _visit_ms1(source, reader, extract)
    _require_ms1(ms1_count)
    _ensure_revision(revision)
    return pd.DataFrame(rows, columns=TRACE_COLUMNS)
