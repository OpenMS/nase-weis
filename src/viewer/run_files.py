"""The spectra of one acquisition (an mzML file) in the columns the library's
MS data components read (``src/viz/ms_data``)."""
from __future__ import annotations

import re
from pathlib import Path

import polars as pl

PROTON = 1.007276
_SCAN = re.compile(r"scan=(\d+)")


def precursor_mass(mz: float, charge: int) -> float:
    """Neutral mass of a precursor; RNA is measured in negative mode, so a
    negative charge adds the protons back."""
    if not charge:
        return mz
    if charge < 0:
        return (mz + PROTON) * abs(charge)
    return (mz - PROTON) * charge


def read_tables(path: Path) -> dict[str, pl.DataFrame]:
    """``scans`` (scan table, with each spectrum's native id), ``peaks``
    (every peak by spectrum), ``ms1`` (the MS1 peak map) and ``chromatogram``
    (TIC and base peak of the MS1 spectra)."""
    import numpy as np
    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    scans, peaks_mz, peaks_int, peaks_idx, chrom = [], [], [], [], []
    for index, spec in enumerate(exp.getSpectra()):
        mz, intensity = spec.get_peaks()
        level, rt = spec.getMSLevel(), spec.getRT()
        native = spec.getNativeID()
        if isinstance(native, bytes):
            native = native.decode()
        found = _SCAN.search(native or "")
        precursor = float("nan")
        if spec.getPrecursors():
            p = spec.getPrecursors()[0]
            precursor = precursor_mass(p.getMZ(), p.getCharge())
        scans.append({
            "index": index, "Scan": int(found.group(1)) if found else index, "MSLevel": level,
            "RT": rt, "PrecursorMass": precursor, "#Peaks": len(mz), "native_id": native or "",
        })
        peaks_mz.append(mz)
        peaks_int.append(intensity)
        peaks_idx.append(np.full(len(mz), index, dtype=np.int64))
        if level == 1:
            chrom.append((rt, float(intensity.sum()) if len(intensity) else 0.0,
                          float(intensity.max()) if len(intensity) else 0.0, index))
    levels = {s["index"]: s["MSLevel"] for s in scans}
    rts = {s["index"]: s["RT"] for s in scans}

    def concat(parts, dtype):
        return np.concatenate(parts).astype(dtype) if parts else np.array([], dtype=dtype)

    peaks = pl.DataFrame({
        "index": concat(peaks_idx, np.int64), "mz": concat(peaks_mz, np.float64),
        "intensity": concat(peaks_int, np.float64),
    })
    ms1 = (peaks.with_columns(pl.col("index").replace_strict(levels, default=0).alias("level"))
           .filter(pl.col("level") == 1)
           .select(pl.col("index").replace_strict(rts, default=0.0).alias("rt"),
                   pl.col("mz").alias("mass"), "intensity", pl.col("index").alias("scan_idx")))
    chromatogram = pl.DataFrame(
        [{"x": rt, "y": y, "series": name, "index": i}
         for rt, tic, bpc, i in chrom for name, y in (("TIC", tic), ("Base peak", bpc))],
        schema={"x": pl.Float64, "y": pl.Float64, "series": pl.Utf8, "index": pl.Int64},
    )
    return {"scans": pl.DataFrame(scans, schema={
        "index": pl.Int64, "Scan": pl.Int64, "MSLevel": pl.Int64, "RT": pl.Float64,
        "PrecursorMass": pl.Float64, "#Peaks": pl.Int64, "native_id": pl.Utf8,
    }), "peaks": peaks, "ms1": ms1, "chromatogram": chromatogram}


def spectra(peaks: pl.DataFrame) -> pl.DataFrame:
    """One row per spectrum in the columns of the library's raw spectrum: its
    peaks as ``MonoMass_Anno``/``SumIntensity_Anno``, empty mass lists."""
    return (peaks.group_by("index", maintain_order=True)
            .agg(pl.col("mz").alias("MonoMass_Anno"), pl.col("intensity").alias("SumIntensity_Anno"))
            .with_columns(pl.lit([], dtype=pl.List(pl.Float64)).alias("MonoMass"),
                          pl.lit([], dtype=pl.List(pl.Float64)).alias("SumIntensity")))


def tables(path: Path) -> dict[str, pl.DataFrame]:
    """``read_tables``, cached while the file is unchanged."""
    import streamlit as st

    @st.cache_data(show_spinner="Reading the spectra of this acquisition...", max_entries=4)
    def read(path: str, version: int):
        return read_tables(Path(path))

    return read(str(path), Path(path).stat().st_mtime_ns)
