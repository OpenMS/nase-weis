"""The components of Explore results, one ``render(ctx, src)`` each.

Two families: the result sections NASE-weis has always shown (summary,
matches, the spectrum viewer, coverage, quantification, downloads), drawn by
the workflow's own section methods, and OpenMS-Insight components from the
visualization library (``src/viz``), linked through the experiment's state so
a picked spectrum or oligonucleotide stays picked across components.

A result section uses Streamlit widget keys of its run, so it is drawn once per
page: in the first experiment, and not again when "All results" already holds
it (``claim``).
"""
from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import streamlit as st

from src.naseweis_viewers import viewer_key
from src.viewer import run_files
from src.viewer.sources import Source

SPECTRUM, PSM, OLIGO = "scan", "psm", "oligo"
_drawn: set[str] = set()

# The result sections the "All results" component holds
SECTIONS_SINGLE = ("summary", "matches", "spectrum_viewer", "coverage", "quantification", "downloads")
SECTIONS_EXPERIMENT = ("experiment_summary", "experiment_coverage", "experiment_quantification", "sample_details")


def reset() -> None:
    """Forget what the last rerun drew (called at the top of the page)."""
    _drawn.clear()


def claim(ctx, *names: str) -> bool:
    """Whether result sections ``names`` may be drawn now; says why not."""
    if ctx.index:
        st.info("Result sections are shown for the first run. Pick this run above to see them.")
        return False
    if _drawn.intersection(names):
        st.caption("Shown in All results.")
        return False
    _drawn.update(names)
    return True


def key(ctx, src: Source, name: str) -> str:
    return f"nase-{name}-{viewer_key(src.dataset)}-e{ctx.index}"


def zoom(cls, name: str) -> dict:
    from src.layout_manager import zoom_identifier

    if "zoom_identifier" not in inspect.signature(cls.__init__).parameters:
        return {}
    return {"zoom_identifier": zoom_identifier(name)}


def _loaded(src: Source):
    cached = getattr(src, "_loaded", None)
    if cached is None:
        cached = src.wf.load_single()
        src._loaded = cached
    return cached


# --- result sections -------------------------------------------------------

def all_results(ctx, src: Source) -> None:
    """Every result section, in tabs, as Explore results always showed them."""
    from src.naseweis_experiment_ui import experiment_context

    try:
        experiment = src.sample is None and experiment_context(src.wf) is not None
    except (OSError, ValueError):
        experiment = src.sample is None
    names = SECTIONS_EXPERIMENT if experiment else SECTIONS_SINGLE
    if not claim(ctx, *names):
        return
    if src.sample is None:
        src.wf.results()
    else:
        src.wf._single_results()


def _single(section: str):
    def render(ctx, src: Source) -> None:
        if not claim(ctx, section):
            return
        loaded = _loaded(src)
        if loaded is not None:
            getattr(src.wf, f"{section}_section")(loaded)
    return render


summary = _single("summary")
matches = _single("matches")
coverage = _single("coverage")


def downloads(ctx, src: Source) -> None:
    if not claim(ctx, "downloads"):
        return
    loaded = _loaded(src)
    if loaded is not None:
        src.wf.export_section(loaded)


def spectrum_viewer(ctx, src: Source) -> None:
    if not claim(ctx, "spectrum_viewer"):
        return
    loaded = _loaded(src)
    if loaded is not None:
        src.wf.spectra_section(loaded)


def quantification(ctx, src: Source) -> None:
    if claim(ctx, "quantification"):
        src.wf.quantification_results()


def _context(src: Source):
    from src.naseweis_experiment_ui import experiment_context

    try:
        return experiment_context(src.wf)
    except (OSError, ValueError) as error:
        st.error(f"Could not read this experiment: {error}")
        return None


def _experiment(section: str):
    def render(ctx, src: Source) -> None:
        if not claim(ctx, section):
            return
        context = _context(src)
        if context is not None:
            from src import naseweis_experiment_ui as ui

            getattr(ui, f"{section}_section")(context)
    return render


experiment_summary = _experiment("experiment_summary")
experiment_coverage = _experiment("experiment_coverage")
experiment_quantification = _experiment("experiment_quantification")


def sample_details(ctx, src: Source) -> None:
    if not claim(ctx, "sample_details"):
        return
    context = _context(src)
    if context is not None:
        from src.naseweis_experiment_ui import sample_details_section

        sample_details_section(context, src.wf._render_sample_results)


# --- OpenMS-Insight components: one acquisition ------------------------------

def _spectra(ctx, src: Source, title: str):
    path = src.mzml()
    if path is None:
        st.markdown(f"**{title}**")
        st.info("The spectra of this acquisition are not available: " + "; ".join(src.notes))
        return None
    return path


def _version(path: Path) -> str:
    return str(Path(path).stat().st_mtime_ns)


def file_scans(ctx, src: Source) -> None:
    from openms_insight import Table

    from src.viz.ms_data import scan_table
    from src.viz.plumbing import component

    path = _spectra(ctx, src, "Scans")
    if path is None:
        return
    table = component(Table, src.cache_dir, f"scans_{_version(path)}",
                      lambda: {"data": run_files.tables(path)["scans"].drop("native_id").lazy()},
                      **scan_table.settings(SPECTRUM, clears=(), peaks=True, title="Scans"))
    table(state_manager=ctx.state, height=376, key=key(ctx, src, "scans"))


def file_spectrum(ctx, src: Source) -> None:
    from openms_insight import DeconvolvedSpectrum

    from src.common.plot_style import spectrum_styling
    from src.viz.ms_data import spectrum
    from src.viz.plumbing import component, only_accepted

    path = _spectra(ctx, src, "Spectrum")
    if path is None:
        return
    if ctx.state.get_selection(SPECTRUM) is None:
        st.markdown("**Spectrum**")
        st.info("Pick a spectrum in the scan table, the matches, the chromatogram or the peak map.")
        return
    settings = spectrum.settings(SPECTRUM, styling=spectrum_styling(), signal_column=None, title="Spectrum")
    plot = component(DeconvolvedSpectrum, src.cache_dir, f"spectra_{_version(path)}",
                     lambda: {"data": run_files.spectra(run_files.tables(path)["peaks"]).lazy()},
                     **only_accepted(DeconvolvedSpectrum, **settings),
                     **zoom(DeconvolvedSpectrum, "file_spectrum"))
    plot(state_manager=ctx.state, height=360, key=key(ctx, src, "spectrum"))


def file_chromatogram(ctx, src: Source) -> None:
    from openms_insight import SeriesPlot

    from src.viz.ms_data import chromatogram
    from src.viz.plumbing import component, only_accepted

    path = _spectra(ctx, src, "Chromatogram (MS1)")
    if path is None:
        return
    data = run_files.tables(path)["chromatogram"]
    if data.is_empty():
        st.markdown("**Chromatogram (MS1)**")
        st.info("This acquisition holds no MS1 spectra.")
        return
    settings = chromatogram.settings(["TIC", "Base peak"], level=1, spectrum=SPECTRUM)
    plot = component(SeriesPlot, src.cache_dir, f"chromatogram_{_version(path)}",
                     lambda: {"data": data.lazy()},
                     **only_accepted(SeriesPlot, **settings), **zoom(SeriesPlot, "file_chromatogram"))
    plot(state_manager=ctx.state, height=320, key=key(ctx, src, "chromatogram"))


def file_peak_map(ctx, src: Source) -> None:
    from openms_insight import Heatmap

    from src.viz.ms_data import raw_heatmap
    from src.viz.plumbing import component, only_accepted

    path = _spectra(ctx, src, "Peak map (MS1)")
    if path is None:
        return
    data = run_files.tables(path)["ms1"]
    if data.is_empty():
        st.markdown("**Peak map (MS1)**")
        st.info("This acquisition holds no MS1 spectra.")
        return
    settings = raw_heatmap.settings(title="Peak map (MS1)", spectrum=SPECTRUM,
                                    zoom_identifier=f"file_peak_map_zoom_{viewer_key(src.dataset)}")
    plot = component(Heatmap, src.cache_dir, f"ms1_{_version(path)}", lambda: {"data": data.lazy()},
                     **only_accepted(Heatmap, **settings))
    plot(state_manager=ctx.state, height=420, key=key(ctx, src, "peak_map"))


def psm_frame(tables: dict, matches: pd.DataFrame, scans: pl.DataFrame | None) -> pl.DataFrame:
    """The spectrum matches of an acquisition in the library's PSM columns;
    ``scan_id`` is the position of the matched spectrum in the mzML (-1 when
    it is not found there)."""
    from src.naseweis_annotations import native_spectrum_id

    osm = tables.get("OSM", pd.DataFrame())
    if osm.empty:
        return pl.DataFrame()
    positions = {}
    if scans is not None:
        positions = dict(zip(scans["native_id"].to_list(), scans["index"].to_list()))

    def position(ref):
        try:
            return positions.get(native_spectrum_id(ref), -1)
        except ValueError:
            return -1

    accession = matches["Accession"] if "Accession" in matches else pd.Series([""] * len(osm))
    frame = pd.DataFrame({
        "id_idx": np.arange(len(osm), dtype=np.int64),
        "file_index": 0,
        "scan_id": [position(r) for r in osm.get("spectra_ref", [""] * len(osm))],
        "sequence": osm["sequence"].astype(str).to_numpy(),
        "charge": pd.to_numeric(osm.get("charge"), errors="coerce").fillna(0).astype(int).to_numpy(),
        "mz": pd.to_numeric(osm.get("exp_mass_to_charge"), errors="coerce").to_numpy(),
        "rt": pd.to_numeric(osm.get("retention_time"), errors="coerce").to_numpy(),
        "score": pd.to_numeric(osm.get("search_engine_score[1]"), errors="coerce").to_numpy(),
        "decoy": False,
        "protein_accession": accession.fillna("").astype(str).to_numpy(),
        "contaminant": False,
        "proforma": osm["sequence"].astype(str).to_numpy(),
        "filename": "",
    })
    return pl.from_pandas(frame)


def match_table(ctx, src: Source) -> None:
    """The spectrum matches as an Insight table; a click picks the spectrum."""
    from openms_insight import Table

    from src.viz.identification import psm_table
    from src.viz.plumbing import component

    loaded = _loaded(src)
    if loaded is None:
        return
    path = src.mzml()
    scans = run_files.tables(path)["scans"] if path is not None else None
    data = psm_frame(loaded["tables"], loaded["matches"], scans)
    if data.is_empty():
        st.markdown("**Spectrum matches**")
        st.info("No matches passed this search's filters.")
        return
    columns = [dict(c, title="Oligonucleotide") if c.get("field") == "sequence" else
               dict(c, title="RNA accessions") if c.get("field") == "protein_accession" else c
               for c in psm_table.columns() if c.get("field") not in ("decoy", "contaminant")]
    stamp = str(Path(loaded["file"]).stat().st_mtime_ns)
    table = component(Table, src.cache_dir, f"psms_{stamp}_{path is not None}", lambda: {"data": data.lazy()},
                      **psm_table.settings(ascending=False, file="file", spectrum=SPECTRUM, psm=PSM,
                                           column_definitions=columns, title="Spectrum matches"))
    table(state_manager=ctx.state, height=376, key=key(ctx, src, "psms"))


def _histogram(values, bins: int = 30) -> pd.DataFrame:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return pd.DataFrame(columns=["x", "y"])
    counts, edges = np.histogram(values, bins=bins)
    return pd.DataFrame({"x": (edges[:-1] + edges[1:]) / 2, "y": counts.astype(float)})


def _series(ctx, src, name, data: pd.DataFrame, settings: dict) -> None:
    from openms_insight import SeriesPlot

    from src.viz.plumbing import fingerprint, only_accepted

    frame = pl.from_pandas(data)
    SeriesPlot(cache_id=f"{name}-{fingerprint(frame, sorted(settings.items(), key=lambda kv: kv[0]))}",
               cache_path=str(src.cache_dir), data=frame.lazy(), **zoom(SeriesPlot, name),
               **only_accepted(SeriesPlot, **settings))(
        key=key(ctx, src, name), state_manager=ctx.state, height=320)


def score_distribution(ctx, src: Source) -> None:
    from src.viz.identification import score_distribution as viz

    loaded = _loaded(src)
    if loaded is None:
        return
    osm = loaded["tables"].get("OSM", pd.DataFrame())
    data = _histogram(osm.get("search_engine_score[1]", []))
    if data.empty:
        st.markdown("**Search scores**")
        st.info("No matches passed this search's filters.")
        return
    data["kind"] = "target"
    _series(ctx, src, "score_distribution", data,
            viz.settings(1, title="Search scores", x_label="Search score", y_label="Matches",
                         fixed_y_range=False))


def precursor_error(ctx, src: Source) -> None:
    from src.viz.identification import precursor_mass_error as viz

    loaded = _loaded(src)
    if loaded is None:
        return
    osm = loaded["tables"].get("OSM", pd.DataFrame())
    if osm.empty:
        st.markdown("**Precursor mass error**")
        st.info("No matches passed this search's filters.")
        return
    exp = pd.to_numeric(osm.get("exp_mass_to_charge"), errors="coerce")
    calc = pd.to_numeric(osm.get("calc_mass_to_charge"), errors="coerce")
    error = (exp - calc) / calc * 1e6
    offset = pd.to_numeric(osm.get("opt_global_isotope_offset", 0), errors="coerce")
    if isinstance(offset, pd.Series):
        error = error[offset.fillna(0).eq(0)]
    data = _histogram(error)
    if data.empty:
        st.markdown("**Precursor mass error**")
        st.info("No precursor errors in this acquisition.")
        return
    st.caption(f"Median {np.nanmedian(error):+.2f} ppm, matches without an isotope offset.")
    _series(ctx, src, "precursor_error", data, viz.settings(noun="matches", title="Precursor mass error"))


# --- OpenMS-Insight components: an experiment ------------------------------

def _comparison(src: Source):
    """(abundances, groups) of a finished comparison; None with a note."""
    from src.naseweis_comparison import load_comparison
    from src.naseweis_experiment_ui import comparison_ready

    context = _context(src)
    if context is None:
        return None
    if not comparison_ready(context):
        st.info("This experiment has no finished comparison of its samples.")
        return None
    try:
        loaded = load_comparison(src.wf.results_dir)
    except (OSError, ValueError) as error:
        st.error(f"Could not read this experiment's comparison: {error}")
        return None
    design = loaded["design"]
    groups = {s: (str(design.loc[s, "condition"]) or "Unassigned") for s in design.index}
    names = {s: str(design.loc[s, "file"]) for s in design.index}
    return loaded["abundances"], groups, names


def _abundance_rows(abundances: pd.DataFrame) -> pd.DataFrame:
    frame = abundances.reset_index().rename(columns={abundances.index.name or "index": "label"})
    frame.insert(0, "oligo_idx", np.arange(len(frame), dtype=np.int64))
    return frame


def abundance_table(ctx, src: Source) -> None:
    """Normalized abundance of every oligonucleotide per sample; a click picks it."""
    from openms_insight import Table

    from src.viz.plumbing import fingerprint

    loaded = _comparison(src)
    if loaded is None:
        return
    abundances, _, names = loaded
    frame = pl.from_pandas(_abundance_rows(abundances))
    columns = [{"field": "label", "title": "RNA oligonucleotide"}] + [
        {"field": s, "title": names.get(s, s), "sorter": "number",
         "formatter": "money", "formatterParams": {"thousand": False, "precision": 4}}
        for s in abundances.columns
    ]
    Table(cache_id=f"abundances-{fingerprint(frame)}", cache_path=str(src.cache_dir), data=frame.lazy(),
          column_definitions=columns, index_field="oligo_idx", interactivity={OLIGO: "oligo_idx"},
          title="Abundances")(state_manager=ctx.state, height=376, key=key(ctx, src, "abundances"))


def abundance_profile(ctx, src: Source) -> None:
    from src.viz.quantification.abundance_profile import show_profile

    loaded = _comparison(src)
    if loaded is None:
        return
    abundances, groups, names = loaded
    rows = _abundance_rows(abundances)
    long = rows.melt(id_vars=["oligo_idx", "label"], var_name="sample", value_name="value")
    long["group"] = long["sample"].map(groups)
    long["sample"] = long["sample"].map(lambda s: names.get(s, s))
    if ctx.state.get_selection(OLIGO) is None:
        st.markdown("**Abundance profile**")
        st.info("Pick an oligonucleotide in the abundance table.")
        return
    show_profile(long, "oligo_idx", OLIGO, "nase-abundance-profile", cache_path=src.cache_dir,
                 key=key(ctx, src, "profile"), state_manager=ctx.state, height=420,
                 table_args={"hover_column": "label"}, title="Abundance profile",
                 y_label="Normalized abundance")


def zscore_heatmap(ctx, src: Source) -> None:
    from src.viz.quantification.zscore_heatmap import show_zscore_heatmap

    loaded = _comparison(src)
    if loaded is None:
        return
    abundances, groups, names = loaded
    matrix = np.log2(abundances.where(abundances > 0))
    matrix.index = matrix.index.astype(str)
    matrix.index.name = "label"
    matrix = matrix.rename(columns=names)
    shown = show_zscore_heatmap(matrix, "label", "nase-zscore", cache_path=src.cache_dir,
                                key=key(ctx, src, "zscores"), state_manager=ctx.state, top=50,
                                groups={names.get(s, s): g for s, g in groups.items()}, height=420,
                                title="Abundance z-scores", y_label="RNA oligonucleotides")
    if shown is None:
        st.markdown("**Abundance z-scores**")
        st.info("No oligonucleotide differs between the samples.")
