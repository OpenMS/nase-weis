"""The components of Explore results, one ``render(ctx, src)`` each.

Two families: the result sections NASE-weis has always shown (summary,
matches, the spectrum viewer, coverage, quantification, downloads), drawn by
the workflow's own section methods, and OpenMS-Insight components from the
visualization library (``src/viz``), linked through the experiment's state so
a picked spectrum or oligonucleotide stays picked across components.

A result section uses Streamlit widget keys of its run, so it is drawn once per
page: in the first experiment, and once when a view holds it twice (``claim``).
"""
from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import streamlit as st

from src.naseweis_viewers import viewer_key
from src.viewer import run_files, teaching
from src.viewer.sources import Source

SPECTRUM, PSM, OLIGO, RNA = "scan", "psm", "oligo", "rna"
# The picked oligonucleotide (by sequence) and the picked peak of a fragment spectrum
OLIGO_ID, FRAGMENT = "oligo_id", "fragment_peak"
_drawn: set[str] = set()

# The result sections of a run and of an experiment
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
        st.caption("Shown above in this view.")
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

TITLES = {"matches": "Matches", "spectrum_viewer": "Spectra", "quantification": "Quantification",
          "downloads": "Export", "experiment_summary": "Experiment summary",
          "experiment_quantification": "Quantification", "sample_details": "Sample details"}


def _titled(section: str, draw) -> None:
    """Draw a result section under the heading Explore results always gave it."""
    title = TITLES.get(section)
    if title is None:
        draw()
        return
    with st.container(border=True):
        st.subheader(title, help=teaching.component(section))
        draw()


def _heading(title: str, name: str) -> None:
    """The title an Insight component would draw, for a component that shows a note instead."""
    st.markdown(f"**{title}**", help=teaching.component(name))


def _single(section: str):
    def render(ctx, src: Source) -> None:
        if not claim(ctx, section):
            return
        loaded = _loaded(src)
        if loaded is not None:
            _titled(section, lambda: getattr(src.wf, f"{section}_section")(loaded))
    return render


summary = _single("summary")
matches = _single("matches")


# --- identification: the RNA table and the coverage map it picks -----------

def _rna_groups(src: Source):
    """``[(group, mzTab tables)]`` of the shown run: the acquisition, or every
    acquisition of an experiment; None when there are none to read."""
    cached = getattr(src, "_rna_groups", None)
    if cached is not None:
        return cached
    if src.sample is None and _is_experiment(src):
        from src.naseweis_result_plots_ui import experiment_coverage_tables

        context = _context(src)
        groups = experiment_coverage_tables(context["root"], context["samples"]) if context else []
    else:
        loaded = _loaded(src)
        groups = [("This acquisition", loaded["tables"])] if loaded is not None else []
    src._rna_groups = groups
    return groups


def _is_experiment(src: Source) -> bool:
    from src.viewer import sources

    return sources.kind_of(src.wf, src.dataset) == sources.EXPERIMENT


def _rnas(src: Source) -> pd.DataFrame:
    cached = getattr(src, "_rnas", None)
    if cached is None:
        from src.naseweis_result_plots_ui import rna_summary

        cached = rna_summary(_rna_groups(src))
        src._rnas = cached
    return cached


def picked_rna(ctx, src: Source) -> str | None:
    """The accession of the RNA picked in the RNA table, None while none is."""
    picked = ctx.state.get_selection(RNA)
    rows = _rnas(src)
    if picked is None or rows.empty:
        return None
    match = rows.loc[rows["rna_idx"] == int(picked), "accession"]
    return str(match.iloc[0]) if len(match) else None


def rna_table(ctx, src: Source) -> None:
    """The identified RNAs as an Insight table; a click shows that RNA's coverage."""
    from openms_insight import Table

    from src.viz.plumbing import fingerprint

    rows = _rnas(src)
    if rows.empty:
        _heading("Identified RNAs", "rna_table")
        st.info("No identified oligonucleotide maps to an RNA of the database.")
        return
    frame = pl.from_pandas(rows)
    columns = [
        {"field": "accession", "title": "RNA"},
        {"field": "length", "title": "Length", "sorter": "number"},
        {"field": "oligonucleotides", "title": "Oligonucleotides", "sorter": "number"},
        {"field": "matches", "title": "Spectrum matches", "sorter": "number"},
        {"field": "coverage", "title": "Coverage (%)", "sorter": "number"},
    ]
    Table(cache_id=f"rnas-{fingerprint(frame)}", cache_path=str(src.cache_dir), data=frame.lazy(),
          column_definitions=columns, index_field="rna_idx", interactivity={RNA: "rna_idx"},
          title="Identified RNAs")(state_manager=ctx.state, height=300, key=key(ctx, src, "rnas"),
                                   help=teaching.component("rna_table"))


def coverage(ctx, src: Source) -> None:
    if not claim(ctx, "coverage"):
        return
    loaded = _loaded(src)
    if loaded is None:
        return
    from src.naseweis_result_plots_ui import render_coverage_groups

    render_coverage_groups([("This acquisition", loaded["tables"])], f"coverage-{src.wf.results_dir}",
                           picked_rna(ctx, src), help=teaching.component("coverage"))


def downloads(ctx, src: Source) -> None:
    if not claim(ctx, "downloads"):
        return
    loaded = _loaded(src)
    if loaded is not None:
        _titled("downloads", lambda: src.wf.export_section(loaded))


def spectrum_viewer(ctx, src: Source) -> None:
    if not claim(ctx, "spectrum_viewer"):
        return
    loaded = _loaded(src)
    if loaded is not None:
        _titled("spectrum_viewer", lambda: src.wf.spectra_section(loaded))


def quantification(ctx, src: Source) -> None:
    if claim(ctx, "quantification"):
        _titled("quantification", src.wf.quantification_results)


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

            _titled(section, lambda: getattr(ui, f"{section}_section")(context))
    return render


experiment_summary = _experiment("experiment_summary")


def experiment_coverage(ctx, src: Source) -> None:
    if not claim(ctx, "experiment_coverage"):
        return
    context = _context(src)
    if context is not None:
        from src.naseweis_result_plots_ui import render_experiment_coverage

        render_experiment_coverage(context["root"], context["samples"], context["prefix"],
                                   picked_rna(ctx, src), help=teaching.component("experiment_coverage"))


def experiment_quantification(ctx, src: Source) -> None:
    if not claim(ctx, "experiment_quantification"):
        return
    context = _context(src)
    if context is not None:
        from src.naseweis_experiment_ui import experiment_quantification_section

        _titled("experiment_quantification",
                lambda: experiment_quantification_section(context, controls=False))


def sample_details(ctx, src: Source) -> None:
    if not claim(ctx, "sample_details"):
        return
    context = _context(src)
    if context is not None:
        from src.naseweis_experiment_ui import sample_details_section

        _titled("sample_details", lambda: sample_details_section(context, src.wf._render_sample_results))


# --- OpenMS-Insight components: one acquisition ------------------------------

def _spectra(ctx, src: Source, title: str, name: str):
    path = src.mzml()
    if path is None:
        _heading(title, name)
        st.info("The spectra of this acquisition are not available: " + "; ".join(src.notes))
        return None
    if _readable(path) is None:
        _heading(title, name)
        st.error(f"Could not read the spectra of {Path(path).name}; the file may be damaged or not mzML.")
        return None
    return path


def _readable(path):
    """The acquisition's tables, or None when the file cannot be read."""
    try:
        return run_files.tables(path)
    except Exception:  # pylint: disable=broad-except  # pyOpenMS raises plain RuntimeError
        return None


def _version(path: Path) -> str:
    return str(Path(path).stat().st_mtime_ns)


def file_scans(ctx, src: Source) -> None:
    from openms_insight import Table

    from src.viz.ms_data import scan_table
    from src.viz.plumbing import component

    path = _spectra(ctx, src, "Scans", "file_scans")
    if path is None:
        return
    table = component(Table, src.cache_dir, f"scans_{_version(path)}",
                      lambda: {"data": run_files.tables(path)["scans"].drop("native_id").lazy()},
                      **scan_table.settings(SPECTRUM, clears=(), peaks=True, title="Scans"))
    table(state_manager=ctx.state, height=376, key=key(ctx, src, "scans"), help=teaching.component("file_scans"))


def file_spectrum(ctx, src: Source) -> None:
    from openms_insight import DeconvolvedSpectrum

    from src.common.plot_style import spectrum_styling
    from src.viz.ms_data import spectrum
    from src.viz.plumbing import component, only_accepted

    path = _spectra(ctx, src, "Spectrum", "file_spectrum")
    if path is None:
        return
    if ctx.state.get_selection(SPECTRUM) is None:
        _heading("Spectrum", "file_spectrum")
        st.info("Pick a spectrum in the scan table, the matches, the chromatogram or the peak map.")
        return
    settings = spectrum.settings(SPECTRUM, styling=spectrum_styling(), signal_column=None, title="Spectrum")
    plot = component(DeconvolvedSpectrum, src.cache_dir, f"spectra_{_version(path)}",
                     lambda: {"data": run_files.spectra(run_files.tables(path)["peaks"]).lazy()},
                     **only_accepted(DeconvolvedSpectrum, **settings),
                     **zoom(DeconvolvedSpectrum, "file_spectrum"))
    plot(state_manager=ctx.state, height=360, key=key(ctx, src, "spectrum"), help=teaching.component("file_spectrum"))


def file_chromatogram(ctx, src: Source) -> None:
    from openms_insight import SeriesPlot

    from src.viz.ms_data import chromatogram
    from src.viz.plumbing import component, only_accepted

    path = _spectra(ctx, src, "Chromatogram (MS1)", "file_chromatogram")
    if path is None:
        return
    data = run_files.tables(path)["chromatogram"]
    if data.is_empty():
        _heading("Chromatogram (MS1)", "file_chromatogram")
        st.info("This acquisition holds no MS1 spectra.")
        return
    settings = chromatogram.settings(["TIC", "Base peak"], level=1, spectrum=SPECTRUM)
    plot = component(SeriesPlot, src.cache_dir, f"chromatogram_{_version(path)}",
                     lambda: {"data": data.lazy()},
                     **only_accepted(SeriesPlot, **settings), **zoom(SeriesPlot, "file_chromatogram"))
    plot(state_manager=ctx.state, height=320, key=key(ctx, src, "chromatogram"), help=teaching.component("file_chromatogram"))


def file_peak_map(ctx, src: Source) -> None:
    from openms_insight import Heatmap

    from src.viz.ms_data import raw_heatmap
    from src.viz.plumbing import component, only_accepted

    path = _spectra(ctx, src, "Peak map (MS1)", "file_peak_map")
    if path is None:
        return
    data = run_files.tables(path)["ms1"]
    if data.is_empty():
        _heading("Peak map (MS1)", "file_peak_map")
        st.info("This acquisition holds no MS1 spectra.")
        return
    settings = raw_heatmap.settings(title="Peak map (MS1)", spectrum=SPECTRUM,
                                    zoom_identifier=f"file_peak_map_zoom_{viewer_key(src.dataset)}")
    plot = component(Heatmap, src.cache_dir, f"ms1_{_version(path)}", lambda: {"data": data.lazy()},
                     **only_accepted(Heatmap, **settings))
    plot(state_manager=ctx.state, height=420, key=key(ctx, src, "peak_map"), help=teaching.component("file_peak_map"))


def psm_frame(tables: dict, matches: pd.DataFrame, scans: pl.DataFrame | None,
              oligo_ids: dict | None = None) -> pl.DataFrame:
    """The spectrum matches of an acquisition in the library's PSM columns;
    ``scan_id`` is the position of the matched spectrum in the mzML (-1 when
    it is not found there), ``oligo_id`` the oligonucleotide's number in
    ``oligo_ids`` (-1 without one) and ``native_id`` the spectrum's native ID."""
    from src.naseweis_annotations import native_spectrum_id

    osm = tables.get("OSM", pd.DataFrame())
    if osm.empty:
        return pl.DataFrame()
    positions = {}
    if scans is not None:
        positions = dict(zip(scans["native_id"].to_list(), scans["index"].to_list()))

    def native(ref):
        try:
            return native_spectrum_id(ref)
        except ValueError:
            return ""

    natives = [native(r) for r in osm.get("spectra_ref", [""] * len(osm))]
    sequences = osm["sequence"].astype(str)
    accession = matches["Accession"] if "Accession" in matches else pd.Series([""] * len(osm))
    frame = pd.DataFrame({
        "id_idx": np.arange(len(osm), dtype=np.int64),
        "file_index": 0,
        "scan_id": [positions.get(n, -1) if n else -1 for n in natives],
        "sequence": sequences.to_numpy(),
        "charge": pd.to_numeric(osm.get("charge"), errors="coerce").fillna(0).astype(int).to_numpy(),
        "mz": pd.to_numeric(osm.get("exp_mass_to_charge"), errors="coerce").to_numpy(),
        "rt": pd.to_numeric(osm.get("retention_time"), errors="coerce").to_numpy(),
        "score": pd.to_numeric(osm.get("search_engine_score[1]"), errors="coerce").to_numpy(),
        "decoy": False,
        "protein_accession": accession.fillna("").astype(str).to_numpy(),
        "contaminant": False,
        "proforma": sequences.to_numpy(),
        "filename": "",
        "oligo_id": np.array([(oligo_ids or {}).get(q, -1) for q in sequences], dtype=np.int64),
        "native_id": natives,
    })
    return pl.from_pandas(frame)


def _psms(src: Source):
    """``(matches, mzML tables or None, mzTab stamp)`` of the shown acquisition, None without results."""
    cached = getattr(src, "_psms", None)
    if cached is not None:
        return cached
    loaded = _loaded(src)
    if loaded is None:
        return None
    path = src.mzml()
    tables = _readable(path) if path is not None else None
    scans = tables["scans"] if tables is not None else None
    data = psm_frame(loaded["tables"], loaded["matches"], scans, _oligo_ids(src))
    src._psms = (data, tables, f"{Path(loaded['file']).stat().st_mtime_ns}_{tables is not None}")
    return src._psms


def _match_columns() -> list[dict]:
    from src.viz.identification import psm_table

    return [dict(c, title="Oligonucleotide") if c.get("field") == "sequence" else
            dict(c, title="RNA accessions") if c.get("field") == "protein_accession" else c
            for c in psm_table.columns() if c.get("field") not in ("decoy", "contaminant")]


def _match_table(ctx, src: Source, name: str, title: str, **extra) -> None:
    from openms_insight import Table

    from src.viz.identification import psm_table
    from src.viz.plumbing import component

    psms = _psms(src)
    if psms is None:
        return
    data, _tables, stamp = psms
    if data.is_empty():
        _heading(title, name)
        st.info("No matches passed this search's filters.")
        return
    table = component(Table, src.cache_dir, f"{name}_v2_{stamp}", lambda: {"data": data.lazy()},
                      **psm_table.settings(ascending=False, file="file", spectrum=SPECTRUM, psm=PSM,
                                           column_definitions=_match_columns(), title=title,
                                           on_select={FRAGMENT: None}, **extra))
    table(state_manager=ctx.state, height=376, key=key(ctx, src, name), help=teaching.component(name))


def match_table(ctx, src: Source) -> None:
    """Every spectrum match as an Insight table; a click picks the match and its spectrum."""
    _match_table(ctx, src, "match_table", "Spectrum matches")


def oligo_matches(ctx, src: Source) -> None:
    """The spectrum matches of the oligonucleotide picked in the oligonucleotide table."""
    if ctx.state.get_selection(OLIGO_ID) is None:
        _heading("Matches of the oligonucleotide", "oligo_matches")
        st.info("Pick an oligonucleotide in the oligonucleotide table to list its spectrum matches.")
        return
    _match_table(ctx, src, "oligo_matches", "Matches of the oligonucleotide",
                 filters={OLIGO_ID: "oligo_id"})


# --- the oligonucleotides of the picked RNA ----------------------------------

def _oligo_ids(src: Source) -> dict:
    """Every oligonucleotide sequence of the shown results, numbered in sorted order."""
    cached = getattr(src, "_oligo_ids", None)
    if cached is None:
        sequences = set()
        for _group, tables in _rna_groups(src):
            for name in ("OLI", "OSM"):
                table = tables.get(name, pd.DataFrame())
                if not table.empty and "sequence" in table:
                    sequences.update(table["sequence"].dropna().astype(str))
        cached = {sequence: index for index, sequence in enumerate(sorted(sequences))}
        src._oligo_ids = cached
    return cached


def oligo_rows(groups, rnas: pd.DataFrame, oligo_ids: dict) -> pd.DataFrame:
    """One row per oligonucleotide of each RNA of ``rnas`` (``rna_idx`` from it):
    its place in the RNA, the spectrum matches over ``groups`` and the best score."""
    counts, best = {}, {}
    for _group, tables in groups:
        osm = tables.get("OSM", pd.DataFrame())
        if osm.empty or "sequence" not in osm:
            continue
        scores = pd.to_numeric(osm.get("search_engine_score[1]", pd.Series(np.nan, index=osm.index)),
                               errors="coerce")
        for sequence, score in zip(osm["sequence"].astype(str), scores):
            counts[sequence] = counts.get(sequence, 0) + 1
            if np.isfinite(score) and (sequence not in best or score > best[sequence]):
                best[sequence] = float(score)
    rna_ids = dict(zip(rnas["accession"].astype(str), rnas["rna_idx"])) if not rnas.empty else {}
    rows, seen = [], set()
    for _group, tables in groups:
        oli = tables.get("OLI", pd.DataFrame())
        if oli.empty or "sequence" not in oli or "accession" not in oli:
            continue
        starts = pd.to_numeric(oli.get("start", pd.Series(np.nan, index=oli.index)), errors="coerce")
        ends = pd.to_numeric(oli.get("end", pd.Series(np.nan, index=oli.index)), errors="coerce")
        for accession, sequence, start, end in zip(oli["accession"].astype(str), oli["sequence"].astype(str),
                                                   starts, ends):
            if accession not in rna_ids or (accession, sequence, start, end) in seen:
                continue
            seen.add((accession, sequence, start, end))
            rows.append({"rna_idx": int(rna_ids[accession]), "accession": accession,
                         "oligo_id": int(oligo_ids.get(sequence, -1)), "sequence": sequence,
                         "start": None if pd.isna(start) else int(start),
                         "end": None if pd.isna(end) else int(end),
                         "matches": counts.get(sequence, 0), "score": best.get(sequence)})
    frame = pd.DataFrame(rows, columns=["rna_idx", "accession", "oligo_id", "sequence", "start", "end",
                                        "matches", "score"])
    frame.insert(0, "row_idx", np.arange(len(frame), dtype=np.int64))
    return frame


def oligo_table(ctx, src: Source) -> None:
    """The oligonucleotides of the picked RNA; a click lists their spectrum matches."""
    from openms_insight import Table

    from src.viz.columns import fixed
    from src.viz.plumbing import fingerprint

    rows = oligo_rows(_rna_groups(src), _rnas(src), _oligo_ids(src))
    if rows.empty:
        _heading("Oligonucleotides", "oligo_table")
        st.info("No identified oligonucleotide maps to an RNA of the database.")
        return
    if ctx.state.get_selection(RNA) is None:
        _heading("Oligonucleotides", "oligo_table")
        st.info("Pick an RNA in the RNA table to list its oligonucleotides.")
        return
    frame = pl.from_pandas(rows)
    columns = [
        {"field": "sequence", "title": "Oligonucleotide"},
        {"field": "start", "title": "Start", "sorter": "number"},
        {"field": "end", "title": "End", "sorter": "number"},
        {"field": "matches", "title": "Spectrum matches", "sorter": "number"},
        {"field": "score", "title": "Best score", "sorter": "number", **fixed(3)},
    ]
    Table(cache_id=f"oligos-{fingerprint(frame)}", cache_path=str(src.cache_dir), data=frame.lazy(),
          column_definitions=columns, index_field="row_idx", filters={RNA: "rna_idx"},
          interactivity={OLIGO_ID: "oligo_id"}, initial_sort=[{"column": "start", "dir": "asc"}],
          title="Oligonucleotides")(state_manager=ctx.state, height=300, key=key(ctx, src, "oligos"),
                                     help=teaching.component("oligo_table"))


# --- the fragment ions of the picked match -----------------------------------

def _idxml(src: Source) -> Path | None:
    folder = Path(src.wf.results_dir, "idxml_results")
    for name in ("search.idXML", "id_out"):
        if Path(folder, name).is_file():
            return Path(folder, name)
    return None


def picked_fragments(ctx, src: Source):
    """``(match, peaks, assignments, note)`` of the match picked in a match table:
    the match's row, its spectrum's peaks (pandas, mzML order), NASE's fragment
    assignments placed on them (empty when there are none) and why something
    is missing; None while no match is picked."""
    from src.naseweis_annotations import fragment_assignment_table, match_annotations
    from src.naseweis_viewers import VIEW_ERRORS, cached_annotations

    picked = ctx.state.get_selection(PSM)
    psms = _psms(src)
    if picked is None or psms is None or psms[0].is_empty():
        return None
    data, tables, _stamp = psms
    found = data.filter(pl.col("id_idx") == int(picked))
    if found.is_empty():
        return None
    match = found.row(0, named=True)
    empty = pd.DataFrame(columns=["mz", "intensity"])
    if tables is None:
        return match, empty, pd.DataFrame(), "The spectra of this acquisition are not available: " + "; ".join(src.notes)
    if int(match["scan_id"]) < 0:
        return match, empty, pd.DataFrame(), (
            "The matched spectrum is not in this acquisition's mzML; no approximate scan was substituted.")
    peaks = (tables["peaks"].filter(pl.col("index") == int(match["scan_id"]))
             .select("mz", "intensity").to_pandas())
    idxml = _idxml(src)
    if idxml is None:
        return match, peaks, pd.DataFrame(), "No fragment assignments were saved for this run."
    try:
        stat = idxml.stat()
        annotations = match_annotations(
            cached_annotations((str(idxml.resolve()), stat.st_mtime_ns, stat.st_size, stat.st_ctime_ns,
                                stat.st_ino, stat.st_dev)),
            match["native_id"], str(match["sequence"]), int(match["charge"]))
    except VIEW_ERRORS as error:
        return match, peaks, pd.DataFrame(), f"Fragment assignments are unavailable: {error}"
    if annotations.empty:
        return match, peaks, pd.DataFrame(), (
            "No saved fragment assignments match this spectrum, oligonucleotide and charge.")
    return match, peaks, fragment_assignment_table(peaks, annotations), None


def _fragment_note(ctx, src: Source, title: str, name: str):
    if _spectra(ctx, src, title, name) is None:
        return None
    picked = picked_fragments(ctx, src)
    if picked is None:
        _heading(title, name)
        st.info("Pick a spectrum match in a match table to see its fragment ions.")
    return picked


def match_spectrum(ctx, src: Source) -> None:
    """The picked match's spectrum, the peaks NASE assigned labelled by their fragment ion."""
    from openms_insight import LinePlot

    from src.viz.identification import fragment_spectrum
    from src.viz.plumbing import component, only_accepted

    picked = _fragment_note(ctx, src, "Fragment spectrum", "match_spectrum")
    if picked is None:
        return
    match, peaks, assignments, note = picked
    if peaks.empty:
        _heading("Fragment spectrum", "match_spectrum")
        st.info(note or "The matched spectrum holds no peaks.")
        return
    ions: dict[int, list[str]] = {}
    if not assignments.empty:
        for position, ion in zip(assignments["peak_index"], assignments["ion"]):
            if not pd.isna(position):
                ions.setdefault(int(position), []).append(str(ion))
    labels = {position: "; ".join(dict.fromkeys(names)) for position, names in ions.items()}
    psm = int(match["id_idx"])
    data = fragment_spectrum.rows(zip(peaks["mz"], peaks["intensity"]), psm, labels)
    title = f"Fragment spectrum: {match['sequence']} ({int(match['charge']):+d})"
    settings = fragment_spectrum.settings(spectrum=PSM, peak=FRAGMENT, title=title)
    plot = component(LinePlot, src.cache_dir, f"fragments_{_psms(src)[2]}_{_idxml_stamp(src)}_{psm}",
                     lambda: {"data": data.lazy()}, **only_accepted(LinePlot, **settings))
    plot(state_manager=ctx.state, height=360, key=key(ctx, src, "fragments"),
         help=teaching.component("match_spectrum"))
    if note:
        st.caption(note)


def _idxml_stamp(src: Source) -> str:
    idxml = _idxml(src)
    return str(idxml.stat().st_mtime_ns) if idxml is not None else "none"


def fragment_rows(assignments: pd.DataFrame, psm: int) -> pl.DataFrame:
    """The library fragment table's columns from NASE's assignments of match ``psm``."""
    from src.viz.identification.fragment_spectrum import ion_series

    observed = pd.to_numeric(assignments["observed_mz"], errors="coerce")
    expected = pd.to_numeric(assignments["mz"], errors="coerce")
    frame = pd.DataFrame({
        "fragment": np.arange(len(assignments), dtype=np.int64),
        "spectrum": psm,
        "peak": assignments["peak_index"].fillna(-1).astype("int64").to_numpy(),
        "ion": assignments["ion"].astype(str).to_numpy(),
        "series": [ion_series(ion) or "" for ion in assignments["ion"]],
        "mz": observed.to_numpy(),
        "expected_mz": expected.to_numpy(),
        "error_ppm": ((observed - expected) / expected * 1e6).to_numpy(),
        "intensity": pd.to_numeric(assignments["relative_intensity"], errors="coerce").to_numpy(),
    })
    return pl.from_pandas(frame)


def fragment_table(ctx, src: Source) -> None:
    """NASE's fragment assignments of the picked match; a click picks the peak in the fragment spectrum."""
    from openms_insight import Table

    from src.viz.identification import fragment_table as library
    from src.viz.plumbing import component

    picked = _fragment_note(ctx, src, "Fragment ions", "fragment_table")
    if picked is None:
        return
    match, _peaks, assignments, note = picked
    if assignments.empty:
        _heading("Fragment ions", "fragment_table")
        st.info(note or "NASE assigned no fragment ions in this spectrum.")
        return
    psm = int(match["id_idx"])
    data = fragment_rows(assignments, psm)
    table = component(Table, src.cache_dir, f"fragment_ions_v2_{_psms(src)[2]}_{_idxml_stamp(src)}_{psm}",
                      lambda: {"data": data.lazy()}, **library.settings(spectrum=PSM, peak=FRAGMENT, default_row_selects=False))
    table(state_manager=ctx.state, height=300, key=key(ctx, src, "fragment_ions"),
          help=teaching.component("fragment_table"))
    unplaced = int((data["peak"] < 0).sum())
    if unplaced:
        st.caption(f"{unplaced} assigned fragments match no single peak of the acquired spectrum "
                   "(for example after deisotoping or decharging); they are listed without one.")


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
        key=key(ctx, src, name), state_manager=ctx.state, height=320, help=teaching.component(name))


def score_distribution(ctx, src: Source) -> None:
    from src.viz.identification import score_distribution as viz

    loaded = _loaded(src)
    if loaded is None:
        return
    osm = loaded["tables"].get("OSM", pd.DataFrame())
    data = _histogram(osm.get("search_engine_score[1]", []))
    if data.empty:
        _heading("Search scores", "score_distribution")
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
        _heading("Precursor mass error", "precursor_error")
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
        _heading("Precursor mass error", "precursor_error")
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
    from src.naseweis_comparison import normalized
    from src.naseweis_experiment_ui import _conditions, comparison_choices

    design = loaded["design"]
    groups = {s: (str(design.loc[s, "condition"]) or "Unassigned") for s in design.index}
    names = {s: str(design.loc[s, "file"]) for s in design.index}
    saved = loaded["metadata"].get("normalization", "none")
    _, _, picked = comparison_choices(context["prefix"], _conditions(design), saved)
    abundances = loaded["abundances"] if picked == saved else normalized(loaded["raw"], picked)
    return abundances, groups, names


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
          title="Abundances")(state_manager=ctx.state, height=376, key=key(ctx, src, "abundances"), help=teaching.component("abundance_table"))


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
        _heading("Abundance profile", "abundance_profile")
        st.info("Pick an oligonucleotide in the abundance table.")
        return
    show_profile(long, "oligo_idx", OLIGO, "nase-abundance-profile", cache_path=src.cache_dir,
                 key=key(ctx, src, "profile"), state_manager=ctx.state, height=420,
                 table_args={"hover_column": "label"}, title="Abundance profile",
                 help=teaching.component("abundance_profile"),
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
                                title="Abundance z-scores", y_label="RNA oligonucleotides",
                                help=teaching.component("zscore_heatmap"))
    if shown is None:
        _heading("Abundance z-scores", "zscore_heatmap")
        st.info("No oligonucleotide differs between the samples.")
