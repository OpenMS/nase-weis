"""Workspace-aware viewers using the template's existing pyOpenMS-viz backend.

Only display overviews are binned. Spectrum peaks and XICs retain their native
values, and result viewers resolve acquisitions from the selected run snapshot.
"""

import hashlib
import json
import math
import re
from pathlib import Path

import pandas as pd
import streamlit as st

from src.common.common import follow_app_theme
from src.naseweis_annotations import (
    annotate_peaks,
    fragment_assignment_table,
    match_annotations,
    native_spectrum_id,
    read_annotations,
)
from src.naseweis_msdata import (
    file_revision,
    ms1_overview,
    scan_index,
    spectrum_peaks,
    xic,
)
from src.naseweis_plots import (
    chromatogram_figure,
    feature_map_figure,
    ion_series,
    peakmap_figure,
    spectrum_figure,
)

VIEW_ERRORS = (OSError, ValueError, RuntimeError)


def viewer_key(value) -> str:
    """Keep file/run state independent without exposing paths as widget keys."""
    return hashlib.sha256(str(value).encode()).hexdigest()[:16]


def read_run_snapshot(path: Path) -> dict:
    """Read this run's scientific provenance, never the live configuration."""
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ValueError("This run's parameter snapshot is unavailable.") from error
    if not isinstance(value, dict):
        raise ValueError("This run's parameter snapshot is invalid.")  # noqa: TRY004
    return value


def acquisition_revision(selected, workflow_dir: Path, online: bool = False):
    """Respect workspace boundaries online and registered external files locally."""
    if not isinstance(selected, str) or not selected:
        raise ValueError("This run has no saved acquisition path.")
    path = Path(selected).resolve()
    if online and not path.is_relative_to(workflow_dir.parent.resolve()):
        raise ValueError("The saved acquisition is outside this online workspace.")
    return file_revision(path)


@st.cache_data(show_spinner=False, max_entries=8)
def cached_fingerprint(revision):
    from src.naseweis_provenance import acquisition_fingerprint

    return acquisition_fingerprint(revision[0])


def saved_acquisition_revision(
    params: dict, workflow_dir: Path, online=False, *, verify_content=True
):
    """Do not show replacement bytes as evidence for a previously searched file."""
    revision = acquisition_revision(params.get("mzML-files"), workflow_dir, online)
    if "_acquisition_fingerprint" not in params:
        from src.naseweis_provenance import (
            acquisition_fingerprint,
            requires_content_verification,
        )

        if verify_content and requires_content_verification():
            return (*revision, acquisition_fingerprint(revision[0])["sha256"])
        return revision
    recorded = params["_acquisition_fingerprint"]
    if not isinstance(recorded, dict):
        raise ValueError("The saved acquisition fingerprint is invalid.")  # noqa: TRY004
    checksum, size = recorded.get("sha256"), recorded.get("size")
    if (
        not isinstance(checksum, str)
        or re.fullmatch(r"[0-9a-fA-F]{64}", checksum) is None
        or isinstance(size, bool)
        or not isinstance(size, int)
        or size < 0
    ):
        raise ValueError("The saved acquisition fingerprint is invalid.")
    actual = None
    if verify_content and size == revision[2]:
        from src.naseweis_provenance import (
            acquisition_fingerprint,
            requires_content_verification,
        )

        actual = (
            acquisition_fingerprint(revision[0])
            if requires_content_verification()
            else cached_fingerprint(revision)
        )
    if size != revision[2] or (
        actual is not None and actual["sha256"] != checksum.lower()
    ):
        raise ValueError(
            "The acquisition has changed since this run was searched. Restore the original mzML; a replacement file cannot be used as this run's evidence."
        )
    return (*revision, actual["sha256"]) if actual is not None else revision


def _legacy_identity_notice(params):
    if "_acquisition_fingerprint" not in params:
        st.warning(
            "Acquisition identity was not recorded for this older run. Verify that this is its original mzML before interpreting the evidence."
        )


def input_acquisitions(workflow_dir: Path, params: dict, online=False) -> list[str]:
    """Find uploaded/registered inputs without changing the search selection."""
    folders = [
        workflow_dir / "input-files" / "mzML-files",
        workflow_dir.parent / "mzML-files",  # Earlier template workspaces.
    ]
    candidates = set()
    selected = params.get("mzML-files")
    if isinstance(selected, str) and selected:
        candidates.add(selected)
    elif isinstance(selected, list):
        candidates.update(value for value in selected if isinstance(value, str) and value)
    for folder in folders:
        if not folder.is_dir():
            continue
        candidates.update(
            str(path)
            for path in folder.iterdir()
            if path.suffix.lower() in {".mzml", ".raw", ".d"}
        )
        registered = folder / "external_files.txt"
        if registered.is_file():
            candidates.update(
                line.strip()
                for line in registered.read_text().splitlines()
                if line.strip()
            )
    result = set()
    for candidate in candidates:
        path = Path(candidate).resolve()
        if path.suffix.lower() not in {".mzml", ".raw", ".d"} or not path.exists():
            continue
        if online and not path.is_relative_to(workflow_dir.parent.resolve()):
            continue
        result.add(str(path))
    return sorted(result, key=lambda value: (Path(value).name.casefold(), value))


@st.cache_data(show_spinner=False, max_entries=4)
def cached_scan_index(revision):
    return scan_index(revision[0])


@st.cache_data(show_spinner=False, max_entries=4)
def cached_overview(revision):
    return ms1_overview(revision[0])


@st.cache_data(show_spinner=False, max_entries=32)
def cached_spectrum(revision, index):
    return spectrum_peaks(revision[0], int(index))


@st.cache_data(show_spinner=False, max_entries=16)
def cached_xic(revision, mz, tolerance, unit, rt_range=None):
    return xic(revision[0], mz, tolerance, unit, rt_range)


@st.cache_data(show_spinner=False, max_entries=8)
def cached_annotations(revision):
    return read_annotations(revision[0])


def show_viewer_plot(figure, key: str, filename: str) -> None:
    """Keep zoom/pan/export controls and provide a self-contained offline plot."""
    follow_app_theme(figure)
    st.plotly_chart(
        figure,
        key=key,
        use_container_width=True,
        config={
            "displaylogo": False,
            "scrollZoom": True,
            "toImageButtonOptions": {"filename": filename, "format": "svg"},
        },
    )
    st.download_button(
        "Download interactive plot (HTML)",
        figure.to_html(full_html=True, include_plotlyjs=True),
        f"{filename}.html",
        "text/html",
        key=f"{key}-html",
    )


def render_raw_viewer(workflow_dir: Path, params: dict, online=False) -> None:
    files = input_acquisitions(workflow_dir, params, online)
    if not files:
        st.info(
            "Upload an acquisition or load the NME1 example to explore its spectra and chromatograms."
        )
        st.page_link(
            "content/topp_workflow_file_upload.py",
            label="Add data",
            icon=":material/upload_file:",
        )
        return
    selected = params.get("mzML-files")
    if isinstance(selected, list):
        selected = selected[0] if selected else None
    path = st.selectbox(
        "Acquisition",
        files,
        index=files.index(selected) if selected in files else 0,
        format_func=lambda value: f"{Path(value).name} · {Path(value).parent.name}",
        key=f"naseweis-raw-file-{viewer_key(workflow_dir)}",
    )
    try:
        revision = acquisition_revision(path, workflow_dir, online)
    except VIEW_ERRORS as error:
        st.info(str(error))
        return
    key = f"naseweis-raw-{viewer_key(revision)}"
    if st.button("Load acquisition", key=f"{key}-load", type="primary"):
        st.session_state[f"{key}-loaded"] = True
    if not st.session_state.get(f"{key}-loaded", False):
        st.caption(
            "Load the selected acquisition to browse its MS1 peak map, individual spectra and chromatograms."
        )
        return
    try:
        with st.spinner("Reading acquisition…"):
            from src.naseweis_provenance import (
                acquisition_fingerprint,
                requires_content_verification,
            )

            if requires_content_verification():
                revision = (*revision, acquisition_fingerprint(revision[0])["sha256"])
            scans = cached_scan_index(revision)
        if scans.empty:
            st.info("This acquisition contains no spectra.")
            return
        peak_tab, spectrum_tab, chromatogram_tab = st.tabs(
            ["MS1 peak map", "Spectra", "Chromatograms"]
        )
        try:
            with st.spinner("Preparing the MS1 overview…"):
                points, traces = cached_overview(revision)
            ms1_error = None
        except VIEW_ERRORS as error:
            points, traces, ms1_error = pd.DataFrame(), pd.DataFrame(), str(error)
        with peak_tab:
            if ms1_error:
                st.info(ms1_error)
            elif not points.empty:
                marginals = st.checkbox("Show marginal signal", key=f"{key}-marginals")
                show_viewer_plot(
                    peakmap_figure(points, Path(path).name, marginals=marginals),
                    f"{key}-peakmap",
                    "naseweis-ms1-peakmap",
                )
                st.caption(
                    "The overview sums intensities within display bins. Individual spectra and XICs use the original peak data. Zoom, pan and hover using the plot controls."
                )
            else:
                st.info("No MS1 peaks to display.")
        with spectrum_tab:
            levels = ["All"] + [
                str(value) for value in sorted(scans["ms_level"].unique())
            ]
            level = st.selectbox("MS level", levels, key=f"{key}-level")
            query = st.text_input("Find spectrum by native ID", key=f"{key}-query")
            visible = (
                scans if level == "All" else scans.loc[scans["ms_level"] == int(level)]
            )
            if query:
                visible = visible.loc[
                    visible["native_id"].str.contains(query, regex=False, case=False)
                ]
            if visible.empty:
                st.info("No spectra match these filters.")
            else:
                if len(visible) > 1000:
                    st.caption(
                        "Showing the first 1,000 matches. Refine the native-ID filter to find another scan."
                    )
                choices = visible.head(1000).set_index("index", drop=False)
                index = st.selectbox(
                    "Spectrum",
                    list(choices.index),
                    format_func=lambda value: (
                        f"{choices.loc[value, 'native_id']} · MS{choices.loc[value, 'ms_level']} · {choices.loc[value, 'rt']:.2f} s"
                    ),
                    key=f"{key}-spectrum-{viewer_key((level, query))}",
                )
                peaks = cached_spectrum(revision, index)
                if peaks.empty:
                    st.info("This spectrum has no peaks.")
                else:
                    show_viewer_plot(
                        spectrum_figure(peaks, str(choices.loc[index, "native_id"])),
                        f"{key}-scan-plot",
                        "naseweis-spectrum",
                    )
                with st.expander("Scan metadata"):
                    st.dataframe(choices, hide_index=True, use_container_width=True)
        with chromatogram_tab:
            if ms1_error:
                st.info(ms1_error)
            else:
                selected_traces = st.multiselect(
                    "Signal traces",
                    ["TIC", "BPC"],
                    default=["TIC", "BPC"],
                    key=f"{key}-traces",
                )
                signals = traces.loc[traces["trace"].isin(selected_traces)]
                if st.checkbox("Add an extracted ion chromatogram", key=f"{key}-xic"):
                    mz, tolerance, unit = _xic_controls(key, 500.0, 10.0, "ppm")
                    signals = pd.concat(
                        [signals, cached_xic(revision, mz, tolerance, unit)],
                        ignore_index=True,
                    )
                if not signals.empty:
                    show_viewer_plot(
                        chromatogram_figure(signals, Path(path).name),
                        f"{key}-chromatograms",
                        "naseweis-chromatograms",
                    )
                else:
                    st.info("Choose a signal trace to display.")
    except VIEW_ERRORS as error:
        st.error(f"Could not display this acquisition: {error}")


def _xic_controls(key, mz, tolerance, unit):
    columns = st.columns(3)
    mz = columns[0].number_input(
        "XIC m/z",
        min_value=0.00001,
        value=float(mz),
        format="%.5f",
        key=f"{key}-xic-mz",
    )
    selected_unit = columns[2].selectbox(
        "Tolerance unit",
        ["ppm", "Da"],
        index=0 if unit == "ppm" else 1,
        key=f"{key}-xic-unit",
    )
    default = (
        tolerance
        if selected_unit == unit
        else (tolerance * mz / 1e6 if unit == "ppm" else tolerance / mz * 1e6)
    )
    tolerance = columns[1].number_input(
        f"Tolerance (± {selected_unit})",
        min_value=0.000001,
        value=float(default),
        format="%.6f",
        key=f"{key}-xic-tolerance-{selected_unit}",
    )
    return mz, tolerance, selected_unit


def _fragment_assignment_controls(peaks, annotations, key):
    """Link visible native assignments to exact peaks in the acquired scan."""
    series = annotations["ion"].map(lambda ion: ion_series(ion) or "Other")
    options = list(dict.fromkeys(series))
    selected_series = st.multiselect(
        "Ion series", options, default=options, key=f"{key}-series",
        help="Filter the displayed native assignments and their peak labels. Search results are unchanged.",
    )
    shown = annotations.loc[series.isin(selected_series)]
    labelled, _ = annotate_peaks(peaks, shown)
    table = fragment_assignment_table(peaks, shown)
    st.subheader("Assigned fragment ions")
    st.caption(
        "Select rows using the checkboxes to highlight their peaks in gold. "
        "Intensities below come from the original acquired scan."
    )
    if table.empty:
        st.info("No assigned ions in the selected series. Choose an ion series to show its assignments.")
        return labelled, []
    located = table["match_status"].eq("matched")
    st.caption(
        f"{len(table)} assignments · {table.loc[located, 'peak_index'].nunique()} located raw peaks · "
        f"{int((~located).sum())} unplaced assignments"
    )
    # Row positions belong to the displayed table. Changing its filter or
    # native assignments must not apply a previous selection to different ions.
    selection_key = f"{key}-ions-{viewer_key(table.to_csv(index=False))}"
    event = st.dataframe(
        table,
        column_order=["ion", "fragment_charge", "mz", "observed_mz", "observed_intensity", "relative_intensity", "match_status"],
        column_config={
            "ion": st.column_config.TextColumn("Ion"),
            "fragment_charge": st.column_config.NumberColumn("Charge", format="%d"),
            "mz": st.column_config.NumberColumn("Assigned m/z", format="%.6f"),
            "observed_mz": st.column_config.NumberColumn("Raw m/z", format="%.6f"),
            "observed_intensity": st.column_config.NumberColumn("Raw intensity", format="%.2f"),
            "relative_intensity": st.column_config.NumberColumn("Relative intensity (%)", format="%.2f"),
            "match_status": st.column_config.TextColumn("Raw peak match"),
        },
        hide_index=True, width="stretch", height=300,
        on_select="rerun", selection_mode="multi-row", key=selection_key,
    )
    chosen = table.iloc[event["selection"]["rows"]]
    placed = chosen.loc[chosen["match_status"].eq("matched")]
    if not chosen.empty:
        st.caption("Selected ions: " + ", ".join(chosen["ion"].astype(str)))
    if len(placed) < len(chosen):
        st.info(
            "Some selected assignments have no unique peak in the original scan and cannot be highlighted. "
            "NASE may have scored a deisotoped or decharged spectrum."
        )
    return labelled, placed["observed_mz"].drop_duplicates().tolist()


def render_identification_viewer(
    results_dir: Path,
    snapshot: Path,
    matches: pd.DataFrame,
    workflow_dir: Path,
    online=False,
):
    st.subheader("RNA identification viewer")
    if matches.empty or "spectra_ref" not in matches:
        st.info("No spectrum references are available for this run.")
        return
    try:
        params = read_run_snapshot(snapshot)
        revision = saved_acquisition_revision(
            params, workflow_dir, online, verify_content=False
        )
    except VIEW_ERRORS as error:
        st.info(
            f"Spectrum viewing is unavailable: {error} Restore this run's original mzML acquisition to view its spectra."
        )
        return
    _legacy_identity_notice(params)
    key = f"naseweis-match-{viewer_key((results_dir.resolve(), revision))}"
    if st.button("Open spectrum viewer", key=f"{key}-load", type="primary"):
        st.session_state[f"{key}-loaded"] = True
    if not st.session_state.get(f"{key}-loaded", False):
        st.caption(
            "Inspect the original scan and the RNA fragment assignments saved by NASE. Scan references are matched by exact native ID."
        )
        return
    try:
        # Hash only after the viewer is opened; the Matches tab stays cheap for
        # large acquisitions while every displayed scan still verifies identity.
        revision = saved_acquisition_revision(params, workflow_dir, online)
        query = st.text_input(
            "Find identification", key=f"{key}-find",
            placeholder="RNA sequence or native spectrum ID",
        ).strip().casefold()
        choices = list(range(len(matches)))
        if query:
            choices = [
                position for position in choices
                if query in str(matches.iloc[position].get("sequence", "")).casefold()
                or query in str(matches.iloc[position]["spectra_ref"]).casefold()
            ]
        if not choices:
            st.info("No identifications match this filter. Clear or change the search to choose another identification.")
            return
        index = st.selectbox(
            "Spectrum match",
            choices,
            format_func=lambda value: (
                f"{matches.iloc[value].get('sequence', '')} · z={matches.iloc[value].get('charge', '')} · {matches.iloc[value].get('retention_time', '')} s · {native_spectrum_id(str(matches.iloc[value]['spectra_ref']))}"
            ),
            key=f"{key}-match",
        )
        match = matches.iloc[index]
        native_id = native_spectrum_id(str(match["spectra_ref"]))
        scans = cached_scan_index(revision)
        found = scans.loc[scans["native_id"] == native_id]
        if len(found) != 1 or int(found.iloc[0]["ms_level"]) != 2:
            st.warning(
                "The saved spectrum reference does not identify exactly one scan in this acquisition. No approximate RT or scan-number match was substituted."
            )
            return
        peaks = cached_spectrum(revision, int(found.iloc[0]["index"]))
        st.caption(f"{Path(revision[0]).name} · {native_id}")
        annotations, unmatched = pd.DataFrame(), pd.DataFrame()
        idxml = results_dir / "idxml_results" / "search.idXML"
        if not idxml.is_file():
            idxml = results_dir / "idxml_results" / "id_out"
        if idxml.is_file():
            try:
                stat = idxml.stat()
                annotations = match_annotations(
                    cached_annotations(
                        (
                            str(idxml.resolve()),
                            stat.st_mtime_ns,
                            stat.st_size,
                            stat.st_ctime_ns,
                            stat.st_ino,
                            stat.st_dev,
                        )
                    ),
                    native_id,
                    str(match.get("sequence", "")),
                    int(match.get("charge", 0)),
                )
                peaks, unmatched = annotate_peaks(peaks, annotations)
                if annotations.empty:
                    st.info("No saved fragment assignments match this exact spectrum, RNA sequence and signed charge.")
            except VIEW_ERRORS as error:
                annotations = pd.DataFrame()
                st.warning(f"Native fragment annotations are unavailable: {error}")
        else:
            st.caption(
                "No idXML annotations were saved for this run; the acquired spectrum remains available."
            )
        # Reserve the plot's position above its linked table. The table must
        # execute first so this rerun can render the current row selection.
        plot_area = st.container()
        highlighted_mz = []
        if not annotations.empty:
            match_key = f"{key}-{viewer_key((index, native_id, str(match.get('sequence', '')), str(match.get('charge', ''))))}"
            peaks, highlighted_mz = _fragment_assignment_controls(peaks, annotations, match_key)
        if peaks.empty:
            st.info("The acquired scan contains no peaks.")
        else:
            with plot_area:
                show_viewer_plot(
                    spectrum_figure(peaks, str(match.get("sequence", native_id)), highlighted_mz=highlighted_mz),
                    f"{key}-spectrum",
                    "naseweis-matched-spectrum",
                )
        if not annotations.empty:
            if not unmatched.empty:
                st.info(
                    "Some scored fragments are absent or ambiguous in the acquired scan, for example after deisotoping or decharging. Their native assignments are retained below."
                )
            with st.expander("All native assignment metadata"):
                st.dataframe(annotations, hide_index=True, use_container_width=True)
            st.download_button(
                "Download fragment assignments (TSV)",
                annotations.to_csv(sep="\t", index=False),
                "naseweis-fragments.tsv",
                "text/tab-separated-values",
                key=f"{key}-annotations",
            )
            st.caption(
                "Labels come from NASE's saved assignments. Acquired peak heights are preserved; processed fragment masses are labelled only when they uniquely match an acquired peak."
            )
    except VIEW_ERRORS as error:
        st.error(f"Could not display this spectrum: {error}")


def _feature_range(feature):
    """Use native detected boundaries when available, then pad for context."""
    for lower, upper in (("leftWidth", "rightWidth"), ("rt_start", "rt_end")):
        try:
            left, right = float(feature[lower]), float(feature[upper])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(left) and math.isfinite(right) and left < right:
            padding = max((right - left) / 2, 5.0)
            return (left - padding, right + padding), (left, right)
    rt = float(feature["rt"])
    return (rt - 120.0, rt + 120.0), None


def render_lfq_viewer(
    results_dir: Path,
    snapshot: Path,
    features: pd.DataFrame,
    workflow_dir: Path,
    online=False,
):
    if features.empty or not {"rt", "mz", "intensity"} <= set(features):
        return
    key = f"naseweis-feature-{viewer_key(results_dir.resolve())}"
    if st.button("Open LFQ feature viewer", key=f"{key}-load"):
        st.session_state[f"{key}-loaded"] = True
    if not st.session_state.get(f"{key}-loaded", False):
        return
    try:
        show_viewer_plot(
            feature_map_figure(features, "RNA MS1 features"),
            f"{key}-map",
            "naseweis-lfq-feature-map",
        )
        index = st.selectbox(
            "LFQ feature",
            list(range(len(features))),
            format_func=lambda value: (
                f"{features.iloc[value].get('label', 'Feature')} · z={features.iloc[value].get('charge', '')} · {features.iloc[value]['rt']:.2f} s"
            ),
            key=f"{key}-feature",
        )
        feature = features.iloc[index]
        try:
            params = read_run_snapshot(snapshot)
            revision = saved_acquisition_revision(params, workflow_dir, online)
        except VIEW_ERRORS as error:
            st.info(
                f"The feature map is available, but its acquisition cannot be opened: {error}"
            )
            return
        _legacy_identity_notice(params)
        settings = params.get("FeatureFinderMetaboIdent", {})
        if not isinstance(settings, dict):
            raise ValueError("This run's saved LFQ extraction settings are invalid.")  # noqa: TRY004
        value = settings.get("extract:mz_window", 50.0)
        try:
            if isinstance(value, bool):
                raise TypeError("Boolean extraction window")
            native_window = float(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(
                "This run has no valid saved extraction window."
            ) from error
        if not math.isfinite(native_window) or native_window <= 0:
            raise ValueError("This run has no valid saved extraction window.")
        unit = "ppm" if native_window >= 1 else "Da"
        control_key = f"{key}-{viewer_key((revision, index))}"
        mz, tolerance, unit = _xic_controls(
            control_key, float(feature["mz"]), native_window / 2, unit
        )
        limited = st.checkbox(
            "Focus on the feature's elution window", value=True, key=f"{key}-rt-window"
        )
        rt_range, boundaries = _feature_range(feature)
        with st.spinner("Extracting the MS1 chromatogram…"):
            trace = cached_xic(
                revision, mz, tolerance, unit, rt_range if limited else None
            )
        figure = chromatogram_figure(trace, str(feature.get("label", "RNA feature")))
        if boundaries:
            figure.add_vrect(
                x0=boundaries[0],
                x1=boundaries[1],
                fillcolor="#ffcb4f",
                opacity=0.18,
                line_width=0,
            )
        figure.add_vline(x=float(feature["rt"]), line_dash="dot", line_color="#3157e9")
        show_viewer_plot(figure, f"{key}-xic-plot", "naseweis-lfq-xic")
        st.caption(
            "The XIC sums raw MS1 signal within the displayed tolerance; the shaded interval marks native feature boundaries. LFQ areas can include multiple isotope traces and are not identical to a single-m/z XIC. Adjusting this view does not change saved quantification."
        )
    except VIEW_ERRORS as error:
        st.warning(f"The LFQ viewer could not display this data: {error}")
