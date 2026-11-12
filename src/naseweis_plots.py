"""Pure, offline figures for NASE spectra, raw acquisitions and RNA features.

The mass-spectrometry plots use the application's pinned pyopenms-viz backend.
Plotly annotations work around its incomplete feature and RNA label support.
"""

import re
from html import escape

import numpy as np
import pandas as pd
import plotly.graph_objects as go

ION_COLORS = {
    "a-B": "#0072B2",
    "a": "#56B4E9",
    "b": "#009E73",
    "c": "#E69F00",
    "d": "#CC79A7",
    "w": "#D55E00",
    "x": "#F0E442",
    "y": "#CC0000",
    "z": "#999999",
}
UNKNOWN_ION_COLOR = "#6B7280"
SELECTED_ION_COLOR = "#C58B00"
TRACE_COLORS = {"TIC": "#0072B2", "BPC": "#009E73", "XIC": "#D55E00"}


def _text(value) -> str:
    """Treat annotations as text, including dollar signs used by MathJax."""
    return escape(str(value), quote=True).replace("$", "&#36;")


def _frame(data, columns: tuple[str, ...], context: str) -> pd.DataFrame:
    if not isinstance(data, pd.DataFrame):
        # All invalid viewer inputs share ValueError for the UI's error path.
        raise ValueError(f"{context} must be a table.")  # noqa: TRY004
    missing = set(columns) - set(data.columns)
    if missing:
        raise ValueError(f"{context} is missing columns: {', '.join(sorted(missing))}.")
    if data.empty:
        raise ValueError(f"No {context.lower()} are available to plot.")
    # Backend hover data includes the index; use numeric row positions rather
    # than serializing arbitrary labels supplied as a DataFrame index.
    frame = data.copy(deep=True).reset_index(drop=True)
    for column in columns:
        try:
            values = pd.to_numeric(frame[column], errors="raise")
            finite = np.isfinite(values.to_numpy(dtype=float, na_value=np.nan))
        except (ValueError, TypeError, OverflowError):
            raise ValueError(f"{context} column {column} must contain finite numbers.") from None
        if not finite.all():
            raise ValueError(f"{context} column {column} must contain finite numbers.")
        frame[column] = values
    if "intensity" in frame and (frame["intensity"] < 0).any():
        raise ValueError(f"{context} intensities must not be negative.")
    if "mz" in frame and (frame["mz"] <= 0).any():
        raise ValueError(f"{context} m/z values must be positive.")
    return frame


def ion_series(ion: str) -> str:
    """Identify a native RNA ion series without changing the assignment text."""
    # OpenMS labels a-B ions both as a-B1 and a1-B, depending on their source.
    if re.match(r"^a(?:\d+)?-B", ion):
        return "a-B"
    match = re.match(r"^([abcdwxyz])(?:\d|$|\^|[+-])", ion)
    return match.group(1) if match else ""


def _spectrum(data: pd.DataFrame, context: str) -> pd.DataFrame:
    frame = _frame(data, ("mz", "intensity"), context)
    if frame["intensity"].max() == 0:
        raise ValueError(f"{context} has no positive intensities to normalize.")
    frame["ion"] = frame.get("ion", pd.Series("", index=frame.index)).fillna("").astype(str)
    # The backend automatically serializes native_id when present. Keep only
    # columns owned by this figure so unrelated metadata cannot enter the HTML.
    frame = frame[["mz", "intensity", "ion"]].copy()
    # The backend overwrites a column literally named 'color'; use our own.
    frame["_naseweis_color"] = frame["ion"].map(
        lambda ion: ION_COLORS.get(ion_series(ion), UNKNOWN_ION_COLOR)
    )
    frame["_naseweis_ion"] = frame["ion"].map(_text)
    return frame.sort_values("mz", kind="stable").reset_index(drop=True)


def spectrum_figure(
    peaks: pd.DataFrame, title: str, reference: pd.DataFrame | None = None,
    *, highlighted_mz=None,
) -> go.Figure:
    """Plot unbinned peaks, optionally mirrored, with explicit RNA ion colors.

    Each spectrum is normalized to its own maximum; only the display changes.
    Supplied ion annotations retain their exact peak m/z coordinates.
    Highlights identify exact, unambiguous acquired peak masses; no peak is
    snapped to a nearby mass or relabelled as an assigned ion.
    """
    spectrum = _spectrum(peaks, "Spectrum peaks")
    selected = _highlighted_peaks(spectrum, highlighted_mz)
    mirror = None if reference is None else _spectrum(reference, "Reference peaks")
    figure = spectrum.plot(
        backend="ms_plotly", kind="spectrum", x="mz", y="intensity",
        show_plot=False, title=_text(title), grid=False, bin_peaks=False,
        aggregate_duplicates=False, relative_intensity=True,
        reference_spectrum=mirror, mirror_spectrum=mirror is not None,
        ion_annotation="_naseweis_ion", peak_color="_naseweis_color",
        annotation_color="_naseweis_color", annotate_top_n_peaks=None,
        annotate_mz=False,
    )
    color_names = {color: series for series, color in ION_COLORS.items()}
    for trace in figure.data:
        mirrored = any(value < 0 for value in trace.y if value is not None)
        series = color_names.get(trace.line.color, "Unannotated")
        trace.name = f"{series} (reference)" if mirrored else series
        trace.hovertemplate = "m/z: %{x}<br>Relative intensity: %{y:.2f}%<extra>%{fullData.name}</extra>"
    # pyopenms-viz 1.0.0 returns black annotations even with annotation_color.
    # Place the supplied labels ourselves, including labels on mirror peaks.
    for data, sign in ((spectrum, 1), (mirror, -1)):
        if data is None:
            continue
        maximum = data["intensity"].max()
        for row in data.loc[data["ion"] != ""].itertuples(index=False):
            figure.add_annotation(
                x=row.mz, y=sign * row.intensity / maximum * 100,
                text=_text(row.ion), showarrow=False, xanchor="center",
                yanchor="bottom" if sign > 0 else "top",
                font={"color": ION_COLORS.get(ion_series(row.ion), UNKNOWN_ION_COLOR), "size": 11},
            )
    if not selected.empty:
        height = selected["intensity"] / spectrum["intensity"].max() * 100
        # An open marker preserves the underlying series color and exact stick
        # height. A narrow translucent line makes weak selected ions visible.
        for mz, intensity in zip(selected["mz"], height):
            figure.add_shape(
                type="line", x0=mz, x1=mz, y0=0, y1=intensity,
                line={"color": SELECTED_ION_COLOR, "width": 6},
                opacity=0.35, layer="below",
            )
        figure.add_trace(go.Scatter(
            x=selected["mz"], y=height, mode="markers",
            name="Selected ion", text=selected["_naseweis_ion"],
            marker={"symbol": "circle-open", "size": 14,
                    "color": SELECTED_ION_COLOR, "line": {"width": 3}},
            hovertemplate="%{text}<br>m/z: %{x}<br>Relative intensity: %{y:.2f}%<extra>Selected ion</extra>",
        ))
    figure.update_layout(title=_text(title), xaxis_title="m/z", yaxis_title="Relative intensity (%)")
    return figure


def _highlighted_peaks(spectrum: pd.DataFrame, highlighted_mz) -> pd.DataFrame:
    """Resolve display selections without mass tolerance or changing signal."""
    if highlighted_mz is None:
        return spectrum.iloc[:0]
    try:
        values = list(highlighted_mz)
        if isinstance(highlighted_mz, (str, bytes)) or any(
            isinstance(value, (bool, np.bool_)) for value in values
        ):
            raise ValueError
        masses = np.asarray(values, dtype=float)
        if masses.ndim != 1 or not np.isfinite(masses).all() or (masses <= 0).any():
            raise ValueError
    except (TypeError, ValueError, OverflowError):
        raise ValueError("Highlighted m/z values must be positive finite numbers.") from None
    for mass in masses:
        if (spectrum["mz"] == mass).sum() != 1:
            raise ValueError("Each highlighted m/z must identify exactly one acquired peak.")
    return spectrum.loc[spectrum["mz"].isin(masses)]


def chromatogram_figure(points: pd.DataFrame, title: str) -> go.Figure:
    """Overlay TIC, BPC, XIC or isotope traces without changing intensities."""
    frame = _frame(points, ("rt", "intensity"), "Chromatogram points")
    if "trace" not in frame or frame["trace"].isna().any():
        raise ValueError("Chromatogram points must identify each trace.")
    frame["trace"] = frame["trace"].map(_text)
    frame = frame.sort_values(["trace", "rt"], kind="stable")
    figure = frame.plot(
        backend="ms_plotly", kind="chromatogram", x="rt", y="intensity",
        by="trace", show_plot=False, title=_text(title), grid=False,
        aggregate_duplicates=False, relative_intensity=False,
    )
    # The backend adds an empty drawing-helper trace that has no data to view.
    figure.data = tuple(trace for trace in figure.data if len(trace.x))
    for trace in figure.data:
        if trace.name in TRACE_COLORS:
            trace.line.color = TRACE_COLORS[trace.name]
        trace.hovertemplate = "Retention time: %{x} s<br>Intensity: %{y}<extra>%{fullData.name}</extra>"
    figure.update_layout(xaxis_title="Retention time (s)", yaxis_title="Intensity")
    return figure


def peakmap_figure(
    points: pd.DataFrame, title: str, bounds: pd.DataFrame | None = None,
    marginals: bool = False,
) -> go.Figure:
    """Plot prepared MS1 points with optional actual feature boundaries.

    The caller owns any overview binning. Boundaries are not inferred from
    centroids, extraction windows or instrument tolerances.
    """
    frame = _frame(points, ("rt", "mz", "intensity"), "Peak-map points")
    boxes = None
    if bounds is not None and not isinstance(bounds, pd.DataFrame):
        raise ValueError("Feature boundaries must be a table.")
    if bounds is not None and not bounds.empty:
        boxes = _frame(bounds, ("rt_start", "rt_end", "mz_start", "mz_end"), "Feature boundaries")
        if (boxes["rt_start"] > boxes["rt_end"]).any() or (boxes["mz_start"] > boxes["mz_end"]).any():
            raise ValueError("Feature boundaries must have their lower limits before their upper limits.")
        if (boxes[["mz_start", "mz_end"]] <= 0).any().any():
            raise ValueError("Feature m/z boundaries must be positive.")
    figure = frame.plot(
        backend="ms_plotly", kind="peakmap", x="rt", y="mz", z="intensity",
        show_plot=False, title=_text(title), grid=False, bin_peaks=False,
        aggregate_duplicates=False, z_log_scale=False, relative_intensity=False,
        add_marginals=marginals, marker_size=4,
        xlabel="Retention time (s)", ylabel="m/z", zlabel="Intensity",
    )
    # The first trace is the main peak map, including in the marginal layout.
    figure.data[0].marker.update(
        colorscale="Viridis", showscale=True,
        cmin=0, cmax=float(frame["intensity"].max()),
        colorbar={"title": "Intensity"},
    )
    xaxis, yaxis = figure.data[0].xaxis or "x", figure.data[0].yaxis or "y"
    if boxes is not None:
        for number, (_, row) in enumerate(boxes.iterrows(), 1):
            raw_label = row.get("label", f"Feature {number}")
            label = _text("" if pd.isna(raw_label) else raw_label)
            figure.add_shape(
                type="rect", x0=row["rt_start"], x1=row["rt_end"],
                y0=row["mz_start"], y1=row["mz_end"], xref=xaxis, yref=yaxis,
                line={"color": ION_COLORS["a-B"], "width": 2},
                fillcolor="rgba(0,0,0,0)",
            )
            figure.add_trace(go.Scatter(
                x=[(row["rt_start"] + row["rt_end"]) / 2],
                y=[(row["mz_start"] + row["mz_end"]) / 2],
                xaxis=xaxis, yaxis=yaxis, mode="markers",
                marker={"size": 12, "opacity": 0}, name=label,
                hovertemplate=f"{label}<extra></extra>", showlegend=False,
            ))
    figure.update_layout(title=_text(title))
    return figure


def feature_map_figure(features: pd.DataFrame, title: str) -> go.Figure:
    """Show every LFQ feature independently, with exact native intensity values.

    A direct scatter retains coincident features and their label/charge hover
    data. Constant marker size keeps zero-intensity features visible.
    """
    frame = _frame(features, ("rt", "mz", "intensity", "charge"), "LFQ features")
    if "label" not in frame:
        raise ValueError("LFQ features are missing RNA labels.")
    if (frame["charge"] % 1 != 0).any():
        raise ValueError("LFQ feature charges must be integers.")
    labels = frame["label"].fillna("").map(_text)
    figure = go.Figure(go.Scattergl(
        x=frame["rt"], y=frame["mz"], mode="markers",
        customdata=np.column_stack([labels, frame["charge"], frame["intensity"]]),
        marker={
            "size": 8, "color": frame["intensity"], "colorscale": "Viridis",
            "cmin": 0, "showscale": True, "colorbar": {"title": "Feature intensity"},
        },
        hovertemplate="RNA: %{customdata[0]}<br>Charge: %{customdata[1]}<br>Retention time: %{x} s<br>m/z: %{y}<br>Feature intensity: %{customdata[2]}<extra></extra>",
    ))
    figure.update_layout(
        title=_text(title), xaxis_title="Retention time (s)", yaxis_title="m/z",
        height=500,
    )
    return figure
