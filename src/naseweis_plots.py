"""Pure, offline plot specifications for NASE spectra, acquisitions and RNA features.

Every plot is an OpenMS-Insight component from the template's visualization
library: a ``Plot`` names the component, the rows it draws and its settings.
``src.naseweis_insight.show_plot`` draws it. Nothing here imports Streamlit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from html import escape

import numpy as np
import pandas as pd
import polars as pl

from src.viz.identification import fragment_spectrum
from src.viz.quantification import feature_xic
from src.viz.identification.fragment_spectrum import ION_COLORS, UNKNOWN_ION_COLOR

# The picked colour of the template's plot style
SELECTED_ION_COLOR = "#ff7501"
TRACE_COLORS = {"TIC": "#0072B2", "BPC": "#009E73", "XIC": "#D55E00"}

__all__ = [
    "ION_COLORS", "UNKNOWN_ION_COLOR", "SELECTED_ION_COLOR", "TRACE_COLORS", "Plot",
    "ion_series", "spectrum_plot", "chromatogram_plot", "peakmap_plots", "feature_map_plot",
    "feature_xic_plot", "abundance_bars_plot",
]


@dataclass(frozen=True)
class Plot:
    """One Insight component: its class name in ``openms_insight``, the rows it
    draws, its constructor arguments and the table offered as a download."""

    component: str
    data: pl.DataFrame
    settings: dict
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    height: int = 380


def _text(value) -> str:
    """Treat labels as text, including dollar signs some renderers read as math."""
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
    return fragment_spectrum.ion_series(ion) or ""


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


def spectrum_plot(peaks: pd.DataFrame, title: str, *, highlighted_mz=None) -> Plot:
    """Unbinned peaks as a fragment spectrum, assigned ions labelled and
    coloured by their RNA ion series.

    Intensities are shown relative to the spectrum's maximum; only the display
    changes. Highlights identify exact, unambiguous acquired peak masses; no
    peak is snapped to a nearby mass or relabelled as an assigned ion.
    """
    frame = _frame(peaks, ("mz", "intensity"), "Spectrum peaks")
    if frame["intensity"].max() == 0:
        raise ValueError("Spectrum peaks has no positive intensities to normalize.")
    frame["ion"] = frame.get("ion", pd.Series("", index=frame.index)).fillna("").astype(str)
    frame = frame[["mz", "intensity", "ion"]].sort_values("mz", kind="stable").reset_index(drop=True)
    selected = _highlighted_peaks(frame, highlighted_mz)
    ions = {index: _text(ion) for index, ion in enumerate(frame["ion"]) if ion}
    data = fragment_spectrum.rows(zip(frame["mz"], frame["intensity"]), "spectrum", ions)
    if not selected.empty:
        picked = pl.Series(frame["mz"].isin(selected["mz"]).to_list())
        data = data.with_columns(
            matched=pl.col("matched") | picked,
            color=pl.when(picked).then(pl.lit(SELECTED_ION_COLOR)).otherwise(pl.col("color")),
        )
    settings = fragment_spectrum.settings(filters={}, interactivity={}, title=_text(title))
    table = frame.assign(relative_intensity=data["intensity"].to_list(),
                         selected=frame["mz"].isin(selected["mz"]))
    return Plot("LinePlot", data, settings, table, 380)


def chromatogram_plot(points: pd.DataFrame, title: str) -> Plot:
    """Overlay TIC, BPC, XIC or isotope traces without changing intensities."""
    frame = _frame(points, ("rt", "intensity"), "Chromatogram points")
    if "trace" not in frame or frame["trace"].isna().any():
        raise ValueError("Chromatogram points must identify each trace.")
    frame["trace"] = frame["trace"].astype(str)
    shown = {trace: _text(trace) for trace in frame["trace"].unique()}
    frame = frame.sort_values(["trace", "rt"], kind="stable").reset_index(drop=True)
    traces = list(dict.fromkeys(frame["trace"]))
    palette = ["#CC79A7", "#56B4E9", "#E69F00", "#F0E442", "#999999"]
    colors = {shown[trace]: TRACE_COLORS.get(trace, palette[i % len(palette)])
              for i, trace in enumerate(traces)}
    data = pl.DataFrame({"rt": frame["rt"].astype(float).tolist(),
                         "intensity": frame["intensity"].astype(float).tolist(),
                         "series": [shown[trace] for trace in frame["trace"]]})
    settings = dict(
        x_column="rt", y_column="intensity", series_column="series",
        series_names={shown[t]: shown[t] for t in traces}, series_colors=colors,
        mode="lines", y_range_mode="nonnegative", title=_text(title),
        x_label="Retention time (s)", y_label="Intensity", svg_filename="naseweis-chromatograms",
    )
    return Plot("SeriesPlot", data, settings, frame[["trace", "rt", "intensity"]], 340)


def _marginal(frame: pd.DataFrame, axis: str, label: str, title: str) -> Plot:
    summed = frame.groupby(axis, sort=True)["intensity"].sum().reset_index()
    data = pl.DataFrame({"x": summed[axis].astype(float).tolist(),
                         "y": summed["intensity"].astype(float).tolist(),
                         "series": ["Summed intensity"] * len(summed)})
    settings = dict(
        x_column="x", y_column="y", series_column="series",
        series_names={"Summed intensity": "Summed intensity"},
        series_colors={"Summed intensity": "#0072B2"}, mode="lines",
        y_range_mode="nonnegative", title=title, x_label=label, y_label="Summed intensity",
        svg_filename=f"naseweis-ms1-{axis}-marginal",
    )
    return Plot("SeriesPlot", data, settings, summed, 240)


def peakmap_plots(points: pd.DataFrame, title: str, marginals: bool = False) -> list[Plot]:
    """Prepared MS1 points as a peak map, with the summed signal over
    retention time and m/z below it when ``marginals``.

    The caller owns any overview binning.
    """
    frame = _frame(points, ("rt", "mz", "intensity"), "Peak-map points")
    data = pl.DataFrame({"rt": frame["rt"].astype(float).tolist(),
                         "mz": frame["mz"].astype(float).tolist(),
                         "intensity": frame["intensity"].astype(float).tolist()})
    settings = dict(
        x_column="rt", y_column="mz", intensity_column="intensity", title=_text(title),
        x_label="Retention time (s)", y_label="m/z", intensity_label="Intensity",
        colorscale="Viridis", svg_filename="naseweis-ms1-peakmap",
    )
    plots = [Plot("Heatmap", data, settings, frame[["rt", "mz", "intensity"]], 440)]
    if marginals:
        plots.append(_marginal(frame, "rt", "Retention time (s)", "Signal over retention time"))
        plots.append(_marginal(frame, "mz", "m/z", "Signal over m/z"))
    return plots


def feature_map_plot(features: pd.DataFrame, title: str) -> Plot:
    """Every LFQ feature independently, with exact native intensity values and
    its RNA label and charge on hover."""
    frame = _frame(features, ("rt", "mz", "intensity", "charge"), "LFQ features")
    if "label" not in frame:
        raise ValueError("LFQ features are missing RNA labels.")
    if (frame["charge"] % 1 != 0).any():
        raise ValueError("LFQ feature charges must be integers.")
    labels = frame["label"].fillna("").astype(str)
    hover = [f"RNA: {_text(label)} · charge {int(charge)} · intensity {value:.6g}"
             for label, charge, value in zip(labels, frame["charge"], frame["intensity"])]
    data = pl.DataFrame({"feature": list(range(len(frame))),
                         "rt": frame["rt"].astype(float).tolist(),
                         "mz": frame["mz"].astype(float).tolist(),
                         "intensity": frame["intensity"].astype(float).tolist(),
                         "label": hover})
    settings = dict(
        x_column="rt", y_column="mz", intensity_column="intensity", hover_column="label",
        title=_text(title), x_label="Retention time (s)", y_label="m/z",
        intensity_label="Feature intensity", colorscale="Viridis", min_points=len(frame),
        svg_filename="naseweis-lfq-features",
    )
    table = frame[["rt", "mz", "intensity", "charge"]].assign(label=labels)
    return Plot("Heatmap", data, settings, table, 500)


def feature_xic_plot(trace: pd.DataFrame, label: str, boundaries=None) -> Plot:
    """One LFQ feature's raw MS1 XIC, with the part inside its native feature
    boundaries drawn as its own series. Retention times stay in seconds."""
    frame = _frame(trace, ("rt", "intensity"), "Chromatogram points")
    frame = frame.sort_values("rt", kind="stable").reset_index(drop=True)
    start, end = (boundaries or (None, None))[:2]
    data = feature_xic.rows(_text(label), frame["rt"], frame["intensity"], start, end)
    settings = feature_xic.settings(
        filters={}, interactivity={}, title=_text(label), x_label="Retention time (s)",
        svg_filename="naseweis-lfq-xic",
    )
    table = frame[["rt", "intensity"]].assign(
        in_feature=frame["rt"].between(start, end) if start is not None and end is not None else False
    )
    return Plot("SeriesPlot", data, settings, table, 340)


def abundance_bars_plot(labels, values, title: str, y_label: str) -> Plot:
    """One RNA target's abundance per acquisition as bars; a missing
    measurement stays missing rather than drawn as zero."""
    labels = [str(label) for label in labels]
    numbers = [None if value is None or pd.isna(value) else float(value) for value in values]
    if len(labels) != len(numbers) or not labels:
        raise ValueError("Abundance bars need one value per acquisition.")
    data = pl.DataFrame({"acquisition": labels, "abundance": numbers,
                         "series": ["Abundance"] * len(labels)},
                        schema={"acquisition": pl.Utf8, "abundance": pl.Float64, "series": pl.Utf8})
    settings = dict(
        x_column="acquisition", y_column="abundance", series_column="series",
        series_names={"Abundance": _text(title)}, series_colors={"Abundance": "#29379b"},
        mode="bars", y_range_mode="nonnegative", title=_text(title), x_label="Acquisition",
        y_label=y_label, svg_filename="naseweis-abundance",
    )
    table = pd.DataFrame({"acquisition": labels, "abundance": numbers})
    return Plot("SeriesPlot", data, settings, table, 340)
