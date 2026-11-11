"""Relative RNA abundance figures built from saved condition comparisons.

Figures preserve biological-replicate means produced by ``compare_conditions``.
They do not impute absent measurements, add pseudocounts or infer modifications.
"""

from html import escape
from sys import float_info

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from src.naseweis_results import nucleotide_token, split_sequence

MEASURE_LABELS = {
    "intensity": "Native feature area",
    "model_area": "Fitted model area",
}


def _text(value):
    """Keep result names and imported condition names literal in Plotly HTML."""
    return escape(str(value), quote=True).replace("$", "&#36;")


def paired_condition_rows(contrasts: pd.DataFrame) -> pd.DataFrame:
    """Return a copy containing positive finite means measured in both groups.

    Missing, zero, negative and nonfinite means are excluded from either display
    scale. No value is substituted for an unpaired measurement. Other columns,
    index labels and their order are preserved for plotting or TSV downloads.
    """
    if not isinstance(contrasts, pd.DataFrame):
        raise ValueError("Condition comparisons must be a table.")  # noqa: TRY004
    required = {"numerator_mean", "denominator_mean"}
    if not required <= set(contrasts):
        raise ValueError("Condition comparisons need numerator and denominator means.")
    if not contrasts.index.is_unique or not all(
        isinstance(label, str) and label.strip() for label in contrasts.index
    ):
        raise ValueError("Condition comparisons need unique, nonempty RNA labels.")
    frame = contrasts.copy(deep=True)
    paired = pd.Series(True, index=frame.index)
    for column in ("numerator_mean", "denominator_mean"):
        try:
            values = pd.to_numeric(frame[column], errors="raise").astype("float64")
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(f"Condition comparison {column} must be numeric.") from error
        paired &= np.isfinite(values) & values.gt(0)
        frame[column] = values
    return frame.loc[paired].copy()


def _modified(label):
    # Ignore terminal chemistry; D/I/Y and bracketed bases are native annotated
    # modifications or ambiguity, not proof
    # of a particular chemical identity or uniquely localized modification site.
    return any(nucleotide_token(token) not in {"A", "C", "G", "U", "N"}
               for token in split_sequence(label))


def _stat(row, field):
    """Return finite optional hover metadata without emitting NaN in HTML."""
    value = row.get(field)
    if value is None or pd.isna(value):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"Condition comparison {field} must be numeric.") from error
    return numeric if np.isfinite(numeric) else None


def condition_scatter(
    contrasts: pd.DataFrame,
    numerator: str,
    denominator: str,
    *,
    measure="intensity",
    normalization="total",
    scale="log",
    labels=(),
) -> go.Figure:
    """Plot denominator means on x and numerator means on y without rescaling.

    The inputs are saved biological-replicate summaries, already averaged over
    technical acquisitions. Both axes share a range and scale, with a ``y=x``
    reference line. ``labels`` selects exact RNA identities to label; it changes
    display only. Terminal phosphate/cyclic-phosphate groups remain unmodified
    markers. D/I/Y and bracketed nucleotide annotations use the modified/ambiguous marker.

    ``paired_condition_rows`` exposes exactly which rows can appear. If no RNA
    has a positive finite mean in both conditions, this raises ``ValueError``
    rather than drawing invented values at zero.
    """
    if not all(isinstance(value, str) and value.strip() for value in (numerator, denominator)):
        raise ValueError("Choose named numerator and denominator conditions.")
    if numerator == denominator:
        raise ValueError("Choose two different conditions for the scatter plot.")
    if measure not in MEASURE_LABELS or normalization not in ("none", "total"):
        raise ValueError("Unsupported RNA area measure or normalization for plotting.")
    if scale not in ("log", "linear"):
        raise ValueError("Scatter plot scale must be log or linear.")
    if isinstance(labels, str) or not all(isinstance(label, str) for label in labels):
        raise ValueError("Plot labels must be a collection of RNA identities.")
    selected = set(labels)
    frame = paired_condition_rows(contrasts)
    if frame.empty:
        raise ValueError("No RNA has a positive measured mean in both conditions.")
    modified = pd.Series([_modified(label) for label in frame.index], index=frame.index)
    values = frame[["denominator_mean", "numerator_mean"]].to_numpy()
    minimum, maximum = float(values.min()), float(values.max())
    if scale == "log":
        low, high = np.log10([minimum, maximum])
        padding = max(float(high - low), 1.0) * 0.06
        axis_range = [float(low - padding), float(high + padding)]
        line_limits = [minimum, maximum]
        if minimum == maximum:
            # A single level still gets a visible diagonal; avoid overflow and
            # underflow at the extremes of the native floating-point range.
            line_limits = [
                minimum / 2 if minimum / 2 > 0 else minimum,
                maximum * 2 if maximum <= float_info.max / 2 else maximum,
            ]
    else:
        high = maximum * 1.06 if maximum <= float_info.max / 1.06 else maximum
        axis_range, line_limits = [0.0, high], [0.0, maximum]
    figure = go.Figure()
    figure.add_trace(go.Scatter(
        x=line_limits, y=line_limits, mode="lines", name="Equal abundance (y = x)",
        line={"color": "#9CA3AF", "dash": "dot", "width": 1}, hoverinfo="skip",
    ))
    for is_modified, name, symbol in (
        (False, "Unmodified", "circle-open"),
        (True, "Modified / ambiguous", "circle"),
    ):
        subset = frame.loc[modified == is_modified]
        if subset.empty:
            continue
        names = [_text(label) for label in subset.index]
        customdata = [
            [
                _text(label), _stat(row, "numerator_n"), _stat(row, "denominator_n"),
                _stat(row, "ratio"), _stat(row, "log2_fold_change"),
            ]
            for label, row in subset.iterrows()
        ]
        texts = [name if label in selected else "" for label, name in zip(subset.index, names)]
        figure.add_trace(go.Scatter(
            x=subset["denominator_mean"].tolist(),
            y=subset["numerator_mean"].tolist(),
            mode="markers+text" if any(texts) else "markers", name=name,
            marker={
                "symbol": symbol,
                "color": ["#D55E00" if label in selected else "#29379b" for label in subset.index],
                "line": {"width": 1.5},
                "size": [11 if label in selected else 8 for label in subset.index],
            },
            text=texts, textposition="top center", customdata=customdata,
            hovertemplate=(
                "RNA: %{customdata[0]}<br>"
                f"{_text(denominator)} mean: %{{x:.6g}}<br>"
                f"{_text(numerator)} mean: %{{y:.6g}}<br>"
                f"{_text(denominator)} biological replicates: %{{customdata[2]}}<br>"
                f"{_text(numerator)} biological replicates: %{{customdata[1]}}<br>"
                "Ratio (y / x): %{customdata[3]:.4g}<br>"
                "Log2 fold change: %{customdata[4]:.4g}<extra>%{fullData.name}</extra>"
            ),
        ))
    units = "Relative RNA abundance (total-area normalized)" if normalization == "total" else MEASURE_LABELS[measure]
    figure.update_layout(
        title="Label-free quantification of RNA oligonucleotides",
        xaxis={"title": f"{_text(denominator)} · {units}", "type": scale, "range": axis_range, "dtick": 1 if scale == "log" else None, "constrain": "domain"},
        yaxis={"title": f"{_text(numerator)} · {units}", "type": scale, "range": axis_range, "dtick": 1 if scale == "log" else None, "scaleanchor": "x", "scaleratio": 1, "constrain": "domain"},
        height=650, legend={"orientation": "h", "y": 1.1, "x": 0},
        margin={"l": 85, "r": 35, "t": 115, "b": 85},
    )
    return figure
