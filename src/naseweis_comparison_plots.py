"""Relative RNA abundance plots built from saved condition comparisons.

Figures preserve biological-replicate means produced by ``compare_conditions``.
They do not impute absent measurements, add pseudocounts or infer modifications.
"""

from html import escape

import numpy as np
import pandas as pd
import polars as pl

from src.naseweis_plots import Plot
from src.naseweis_results import nucleotide_token, split_sequence
from src.viz.quantification import condition_scatter as scatter

MEASURE_LABELS = {
    "intensity": "Native feature area",
    "model_area": "Fitted model area",
}


def _text(value):
    """Keep result names and imported condition names literal in plot labels."""
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
) -> Plot:
    """Denominator means on x and numerator means on y, without rescaling, as
    the visualization library's condition scatter (an Insight ``SeriesPlot``).

    The inputs are saved biological-replicate summaries, already averaged over
    technical acquisitions. A ``y=x`` line marks equal abundance; on the log
    scale both axes are log10 means. ``labels`` picks exact RNA identities to
    draw apart (``Selected``); it changes display only. Terminal
    phosphate/cyclic-phosphate groups remain unmodified. D/I/Y and bracketed
    nucleotide annotations are drawn as modified/ambiguous.

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
    log = scale == "log"
    entries = []
    for number, (label, row) in enumerate(frame.iterrows()):
        stats = [_stat(row, field) for field in ("denominator_n", "numerator_n", "ratio", "log2_fold_change")]
        hover = (
            f"RNA: {_text(label)} · {_text(denominator)} mean {row['denominator_mean']:.6g} "
            f"(n={_number(stats[0])}) · {_text(numerator)} mean {row['numerator_mean']:.6g} "
            f"(n={_number(stats[1])}) · ratio {_number(stats[2])} · log2 FC {_number(stats[3])}"
        )
        entries.append({"entry": number, "x": float(row["denominator_mean"]),
                        "y": float(row["numerator_mean"]), "modified": _modified(label),
                        "label": hover})
    data = scatter.rows(entries, "naseweis", log=log)
    picked = [number for number, label in enumerate(frame.index) if label in selected]
    if picked:
        data = data.with_columns(series=pl.when(pl.col("entry").is_in(picked))
                                 .then(pl.lit("Selected")).otherwise(pl.col("series")))
    units = "relative RNA abundance (total-area normalized)" if normalization == "total" else MEASURE_LABELS[measure].lower()
    prefix = "log10 mean" if log else "Mean"
    series = dict(scatter.SERIES) | {"Selected": "#ff7501"}
    names = {"Unmodified": "Unmodified", "Modified": "Modified / ambiguous",
             "Equal": "Equal abundance (y = x)", "Selected": "Labelled"}
    settings = scatter.settings(
        _text(numerator), _text(denominator), log=log, filters={}, interactivity={}, highlight={},
        series_names=names, series_colors=series,
        title="Label-free quantification of RNA oligonucleotides",
        x_label=f"{prefix} {_text(denominator)} · {units}",
        y_label=f"{prefix} {_text(numerator)} · {units}",
        svg_filename="naseweis-lfq-comparison",
    )
    table = frame.assign(modified=[_modified(label) for label in frame.index],
                         labelled=[label in selected for label in frame.index])
    table = table.rename_axis("rna").reset_index()
    return Plot("SeriesPlot", data, settings, table, 650)


def _number(value):
    return "–" if value is None else f"{value:.4g}"
