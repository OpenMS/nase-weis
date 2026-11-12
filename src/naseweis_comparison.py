"""RNA LFQ matrices and relative comparisons for a saved experiment.

Following the recovered RNA KNIME workflow, native adduct forms and charges
are summed by modified RNA label, then samples can be normalized by total area.
Missing measurements stay missing; this module does not match features by RT.
"""

# The public validation contract uses ValueError for actionable UI/worker errors.
# ruff: noqa: TRY004

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.naseweis_lfq import read_lfq_table

MEASURES = {"intensity", "model_area"}
NORMALIZATIONS = {"none", "total"}
DESIGN_COLUMNS = ["file", "path", "condition", "replicate", "fraction", "status"]
COMPARISON_COLUMNS = [
    "numerator_mean",
    "numerator_std",
    "numerator_n",
    "denominator_mean",
    "denominator_std",
    "denominator_n",
    "ratio",
    "log2_fold_change",
]
_SAMPLE_ID = re.compile(r"sample_[0-9]{3,}")
_ADDUCT_SUFFIX = re.compile(r"\+\[[^\[\]]+\]$")


def _positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"Sample {name} must be a positive integer.")
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"Sample {name} must be a positive integer.") from error
    if not np.isfinite(number) or not number.is_integer() or number < 1:
        raise ValueError(f"Sample {name} must be a positive integer.")
    return int(number)


def _sample_design(samples):
    if not isinstance(samples, list) or not samples:
        raise ValueError("Select at least one sample for the LFQ experiment.")
    rows, seen = [], set()
    for sample in samples:
        if not isinstance(sample, dict):
            raise ValueError("Every LFQ sample must have a design record.")
        identity = sample.get("id")
        if not isinstance(identity, str) or not _SAMPLE_ID.fullmatch(identity):
            raise ValueError("LFQ sample IDs must have the controlled form sample_001.")
        if identity in seen:
            raise ValueError(f"Duplicate LFQ sample ID: {identity}.")
        seen.add(identity)
        if not all(
            isinstance(sample.get(key), str) and sample[key]
            for key in ("file", "path", "status")
        ):
            raise ValueError(
                f"Sample {identity} needs its file, path and execution status."
            )
        condition = sample.get("condition", "")
        if not isinstance(condition, str):
            raise ValueError(f"Sample {identity} condition must be text.")
        rows.append(
            {
                "sample_id": identity,
                "file": sample["file"],
                "path": sample["path"],
                "condition": condition.strip(),
                "replicate": _positive_integer(sample.get("replicate"), "replicate"),
                "fraction": _positive_integer(sample.get("fraction", 1), "fraction"),
                "status": sample["status"],
            }
        )
    return pd.DataFrame(rows).set_index("sample_id")[DESIGN_COLUMNS]


def _areas(values, context, *, missing=False):
    """Native areas must be finite/nonnegative; matrix gaps alone may be NaN."""
    try:
        numeric = pd.to_numeric(values, errors="raise").astype("float64")
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{context} contains invalid numeric areas.") from error
    valid = np.isfinite(numeric)
    if missing:
        valid |= numeric.isna()
    if not valid.all() or (numeric < 0).any():
        raise ValueError(f"{context} areas must be finite and nonnegative.")
    return numeric


def _rna_areas(table, identity, measure, formulas, warnings):
    if table.empty:
        return pd.Series(dtype="float64", name=identity)
    if "label" not in table or measure not in table:
        raise ValueError(
            f"Sample {identity} requires native label and {measure} columns; a different area measure cannot be substituted."
        )
    areas = _areas(table[measure], f"Sample {identity} {measure}")
    labels = table["label"]
    if not labels.map(
        lambda label: isinstance(label, str) and bool(label.strip())
    ).all():
        raise ValueError(
            f"Sample {identity} has a missing or invalid native RNA label."
        )
    cleaned = labels.map(lambda label: _ADDUCT_SUFFIX.sub("", label))
    if not cleaned.map(lambda label: bool(label.strip())).all():
        raise ValueError(
            f"Sample {identity} has an empty RNA label after removing its adduct suffix."
        )
    known = table.get("sum_formula", pd.Series("", index=table.index))
    incomplete = False
    for label, formula in zip(labels, known):
        if pd.isna(formula) or formula == "":
            incomplete = True
            continue
        if not isinstance(formula, str):
            raise ValueError(
                f"Sample {identity} has an invalid sum formula for RNA label {label}."
            )
        if label in formulas and formulas[label] != formula:
            raise ValueError(
                f"Conflicting sum formulas for native RNA label {label}; those identities cannot be pooled."
            )
        formulas[label] = formula
    if incomplete:
        warnings.append(
            f"Sample {identity} has incomplete formula metadata; those identities are pooled by native label only."
        )
    # Different adduct labels legitimately have different sum formulas. Check
    # chemistry per full native label before collapsing charge/adduct forms.
    positive = pd.DataFrame({"label": cleaned, "area": areas}).loc[areas > 0]
    with np.errstate(over="ignore", invalid="ignore"):
        pooled = positive.groupby("label", sort=True)["area"].sum()
    if not np.isfinite(pooled).all():
        raise ValueError(
            f"Sample {identity} pooled RNA areas overflowed the numeric range."
        )
    return pooled.rename(identity)


def build_comparison(
    results_dir: Path,
    samples: list[dict],
    *,
    measure="intensity",
    normalization="total",
) -> dict:
    """Build and save raw/relative RNA matrices from finished sample outputs.

    Fitted ``model_area`` is an explicit alternative to native feature intensity.
    Native intensity can contain an imputed estimate after failed model fitting;
    it is not silently substituted for a requested fitted area.
    """
    if (
        not isinstance(measure, str)
        or measure not in MEASURES
        or not isinstance(normalization, str)
        or normalization not in NORMALIZATIONS
    ):
        raise ValueError(
            "LFQ measure must be intensity or model_area; normalization must be total or none."
        )
    design = _sample_design(samples)
    root = Path(results_dir)
    columns, formulas, warnings = [], {}, []
    for identity, sample in design.iterrows():
        if sample["status"] != "finished":
            warnings.append(
                f"Sample {identity} did not finish; its measurements remain missing."
            )
            column = pd.Series(dtype="float64", name=identity)
        else:
            table_path = (
                root
                / "samples"
                / identity
                / "results"
                / "quantification"
                / "features.tsv"
            )
            if not table_path.resolve().is_relative_to(root.resolve()):
                raise ValueError(
                    f"Sample {identity} output is outside this experiment."
                )
            column = _rna_areas(
                read_lfq_table(table_path), identity, measure, formulas, warnings
            )
        columns.append(column)
    raw = pd.concat(columns, axis=1).sort_index()
    raw.index.name, raw.columns.name = "label", "sample_id"
    raw = raw.astype("float64")
    abundances = raw.copy(deep=True)
    factors, totals = {}, {}
    for identity in raw:
        with np.errstate(over="ignore", invalid="ignore"):
            total = float(raw[identity].sum())
        if not np.isfinite(total):
            raise ValueError(
                f"Sample {identity} total RNA area overflowed the numeric range."
            )
        totals[identity] = total
        factor = total if normalization == "total" and total > 0 else None
        if normalization == "none":
            factor = 1.0
        elif factor is None:
            abundances[identity] = np.nan
            warnings.append(
                f"Sample {identity} has no positive quantified area; total normalization is unavailable."
            )
        else:
            abundances[identity] = raw[identity] / factor
        factors[identity] = factor
    metadata = {
        "schema_version": 1,
        "measure": measure,
        "normalization": normalization,
        "aggregation": "sum_by_modified_rna_label_across_charges_and_adducts",
        "adduct_suffix": "trailing +[adduct] only",
        "missing_values": "not_quantified",
        "zero_area_features": "excluded",
        "sample_ids": list(design.index),
        "sample_totals": totals,
        "factors": factors,
        "warnings": warnings,
    }
    output = root / "comparison"
    output.mkdir(parents=True, exist_ok=True)
    raw.to_csv(output / "raw-abundances.tsv", sep="\t", float_format="%.17g")
    abundances.to_csv(output / "abundances.tsv", sep="\t", float_format="%.17g")
    design.to_csv(output / "sample-design.tsv", sep="\t")
    (output / "comparison.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return {
        "raw": raw,
        "abundances": abundances,
        "design": design,
        "metadata": metadata,
    }


def _read_matrix(path):
    try:
        table = pd.read_csv(
            path, sep="\t", index_col=0, keep_default_na=False, dtype=str
        )
    except (ValueError, pd.errors.ParserError) as error:
        raise ValueError(f"Cannot read the saved LFQ matrix {path.name}.") from error
    if (
        table.index.name != "label"
        or not table.index.is_unique
        or any(not str(label) for label in table.index)
    ):
        raise ValueError(f"Invalid RNA identities in saved LFQ matrix {path.name}.")
    for identity in table:
        table[identity] = _areas(
            table[identity].replace("", np.nan),
            f"Saved sample {identity}",
            missing=True,
        )
    table.columns.name = "sample_id"
    return table


def load_comparison(results_dir: Path) -> dict:
    """Restore the same matrices/design/metadata returned by build_comparison."""
    output = Path(results_dir) / "comparison"
    try:
        metadata = json.loads((output / "comparison.json").read_text(encoding="utf-8"))
    except ValueError as error:
        raise ValueError("The saved LFQ comparison metadata is invalid.") from error
    if (
        not isinstance(metadata, dict)
        or metadata.get("schema_version") != 1
        or not isinstance(metadata.get("measure"), str)
        or metadata["measure"] not in MEASURES
        or not isinstance(metadata.get("normalization"), str)
        or metadata["normalization"] not in NORMALIZATIONS
    ):
        raise ValueError(
            "The saved LFQ comparison has an unsupported schema or method."
        )
    try:
        stored_design = pd.read_csv(
            output / "sample-design.tsv", sep="\t", keep_default_na=False, dtype=str
        )
        records = stored_design.rename(columns={"sample_id": "id"}).to_dict("records")
        design = _sample_design(records)
    except (KeyError, ValueError, pd.errors.ParserError) as error:
        raise ValueError("The saved LFQ sample design is invalid.") from error
    raw, abundances = (
        _read_matrix(output / "raw-abundances.tsv"),
        _read_matrix(output / "abundances.tsv"),
    )
    identities = list(design.index)
    if (
        metadata.get("sample_ids") != identities
        or list(raw.columns) != identities
        or list(abundances.columns) != identities
        or not raw.index.equals(abundances.index)
    ):
        raise ValueError(
            "The saved LFQ matrices and sample design do not describe the same experiment."
        )
    return {
        "raw": raw,
        "abundances": abundances,
        "design": design,
        "metadata": metadata,
    }


def compare_conditions(matrix, design, numerator, denominator):
    """Compare biological-replicate means after averaging technical injections.

    The experiment preflight must reject fractionated designs for this initial
    workflow. Standard deviations use biological replicates and ddof=1. Ratios
    and log2 changes require positive means in both conditions; no pseudocount,
    imputation or significance test is applied.
    """
    if not isinstance(matrix, pd.DataFrame) or not isinstance(design, pd.DataFrame):
        raise ValueError(
            "Condition comparison requires an abundance matrix and sample design."
        )
    if (
        not matrix.columns.is_unique
        or not design.index.is_unique
        or not {"condition", "replicate"} <= set(design)
    ):
        raise ValueError(
            "Condition comparison requires a unique sample design with condition and replicate columns."
        )
    if not set(matrix.columns) <= set(design.index):
        raise ValueError("The sample design does not cover every abundance column.")
    aligned = design.loc[matrix.columns].copy()
    aligned["replicate"] = aligned["replicate"].map(
        lambda value: _positive_integer(value, "replicate")
    )
    available = set(aligned["condition"])
    if (
        not numerator
        or not denominator
        or numerator == denominator
        or numerator not in available
        or denominator not in available
    ):
        raise ValueError(
            "Choose two different conditions represented in this experiment."
        )
    numeric = matrix.copy(deep=True)
    for identity in numeric:
        numeric[identity] = _areas(
            numeric[identity], f"Sample {identity}", missing=True
        )
    result = pd.DataFrame(index=matrix.index)
    for prefix, condition in (("numerator", numerator), ("denominator", denominator)):
        members = aligned.loc[aligned["condition"] == condition]
        biological = pd.concat(
            [
                numeric.loc[:, group.index].mean(axis=1).rename(replicate)
                for replicate, group in members.groupby("replicate", sort=True)
            ],
            axis=1,
        )
        result[f"{prefix}_mean"] = biological.mean(axis=1)
        result[f"{prefix}_std"] = biological.std(axis=1, ddof=1)
        result[f"{prefix}_n"] = biological.count(axis=1)
    positive = (result["numerator_mean"] > 0) & (result["denominator_mean"] > 0)
    result["ratio"] = (result["numerator_mean"] / result["denominator_mean"]).where(
        positive
    )
    result["log2_fold_change"] = np.log2(result["ratio"])
    return result[COMPARISON_COLUMNS]
