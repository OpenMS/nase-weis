"""Native RNA fragment assignments for the existing spectrum viewer.

Assignments come from NASE's idXML output, without regenerating theoretical
spectra or guessing spectrum identity from retention times. Peak heights always
come from the spectrum being viewed, because NASE stores scored intensities.
Invalid viewer inputs consistently raise ValueError for the UI's error path.
"""

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

ANNOTATION_COLUMNS = [
    "native_id",
    "sequence",
    "charge",
    "mz",
    "intensity",
    "ion",
    "fragment_charge",
    "identification_index",
    "hit_index",
]


def native_spectrum_id(spectra_ref: str) -> str:
    """Remove one mzTab run prefix while preserving the complete native ID."""
    if not isinstance(spectra_ref, str):
        raise ValueError("A spectrum reference must be text containing a native ID.")  # noqa: TRY004
    return re.sub(r"^ms_run\[[0-9]+\]:", "", spectra_ref, count=1)


def _integer(value, label: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be an integer charge.")  # noqa: TRY004
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{label} must be an integer charge.") from error
    if not math.isfinite(number) or not number.is_integer():
        raise ValueError(f"{label} must be an integer charge.")
    return int(number)


def _finite_number(value, label: str, *, positive=False) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be a finite number.")  # noqa: TRY004
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{label} must be a finite number.") from error
    if not math.isfinite(number) or (positive and number <= 0):
        suffix = "positive finite" if positive else "finite"
        raise ValueError(f"{label} must be a {suffix} number.")
    return number


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    if not isinstance(frame, pd.DataFrame):
        raise ValueError(f"{label} must be a pandas DataFrame.")  # noqa: TRY004
    missing = sorted(required - set(frame.columns))
    if missing or not frame.columns.is_unique:
        detail = (
            f"missing columns {', '.join(missing)}" if missing else "duplicate columns"
        )
        raise ValueError(f"Invalid {label}: {detail}.")


def read_annotations(idxml_path: str | Path) -> pd.DataFrame:
    """Read all native NASE RNA hit annotations with their exact provenance.

    RNA sequences are stored in the hit's ``label`` metadata, not its peptide
    sequence field. The two provenance columns distinguish different native
    records with otherwise identical spectrum, sequence and signed charge.
    """
    import pyopenms as oms

    path = Path(idxml_path)
    if not path.exists():
        raise FileNotFoundError(f"NASE idXML annotations are unavailable: {path}")
    if not path.is_file():
        raise OSError(f"Select a NASE idXML file for RNA annotations: {path}")
    try:
        _, identifications = oms.IdXMLFile().load(str(path))
    except Exception as error:
        raise ValueError(
            f"Cannot read RNA annotations from {path.name}; use an intact NASE idXML output."
        ) from error

    rows = []
    for identification_index, identification in enumerate(identifications):
        for hit_index, hit in enumerate(identification.getHits()):
            annotations = hit.getPeakAnnotations()
            if not annotations:
                continue
            context = f"identification {identification_index + 1}, hit {hit_index + 1}"
            native_id = native_spectrum_id(identification.getSpectrumReference())
            if not native_id:
                raise ValueError(
                    f"Missing native spectrum ID in {context} of {path.name}."
                )
            if not hit.metaValueExists("label"):
                raise ValueError(
                    f"Missing RNA label in {context} of {path.name}; use NASE idXML output."
                )
            sequence = hit.getMetaValue("label")
            if not isinstance(sequence, str) or not sequence.strip():
                raise ValueError(
                    f"Empty or invalid RNA label in {context} of {path.name}."
                )
            if (
                hit.metaValueExists("molecule_type")
                and hit.getMetaValue("molecule_type") != "RNA"
            ):
                raise ValueError(
                    f"The annotated hit in {context} of {path.name} is not RNA."
                )
            charge = _integer(hit.getCharge(), "RNA precursor charge")
            if charge == 0:
                raise ValueError(
                    f"Missing nonzero RNA precursor charge in {context} of {path.name}."
                )
            for annotation in annotations:
                if (
                    not isinstance(annotation.annotation, str)
                    or not annotation.annotation.strip()
                ):
                    raise ValueError(
                        f"Missing fragment ion label in {context} of {path.name}."
                    )
                rows.append(
                    {
                        "native_id": native_id,
                        "sequence": sequence,
                        "charge": charge,
                        "mz": _finite_number(
                            annotation.mz, "Fragment m/z", positive=True
                        ),
                        "intensity": _finite_number(
                            annotation.intensity, "Scored fragment intensity"
                        ),
                        "ion": annotation.annotation,
                        "fragment_charge": _integer(
                            annotation.charge, "Fragment charge"
                        ),
                        "identification_index": identification_index,
                        "hit_index": hit_index,
                    }
                )
    return pd.DataFrame(rows, columns=ANNOTATION_COLUMNS).astype(
        {
            "charge": "int64",
            "mz": "float64",
            "intensity": "float64",
            "fragment_charge": "int64",
            "identification_index": "int64",
            "hit_index": "int64",
        }
    )


def match_annotations(
    annotations: pd.DataFrame, native_id: str, sequence: str, charge
) -> pd.DataFrame:
    """Select the exact spectrum/RNA/signed-charge triple, never a nearby scan.

    Repeated native records are safe only when their complete annotation
    multisets agree. Conflicting records are ambiguous and raise ValueError.
    All annotations, including several assignments to one observed peak, stay
    in the returned table.
    """
    _require_columns(annotations, set(ANNOTATION_COLUMNS[:7]), "RNA annotations")
    native_id = native_spectrum_id(native_id)
    charge = _integer(charge, "Selected precursor charge")
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("Select an RNA sequence label before matching annotations.")
    selected = annotations.loc[
        (annotations["native_id"] == native_id)
        & (annotations["sequence"] == sequence)
        & (annotations["charge"] == charge)
    ].copy()
    provenance = ["identification_index", "hit_index"]
    if set(provenance) <= set(selected.columns):
        signatures = {
            tuple(
                sorted(
                    zip(
                        group["mz"],
                        group["intensity"],
                        group["ion"],
                        group["fragment_charge"],
                    )
                )
            )
            for _, group in selected.groupby(provenance, dropna=False)
        }
        if len(signatures) > 1:
            raise ValueError(
                "Conflicting native annotation hits match this spectrum, RNA sequence "
                "and charge; an unambiguous fragment assignment cannot be selected."
            )
    return selected


def annotate_peaks(
    peaks: pd.DataFrame,
    annotations: pd.DataFrame,
    tolerance_da: float = 1e-5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Label uniquely matching viewed peaks and return unplaced assignments.

    The absolute tolerance accommodates file serialization, not a new fragment
    search. A native assignment matching more than one viewed peak is left
    unplaced. Multiple native ion labels on one unique peak are joined with
    ``; ``. Original rows, ordering, indices and intensities are preserved.
    """
    _require_columns(peaks, {"mz", "intensity"}, "spectrum peaks")
    _require_columns(annotations, {"mz", "ion"}, "RNA annotations")
    tolerance = _finite_number(tolerance_da, "Annotation tolerance")
    if tolerance < 0:
        raise ValueError("Annotation tolerance must be zero or positive in Da.")
    mz_values = np.array(
        [_finite_number(value, "Spectrum m/z", positive=True) for value in peaks["mz"]]
    )
    for value in peaks["intensity"]:
        _finite_number(value, "Spectrum intensity")
    annotated = peaks.copy(deep=True)
    labels = [[] for _ in mz_values]
    order = np.argsort(mz_values, kind="stable")
    sorted_mz = mz_values[order]
    unplaced, statuses = [], []
    for position, (_, annotation) in enumerate(annotations.iterrows()):
        mz = _finite_number(annotation["mz"], "Annotated m/z", positive=True)
        ion = annotation["ion"]
        if not isinstance(ion, str) or not ion.strip():
            raise ValueError("Each RNA fragment assignment must have an ion label.")
        left = np.searchsorted(sorted_mz, mz - tolerance, side="left")
        right = np.searchsorted(sorted_mz, mz + tolerance, side="right")
        if right - left == 1:
            labels[order[left]].append(ion)
        else:
            unplaced.append(position)
            statuses.append(
                "ambiguous peak match" if right > left else "no matching peak"
            )
    annotated["ion"] = ["; ".join(dict.fromkeys(ions)) for ions in labels]
    unmatched = annotations.iloc[unplaced].copy()
    unmatched["match_status"] = statuses
    return annotated, unmatched
