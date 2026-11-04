"""Compare scientific output, excluding runtime paths, timestamps and metadata."""

from pathlib import Path

import pandas as pd

REFERENCE = Path(__file__).parent / "data/naseweis"


def assert_nme1_reference(tables: dict, use_decoys: bool) -> None:
    name = "nme1-fdr" if use_decoys else "nme1-unfiltered"
    for section in ["NUC", "OLI", "OSM"]:
        expected = pd.read_csv(REFERENCE / f"{name}.{section.lower()}.tsv", sep="\t")
        columns = list(expected.columns)
        actual = tables[section][columns].copy()
        expected = expected.sort_values(columns).reset_index(drop=True)
        actual = actual.sort_values(columns).reset_index(drop=True)
        # Allow harmless floating-point rounding, but preserve all rows,
        # assignments, charges, sequences, coordinates, scores and q-values.
        pd.testing.assert_frame_equal(
            actual,
            expected,
            check_dtype=False,
            check_exact=False,
            rtol=1e-9,
            atol=1e-9,
            obj=f"{name} {section}",
        )
