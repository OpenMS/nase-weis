"""Prove the reference comparator catches plausible silent scientific changes."""

import pandas as pd
import pytest
from naseweis_reference import REFERENCE, assert_nme1_reference


@pytest.mark.parametrize(
    "change", ["lost-match", "charge", "score", "q-value", "position", "reference"]
)
def test_scientific_regressions_are_rejected(change):
    tables = {
        section: pd.read_csv(REFERENCE / f"nme1-fdr.{section.lower()}.tsv", sep="\t")
        for section in ["NUC", "OLI", "OSM"]
    }
    if change == "lost-match":
        tables["OSM"] = tables["OSM"].iloc[1:]
    elif change == "charge":
        tables["OSM"].loc[0, "charge"] += 1
    elif change == "score":
        tables["OSM"].loc[0, "search_engine_score[1]"] += 1
    elif change == "q-value":
        tables["OSM"].loc[0, "search_engine_score[2]"] += 0.01
    elif change == "position":
        tables["OLI"].loc[0, "start"] += 1
    else:
        tables["NUC"].loc[0, "opt_sequence"] = "ACGU"
    with pytest.raises(AssertionError):
        assert_nme1_reference(tables, use_decoys=True)


def test_reference_comparison_ignores_only_row_order_and_roundoff():
    tables = {
        section: pd.read_csv(
            REFERENCE / f"nme1-fdr.{section.lower()}.tsv", sep="\t"
        ).iloc[::-1]
        for section in ["NUC", "OLI", "OSM"]
    }
    tables["OSM"]["search_engine_score[1]"] *= 1 + 1e-12
    assert_nme1_reference(tables, use_decoys=True)
