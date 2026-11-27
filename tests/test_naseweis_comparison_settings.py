"""The settings of Edit comparison and the table of identified RNAs."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from src import naseweis_experiment_ui as ui
from src.naseweis_comparison import normalized
from src.naseweis_result_plots_ui import rna_summary


def test_normalized_none_is_a_copy():
    raw = pd.DataFrame({"a": [1.0, 3.0], "b": [2.0, 2.0]})
    out = normalized(raw, "none")
    assert out.equals(raw) and out is not raw


def test_normalized_total_divides_by_column_sums():
    raw = pd.DataFrame({"a": [1.0, 3.0], "b": [0.0, 0.0]})
    out = normalized(raw, "total")
    assert out["a"].tolist() == [0.25, 0.75]
    assert out["b"].isna().all()


def test_normalized_refuses_unknown_methods():
    with pytest.raises(ValueError):
        normalized(pd.DataFrame({"a": [1.0]}), "median-of-ratios")


@pytest.fixture
def state(monkeypatch):
    session = {}
    monkeypatch.setattr(ui, "st", SimpleNamespace(session_state=session))
    return session


def test_comparison_choices_default_to_last_over_first(state):
    assert ui.comparison_choices("p", ["ctrl", "drug"], "none") == ("drug", "ctrl", "none")


def test_comparison_choices_follow_the_dialog_and_reset(state):
    num, den, norm = ui._choice_keys("p")
    state.update({num: "ctrl", den: "drug", norm: "total"})
    assert ui.comparison_choices("p", ["ctrl", "drug"], "none") == ("ctrl", "drug", "total")
    # A pick that names no condition of this experiment falls back
    state[num] = "gone"
    assert ui.comparison_choices("p", ["ctrl", "drug"], "none")[0] == "drug"
    ui.reset_comparison_choices("p")
    assert not state


def test_rna_summary_without_rnas_is_empty():
    table = rna_summary([("all", {})])
    assert table.empty
    assert list(table.columns) == ["rna_idx", "accession", "length", "oligonucleotides",
                                   "matches", "coverage"]


def test_rna_summary_keeps_an_rna_whose_coverage_cannot_be_built():
    tables = {"NUC": pd.DataFrame({"accession": ["RNA1", np.nan, "RNA1", "RNA2"]})}
    table = rna_summary([("a", tables), ("b", tables)])
    assert table["accession"].tolist() == ["RNA1", "RNA2"]
    assert table["rna_idx"].tolist() == [0, 1]
