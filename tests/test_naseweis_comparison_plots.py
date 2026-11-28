"""Condition scatter plots retain saved relative RNA measurements exactly."""

import json

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from src.naseweis_comparison_plots import condition_scatter, paired_condition_rows


@pytest.fixture
def contrasts():
    return pd.DataFrame({
        "numerator_mean": [200.0, 50.0, 400.0, np.nan, 0.0, 13.0],
        "denominator_mean": [100.0, 100.0, 50.0, 10.0, 7.0, np.nan],
        "numerator_n": [3, 2, 3, 0, 1, 1],
        "denominator_n": [3, 3, 2, 1, 1, 0],
        "ratio": [2.0, 0.5, 8.0, np.nan, np.nan, np.nan],
        "log2_fold_change": [1.0, -1.0, 3.0, np.nan, np.nan, np.nan],
    }, index=pd.Index(["ACGUp", "A[m5C]GUp", "pACGUc", "GGA", "CCG", "AGG"], name="label"))


def markers(plot, series=None):
    """The plotted points (not the y = x line), optionally of one series."""
    rows = [r for r in plot.data.to_dicts() if r["series"] != "Equal"]
    return [r for r in rows if series is None or r["series"] == series]


def line(plot):
    return [r for r in plot.data.to_dicts() if r["series"] == "Equal"]


def test_paired_rows_exclude_missing_and_zero_without_modifying_input(contrasts):
    before = contrasts.copy(deep=True)
    paired = paired_condition_rows(contrasts)
    assert paired.index.tolist() == ["ACGUp", "A[m5C]GUp", "pACGUc"]
    assert_frame_equal(paired, contrasts.iloc[:3])
    paired.iloc[0, 0] = 99
    assert_frame_equal(contrasts, before)


def test_scatter_preserves_means_direction_replicates_and_fold_changes(contrasts):
    before = contrasts.copy(deep=True)
    plot = condition_scatter(contrasts, "NCL1", "control")
    assert plot.component == "SeriesPlot"
    unmodified, modified = markers(plot, "Unmodified"), markers(plot, "Modified")
    # Log scale: the axes are log10 of the saved means, nothing else.
    assert [10 ** r["x"] for r in unmodified] == pytest.approx([100, 50])
    assert [10 ** r["y"] for r in unmodified] == pytest.approx([200, 400])
    assert [10 ** r["x"] for r in modified] == pytest.approx([100])
    assert [10 ** r["y"] for r in modified] == pytest.approx([50])
    hover = modified[0]["label"]
    for part in ("A[m5C]GUp", "control mean 100", "(n=3)", "NCL1 mean 50", "(n=2)",
                 "ratio 0.5", "log2 FC -1"):
        assert part in hover
    names = plot.settings["series_names"]
    assert names["Unmodified"] == "Unmodified" and names["Modified"] == "Modified / ambiguous"
    assert plot.settings["x_label"].startswith("log10 mean control ·")
    assert plot.settings["y_label"].startswith("log10 mean NCL1 ·")
    assert "total-area normalized" in plot.settings["x_label"]
    table = plot.table.set_index("rna")
    assert table.loc["A[m5C]GUp", "numerator_mean"] == 50
    assert table.loc["A[m5C]GUp", "denominator_mean"] == 100
    assert bool(table.loc["A[m5C]GUp", "modified"]) is True
    assert_frame_equal(contrasts, before)


def test_equal_abundance_line_spans_both_axes(contrasts):
    plot = condition_scatter(contrasts, "NCL1", "control")
    points = line(plot)
    assert [r["x"] for r in points] == [r["y"] for r in points]
    assert [10 ** r["x"] for r in points] == pytest.approx([50.0, 400.0])
    assert plot.settings["series_names"]["Equal"] == "Equal abundance (y = x)"
    assert plot.settings["series_modes"] == {"Equal": "lines"}


@pytest.mark.parametrize("measure,label", [("intensity", "Native feature area"), ("model_area", "Fitted model area")])
def test_raw_linear_plot_uses_native_area_units(contrasts, measure, label):
    plot = condition_scatter(contrasts, "NCL1", "control", measure=measure, normalization="none", scale="linear")
    assert plot.settings["x_label"] == "Mean control · " + label.lower()
    assert plot.settings["y_label"] == "Mean NCL1 · " + label.lower()
    assert [r["x"] for r in markers(plot)] == [100, 100, 50]
    assert [r["x"] for r in line(plot)] == [50, 400]


def test_highlighting_selects_exact_identity_and_retains_other_rnas(contrasts):
    plot = condition_scatter(contrasts, "NCL1", "control", labels=["A[m5C]GUp", "CCG"])
    selected = markers(plot, "Selected")
    assert len(selected) == 1 and "A[m5C]GUp" in selected[0]["label"]
    assert len(markers(plot)) == 3
    assert list(plot.table["labelled"]) == [False, True, False]
    assert list(plot.table["modified"]) == [False, True, False]


@pytest.mark.parametrize("sequence,expected", [
    ("pACGUc", "Unmodified"), ("ACGUp", "Unmodified"),
    ("[5'-p]ACGU[3'-p]", "Unmodified"),
    ("ACGU[3'-c]", "Unmodified"), ("A[C|m5C]GU", "Modified"),
    ("A[m5C?]GU", "Modified"),
    ("AYGp", "Modified"), ("AIGp", "Modified"),
    ("ADGp", "Modified"), ("*ACGp", "Unmodified"),
    ("[5'-p*]ACGp", "Unmodified"),
])
def test_terminal_groups_are_distinct_from_nucleotide_annotations(sequence, expected):
    frame = pd.DataFrame({"numerator_mean": [2], "denominator_mean": [1]}, index=[sequence])
    assert markers(condition_scatter(frame, "treatment", "control"))[0]["series"] == expected


def test_imported_text_is_escaped_in_hover_axes_and_selected_labels():
    label = "A[<script>$x</script>]G"
    frame = pd.DataFrame({"numerator_mean": [2], "denominator_mean": [1]}, index=[label])
    plot = condition_scatter(frame, "NCL1<script>", "control$x", labels=[label])
    hover = markers(plot)[0]["label"]
    assert "<script>" not in hover
    assert "&lt;script&gt;" in hover
    assert "&#36;" in hover
    assert "control&#36;x" in plot.settings["x_label"]
    assert "NCL1&lt;script&gt;" in hover
    assert "<script>" not in repr(plot.settings)


@pytest.mark.parametrize("scale", ["log", "linear"])
@pytest.mark.parametrize("value", [np.nextafter(0.0, 1.0), 1.0, np.finfo(float).max])
def test_single_equal_point_is_finite_and_safe_to_export(value, scale):
    frame = pd.DataFrame({"numerator_mean": [value], "denominator_mean": [value]}, index=["ACGU"])
    plot = condition_scatter(frame, "NCL1", "control", scale=scale)
    values = [r[axis] for r in plot.data.to_dicts() for axis in ("x", "y")]
    assert all(np.isfinite(values))
    json.loads(json.dumps(plot.data.to_dicts()), parse_constant=lambda v: pytest.fail(v))


def test_nonfinite_metadata_remains_missing_in_hover(contrasts):
    contrasts.loc["ACGUp", "ratio"] = np.inf
    contrasts.loc["ACGUp", "log2_fold_change"] = np.nan
    hover = markers(condition_scatter(contrasts, "NCL1", "control"))[0]["label"]
    assert "ratio – · log2 FC –" in hover
    assert "nan" not in hover.lower() and "inf" not in hover.lower()


@pytest.mark.parametrize("changes", [
    {"numerator_mean": np.inf}, {"denominator_mean": -1},
    {"numerator_mean": 0}, {"denominator_mean": np.nan},
])
def test_unpaired_values_are_excluded_on_both_scales(changes):
    frame = pd.DataFrame({"numerator_mean": [2.0, 4.0], "denominator_mean": [1.0, 2.0]}, index=["ACGU", "GACU"])
    for field, value in changes.items():
        frame.loc["ACGU", field] = value
    assert paired_condition_rows(frame).index.tolist() == ["GACU"]
    for scale in ("log", "linear"):
        plot = condition_scatter(frame, "NCL1", "control", scale=scale)
        assert list(plot.table["rna"]) == ["GACU"]
        point = markers(plot)[0]
        expected = (2, 4) if scale == "linear" else (np.log10(2), np.log10(4))
        assert (point["x"], point["y"]) == pytest.approx(expected)


@pytest.mark.parametrize("frame,reason", [
    (None, "must be a table"), (pd.DataFrame(), "need numerator"),
    (pd.DataFrame({"numerator_mean": [1], "denominator_mean": ["bad"]}, index=["ACGU"]), "must be numeric"),
    (pd.DataFrame({"numerator_mean": [1], "denominator_mean": [1]}, index=[""]), "unique, nonempty"),
    (pd.DataFrame({"numerator_mean": [1, 2], "denominator_mean": [1, 2]}, index=["ACGU", "ACGU"]), "unique, nonempty"),
])
def test_invalid_comparison_input_fails_explicitly(frame, reason):
    with pytest.raises(ValueError, match=reason):
        paired_condition_rows(frame)


@pytest.mark.parametrize("options,reason", [
    ({"scale": "symlog"}, "scale must"), ({"measure": "other"}, "Unsupported"),
    ({"normalization": "other"}, "Unsupported"), ({"labels": "ACGU"}, "collection"),
    ({"labels": [None]}, "collection"),
])
def test_invalid_plot_settings_fail_explicitly(contrasts, options, reason):
    with pytest.raises(ValueError, match=reason):
        condition_scatter(contrasts, "NCL1", "control", **options)


@pytest.mark.parametrize("numerator,denominator", [("", "control"), (None, "control"), ("control", "control")])
def test_conditions_must_be_named_and_different(contrasts, numerator, denominator):
    with pytest.raises(ValueError, match="conditions"):
        condition_scatter(contrasts, numerator, denominator)


def test_no_positive_paired_measurements_are_not_fabricated(contrasts):
    with pytest.raises(ValueError, match="positive measured mean"):
        condition_scatter(contrasts.iloc[3:], "NCL1", "control")


def test_invalid_rna_identity_and_hover_statistics_fail_explicitly(contrasts):
    invalid = contrasts.rename(index={"ACGUp": "not RNA"})
    with pytest.raises(ValueError, match="Invalid RNA sequence"):
        condition_scatter(invalid, "NCL1", "control")
    contrasts["ratio"] = contrasts["ratio"].astype(object)
    contrasts.loc["ACGUp", "ratio"] = "bad"
    with pytest.raises(ValueError, match="ratio must be numeric"):
        condition_scatter(contrasts, "NCL1", "control")
