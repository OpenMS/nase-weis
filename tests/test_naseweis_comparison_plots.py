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


def markers(figure):
    return [trace for trace in figure.data if trace.mode != "lines"]


def test_paired_rows_exclude_missing_and_zero_without_modifying_input(contrasts):
    before = contrasts.copy(deep=True)
    paired = paired_condition_rows(contrasts)
    assert paired.index.tolist() == ["ACGUp", "A[m5C]GUp", "pACGUc"]
    assert_frame_equal(paired, contrasts.iloc[:3])
    paired.iloc[0, 0] = 99
    assert_frame_equal(contrasts, before)


def test_scatter_preserves_means_direction_replicates_and_fold_changes(contrasts):
    before = contrasts.copy(deep=True)
    figure = condition_scatter(contrasts, "NCL1", "control")
    traces = {trace.name: trace for trace in markers(figure)}
    assert set(traces) == {"Unmodified", "Modified / ambiguous"}
    assert list(traces["Unmodified"].x) == [100, 50]
    assert list(traces["Unmodified"].y) == [200, 400]
    assert list(traces["Modified / ambiguous"].x) == [100]
    assert list(traces["Modified / ambiguous"].y) == [50]
    assert traces["Modified / ambiguous"].customdata == (["A[m5C]GUp", 2, 3, 0.5, -1],)
    assert traces["Unmodified"].marker.symbol == "circle-open"
    assert traces["Modified / ambiguous"].marker.symbol == "circle"
    assert figure.layout.xaxis.title.text.startswith("control ·")
    assert figure.layout.yaxis.title.text.startswith("NCL1 ·")
    assert "total-area normalized" in figure.layout.xaxis.title.text
    assert "control mean: %{x" in traces["Unmodified"].hovertemplate
    assert "NCL1 mean: %{y" in traces["Unmodified"].hovertemplate
    assert_frame_equal(contrasts, before)


def test_equal_abundance_line_and_shared_log_axes(contrasts):
    figure = condition_scatter(contrasts, "NCL1", "control")
    line = figure.data[0]
    assert list(line.x) == list(line.y) == [50.0, 400.0]
    assert line.name == "Equal abundance (y = x)"
    assert figure.layout.xaxis.type == figure.layout.yaxis.type == "log"
    assert figure.layout.xaxis.range == figure.layout.yaxis.range
    assert figure.layout.yaxis.scaleanchor == "x"
    assert figure.layout.yaxis.scaleratio == 1


@pytest.mark.parametrize("measure,label", [("intensity", "Native feature area"), ("model_area", "Fitted model area")])
def test_raw_linear_plot_uses_native_area_units(contrasts, measure, label):
    figure = condition_scatter(contrasts, "NCL1", "control", measure=measure, normalization="none", scale="linear")
    assert figure.layout.xaxis.title.text == "control · " + label
    assert figure.layout.yaxis.title.text == "NCL1 · " + label
    assert figure.layout.xaxis.type == figure.layout.yaxis.type == "linear"
    assert list(figure.data[0].x) == list(figure.data[0].y) == [0, 400]
    assert figure.layout.xaxis.range[0] == 0


def test_highlighting_selects_exact_identity_and_retains_other_rnas(contrasts):
    figure = condition_scatter(contrasts, "NCL1", "control", labels=["A[m5C]GUp", "CCG"])
    unmodified, modified = markers(figure)
    assert unmodified.mode == "markers"
    assert modified.mode == "markers+text"
    assert list(modified.text) == ["A[m5C]GUp"]
    assert list(modified.marker.size) == [11]
    assert len(unmodified.x) + len(modified.x) == 3
    assert not figure.layout.annotations


@pytest.mark.parametrize("sequence,expected", [
    ("pACGUc", "Unmodified"), ("ACGUp", "Unmodified"),
    ("[5'-p]ACGU[3'-p]", "Unmodified"),
    ("ACGU[3'-c]", "Unmodified"), ("A[C|m5C]GU", "Modified / ambiguous"),
    ("A[m5C?]GU", "Modified / ambiguous"),
    ("AYGp", "Modified / ambiguous"), ("AIGp", "Modified / ambiguous"),
    ("ADGp", "Modified / ambiguous"), ("*ACGp", "Unmodified"),
    ("[5'-p*]ACGp", "Unmodified"),
])
def test_terminal_groups_are_distinct_from_nucleotide_annotations(sequence, expected):
    frame = pd.DataFrame({"numerator_mean": [2], "denominator_mean": [1]}, index=[sequence])
    assert markers(condition_scatter(frame, "treatment", "control"))[0].name == expected


def test_imported_text_is_escaped_in_hover_axes_and_selected_labels():
    label = "A[<script>$x</script>]G"
    frame = pd.DataFrame({"numerator_mean": [2], "denominator_mean": [1]}, index=[label])
    figure = condition_scatter(frame, "NCL1<script>", "control$x", labels=[label])
    trace = markers(figure)[0]
    assert "<script>" not in trace.text[0]
    assert "&lt;script&gt;" in trace.text[0]
    assert "&#36;" in trace.text[0]
    assert "control&#36;x" in figure.layout.xaxis.title.text
    assert "NCL1&lt;script&gt;" in trace.hovertemplate


@pytest.mark.parametrize("scale", ["log", "linear"])
@pytest.mark.parametrize("value", [np.nextafter(0.0, 1.0), 1.0, np.finfo(float).max])
def test_single_equal_point_has_a_finite_nonzero_diagonal_and_safe_export(value, scale):
    frame = pd.DataFrame({"numerator_mean": [value], "denominator_mean": [value]}, index=["ACGU"])
    figure = condition_scatter(frame, "NCL1", "control", scale=scale)
    diagonal = figure.data[0]
    assert diagonal.x[0] < diagonal.x[1]
    assert np.isfinite(diagonal.x).all()
    assert figure.layout.xaxis.range[0] < figure.layout.xaxis.range[1]
    assert all(np.isfinite(figure.layout.xaxis.range))
    # Plotly serializes optional missing hover fields as null, never NaN/Infinity.
    json.loads(figure.to_json(), parse_constant=lambda value: pytest.fail(value))


def test_nonfinite_metadata_remains_missing_in_hover(contrasts):
    contrasts.loc["ACGUp", "ratio"] = np.inf
    contrasts.loc["ACGUp", "log2_fold_change"] = np.nan
    trace = markers(condition_scatter(contrasts, "NCL1", "control"))[0]
    assert trace.customdata[0][3:] == [None, None]


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
        trace = markers(condition_scatter(frame, "NCL1", "control", scale=scale))[0]
        assert list(trace.x) == [2]
        assert list(trace.y) == [4]


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
