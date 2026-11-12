"""Exercise the installed plotting backend and scientific display contracts."""

import json

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from src.naseweis_plots import (
    ION_COLORS,
    SELECTED_ION_COLOR,
    UNKNOWN_ION_COLOR,
    chromatogram_figure,
    feature_map_figure,
    ion_series,
    peakmap_figure,
    spectrum_figure,
)


@pytest.fixture
def peaks():
    return pd.DataFrame({
        "mz": [100.25, 200.5, 300.75, 400.125],
        "intensity": [10.0, 20.0, 40.0, 0.0],
        "ion": ["a1-B", "d2-", "w3-", ""],
    })


@pytest.fixture
def points():
    return pd.DataFrame({
        "rt": [10.0, 10.0, 20.0, 20.0],
        "mz": [500.0, 500.5, 500.0, 500.5],
        "intensity": [0.0, 20.0, 30.0, 40.0],
    })


@pytest.fixture
def bounds():
    return pd.DataFrame({
        "rt_start": [9.5], "rt_end": [20.5],
        "mz_start": [499.75], "mz_end": [500.75],
        "label": ["[m5C]ACGU"],
    })


def peak_vertices(figure):
    return {(float(x), float(y)) for trace in figure.data for x, y in zip(trace.x, trace.y) if y != 0}


def test_spectrum_uses_exact_masses_relative_intensities_and_rna_colors(peaks):
    original = peaks.copy(deep=True)
    figure = spectrum_figure(peaks, "RNA match")

    assert peak_vertices(figure) == {(100.25, 25), (200.5, 50), (300.75, 100)}
    assert {trace.line.color for trace in figure.data} == {
        ION_COLORS["a-B"], ION_COLORS["d"], ION_COLORS["w"], UNKNOWN_ION_COLOR,
    }
    labels = {annotation.text: annotation for annotation in figure.layout.annotations}
    assert set(labels) == {"a1-B", "d2-", "w3-"}
    for label, series in [("a1-B", "a-B"), ("d2-", "d"), ("w3-", "w")]:
        assert labels[label].font.color == ION_COLORS[series]
    assert labels["d2-"].x == 200.5
    assert labels["d2-"].y == 50
    assert_frame_equal(peaks, original)


def test_mirror_spectrum_keeps_separate_normalization_and_reference_annotations(peaks):
    reference = pd.DataFrame({"mz": [101.25, 201.5], "intensity": [100.0, 50.0], "ion": ["a-B1-", "w2-"]})
    before = reference.copy(deep=True)
    figure = spectrum_figure(peaks, "Mirror", reference)

    assert peak_vertices(figure) == {
        (100.25, 25), (200.5, 50), (300.75, 100), (101.25, -100), (201.5, -50),
    }
    labels = {annotation.text: annotation for annotation in figure.layout.annotations}
    assert labels["a-B1-"].y == -100
    assert labels["w2-"].y == -50
    assert labels["a-B1-"].font.color == ION_COLORS["a-B"]
    assert any(trace.name.endswith("(reference)") for trace in figure.data)
    assert_frame_equal(reference, before)


def test_unannotated_spectrum_has_no_invented_ions(peaks):
    figure = spectrum_figure(peaks.drop(columns="ion"), "Raw MS2")
    assert not figure.layout.annotations
    assert {trace.line.color for trace in figure.data} == {UNKNOWN_ION_COLOR}


def test_selected_ions_emphasize_exact_acquired_peaks_without_changing_signal(peaks):
    before = peaks.copy(deep=True)
    plain = spectrum_figure(peaks, "RNA match")
    selected = spectrum_figure(peaks, "RNA match", highlighted_mz=[200.5, 300.75, 200.5])
    overlay = selected.data[-1]

    assert overlay.name == "Selected ion"
    assert list(overlay.x) == [200.5, 300.75]
    assert list(overlay.y) == [50, 100]
    assert overlay.marker.symbol == "circle-open"
    assert overlay.marker.color == SELECTED_ION_COLOR
    assert json.loads(selected.to_json())["data"][:-1] == json.loads(plain.to_json())["data"]
    assert selected.layout.annotations == plain.layout.annotations
    assert selected.layout.shapes[:-2] == plain.layout.shapes
    assert [(shape.x0, shape.x1, shape.y0, shape.y1) for shape in selected.layout.shapes[-2:]] == [
        (200.5, 200.5, 0, 50), (300.75, 300.75, 0, 100)
    ]
    assert_frame_equal(peaks, before)


def test_selection_highlight_does_not_change_mirror_normalization(peaks):
    reference = pd.DataFrame({"mz": [101.25], "intensity": [1000.0], "ion": ["w1-"]})
    figure = spectrum_figure(peaks, "Mirror", reference, highlighted_mz=[100.25])
    assert list(figure.data[-1].y) == [25]
    assert any(-100 in trace.y for trace in figure.data[:-1])
    assert peak_vertices(figure) == peak_vertices(spectrum_figure(peaks, "Mirror", reference))


def test_zero_height_selected_peak_remains_visible_without_inventing_intensity(peaks):
    figure = spectrum_figure(peaks, "Zero height", highlighted_mz=[400.125])
    assert list(figure.data[-1].x) == [400.125]
    assert list(figure.data[-1].y) == [0]


def test_empty_selection_keeps_existing_figure_unchanged(peaks):
    assert spectrum_figure(peaks, "RNA", highlighted_mz=[]).to_json() == spectrum_figure(peaks, "RNA").to_json()


@pytest.mark.parametrize("selected", [
    [float("nan")], [float("inf")], [0], [-1], [True], ["bad"],
    "100.25", 100.25, [[100.25]], [100.2500001], [999.0],
])
def test_invalid_or_unobserved_highlights_are_not_snapped_to_a_peak(peaks, selected):
    with pytest.raises(ValueError, match="Highlighted|highlighted"):
        spectrum_figure(peaks, "Invalid selection", highlighted_mz=selected)


def test_ambiguous_duplicate_mass_cannot_highlight_an_arbitrary_peak(peaks):
    duplicate = pd.concat([peaks, peaks.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="exactly one acquired peak"):
        spectrum_figure(duplicate, "Ambiguous", highlighted_mz=[100.25])


@pytest.mark.parametrize("label,series", [
    ("a1-B", "a-B"), ("a-B2-", "a-B"), ("d2-", "d"),
    ("w3^2-", "w"), ("y1-H2O", "y"), ("precursor", ""),
])
def test_ion_series_uses_native_rna_labels(label, series):
    assert ion_series(label) == series


def test_offline_export_keeps_selected_ion_evidence_and_escapes_its_label(peaks):
    peaks.loc[1, "ion"] = '</script><img src=x onerror=alert(1)>'
    figure = spectrum_figure(peaks, "RNA", highlighted_mz=[200.5])
    html = figure.to_html(include_plotlyjs=True, full_html=True)
    assert "Selected ion" in html
    assert SELECTED_ION_COLOR in html
    assert '</script><img src=x' not in html
    assert '&lt;/script&gt;' in figure.data[-1].text[0]
    assert 'src="https://cdn.plot.ly' not in html


def test_chromatogram_traces_overlay_without_replacing_or_rescaling_each_other():
    points = pd.DataFrame({
        "rt": [20.0, 10.0, 20.0, 10.0, 20.0, 10.0],
        "intensity": [200.0, 100.0, 40.0, 20.0, 10.0, 0.0],
        "trace": ["TIC", "TIC", "BPC", "BPC", "XIC", "XIC"],
    })
    original = points.copy(deep=True)
    figure = chromatogram_figure(points, "Acquisition chromatograms")
    traces = {trace.name: trace for trace in figure.data}

    assert set(traces) == {"TIC", "BPC", "XIC"}
    for name, intensities in [("TIC", [100, 200]), ("BPC", [20, 40]), ("XIC", [0, 10])]:
        assert list(traces[name].x) == [10, 20]
        assert list(traces[name].y) == intensities
    assert len({trace.line.color for trace in figure.data}) == 3
    assert_frame_equal(points, original)


@pytest.mark.parametrize("marginals", [False, True])
def test_peakmap_keeps_prepared_points_and_exact_feature_boundaries(points, bounds, marginals):
    original_points, original_bounds = points.copy(deep=True), bounds.copy(deep=True)
    figure = peakmap_figure(points, "MS1 map", bounds, marginals)
    main = figure.data[0]

    actual = set(zip(main.x, main.y, main.marker.color))
    assert actual == set(points.itertuples(index=False, name=None))
    assert main.marker.colorscale[0] == (0.0, "#440154")
    assert main.marker.colorscale[-1] == (1.0, "#fde725")
    assert main.marker.showscale is True
    assert main.marker.colorbar.title.text == "Intensity"
    assert (main.marker.cmin, main.marker.cmax) == (0, 40)
    box = figure.layout.shapes[0]
    assert (box.x0, box.x1, box.y0, box.y1) == (9.5, 20.5, 499.75, 500.75)
    assert (box.xref, box.yref) == (main.xaxis or "x", main.yaxis or "y")
    assert figure.data[-1].name == "[m5C]ACGU"
    if marginals:
        chromatogram = next(trace for trace in figure.data if list(trace.x) == [10, 20])
        assert list(chromatogram.y) == [20, 70]
        assert all(trace.marker.showscale is not True for trace in figure.data[1:])
    assert_frame_equal(points, original_points)
    assert_frame_equal(bounds, original_bounds)


def test_peakmap_accepts_absent_or_empty_feature_boundaries(points):
    for boxes in [None, pd.DataFrame()]:
        figure = peakmap_figure(points, "Map", boxes)
        assert len(figure.data) == 1
        assert not figure.layout.shapes


def test_feature_labels_can_be_missing_without_index_dependent_failures(points, bounds):
    bounds.index = ["noninteger"]
    unnamed = peakmap_figure(points, "Map", bounds.drop(columns="label"))
    assert unnamed.data[-1].name == "Feature 1"
    bounds.loc["noninteger", "label"] = None
    assert peakmap_figure(points, "Map", bounds).data[-1].name == ""


@pytest.mark.parametrize("changes,reason", [
    ({"rt_start": 21}, "lower limits"),
    ({"mz_start": 501}, "lower limits"),
    ({"mz_start": 0}, "positive"),
    ({"rt_end": float("inf")}, "finite"),
])
def test_invalid_feature_boundaries_are_rejected(points, bounds, changes, reason):
    for column, value in changes.items():
        bounds[column] = value
    with pytest.raises(ValueError, match=reason):
        peakmap_figure(points, "Map", bounds)


def test_feature_map_preserves_coincident_rows_zero_intensities_and_negative_charge(points):
    features = pd.DataFrame({
        "rt": [20.0, 20.0], "mz": [500.25, 500.25],
        "intensity": [0.0, 100000.5], "label": ["[Am]CGU", "ACGU"], "charge": [-2, -3],
    })
    original = features.copy(deep=True)
    figure = feature_map_figure(features, "LFQ features")
    trace = figure.data[0]

    assert list(trace.x) == [20, 20]
    assert list(trace.y) == [500.25, 500.25]
    assert list(trace.marker.color) == [0, 100000.5]
    assert trace.marker.size > 0
    assert list(trace.customdata[:, 0]) == ["[Am]CGU", "ACGU"]
    assert list(trace.customdata[:, 1]) == [-2, -3]
    assert list(trace.customdata[:, 2]) == [0, 100000.5]
    assert_frame_equal(features, original)


def test_all_zero_feature_intensities_still_have_visible_markers(points):
    features = points.assign(label="RNA", charge=-2, intensity=0)
    figure = feature_map_figure(features, "Zero-area features")
    assert list(figure.data[0].marker.color) == [0, 0, 0, 0]
    assert figure.data[0].marker.size == 8


@pytest.mark.parametrize("figure", [spectrum_figure, chromatogram_figure, peakmap_figure, feature_map_figure])
def test_empty_tables_raise_an_actionable_error(figure):
    columns = ["mz", "intensity", "rt", "trace", "label", "charge"]
    with pytest.raises(ValueError, match="available to plot"):
        figure(pd.DataFrame(columns=columns), "Empty")


@pytest.mark.parametrize("data,reason", [
    (None, "must be a table"),
    (pd.DataFrame({"mz": [1]}), "missing columns"),
    (pd.DataFrame({"mz": ["wrong"], "intensity": [1]}), "finite numbers"),
    (pd.DataFrame({"mz": [1], "intensity": [float("nan")]}), "finite numbers"),
    (pd.DataFrame({"mz": [1], "intensity": [-1]}), "negative"),
    (pd.DataFrame({"mz": [0], "intensity": [1]}), "positive"),
    (pd.DataFrame({"mz": [1], "intensity": [0]}), "positive intensities"),
])
def test_invalid_spectrum_measurements_are_not_silently_changed(data, reason):
    with pytest.raises(ValueError, match=reason):
        spectrum_figure(data, "Invalid")


@pytest.mark.parametrize("kind", ["trace_missing", "trace_null", "label_missing", "charge_fraction", "bounds_wrong_type", "bounds_missing_column"])
def test_invalid_plot_metadata_is_actionable(points, bounds, kind):
    with pytest.raises(ValueError):
        if kind == "trace_missing":
            chromatogram_figure(points, "Invalid")
        elif kind == "trace_null":
            chromatogram_figure(points.assign(trace=None), "Invalid")
        elif kind == "label_missing":
            feature_map_figure(points.assign(charge=-2), "Invalid")
        elif kind == "charge_fraction":
            feature_map_figure(points.assign(charge=-2.5, label="RNA"), "Invalid")
        elif kind == "bounds_wrong_type":
            peakmap_figure(points, "Invalid", bounds={})
        else:
            peakmap_figure(points, "Invalid", bounds.drop(columns="mz_start"))


def test_untrusted_titles_and_annotations_cannot_break_out_of_exported_html(peaks, points, bounds):
    malicious = '</script><script>alert("viewer")</script><img src=x onerror=alert(1)>$x$'
    peaks["ion"] = malicious
    peaks["native_id"] = malicious
    peaks.index = [malicious] * len(peaks)
    points.index = [malicious] * len(points)
    bounds["label"] = malicious
    figures = [
        spectrum_figure(peaks, malicious, peaks),
        chromatogram_figure(points.drop_duplicates("rt").assign(trace=malicious), malicious),
        peakmap_figure(points, malicious, bounds),
        feature_map_figure(points.assign(label=malicious, charge=-2), malicious),
    ]
    for figure in figures:
        html = figure.to_html(include_plotlyjs=True, full_html=True)
        serialized = json.dumps(figure.to_plotly_json(), default=str)
        assert malicious not in html
        assert '</script><script>alert' not in html
        assert '<img src=x onerror' not in serialized
        assert '&lt;/script&gt;' in serialized
        assert 'src="https://cdn.plot.ly' not in html
        assert '&#36;x&#36;' in serialized
