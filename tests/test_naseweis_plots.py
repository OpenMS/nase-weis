"""The Insight plot specifications NASE-weis draws, and their scientific display contracts."""

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from src.naseweis_plots import (
    ION_COLORS,
    SELECTED_ION_COLOR,
    TRACE_COLORS,
    UNKNOWN_ION_COLOR,
    Plot,
    abundance_bars_plot,
    chromatogram_plot,
    feature_map_plot,
    feature_xic_plot,
    ion_series,
    peakmap_plots,
    spectrum_plot,
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


def test_spectrum_uses_exact_masses_relative_intensities_and_rna_colors(peaks):
    original = peaks.copy(deep=True)
    plot = spectrum_plot(peaks, "RNA match")

    assert isinstance(plot, Plot) and plot.component == "LinePlot"
    rows = plot.data.to_dicts()
    assert [(r["mz"], r["intensity"]) for r in rows] == [(100.25, 25), (200.5, 50), (300.75, 100), (400.125, 0)]
    assert [r["ion"] for r in rows] == ["a1-B", "d2-", "w3-", ""]
    assert [r["color"] for r in rows] == [ION_COLORS["a-B"], ION_COLORS["d"], ION_COLORS["w"], None]
    assert [r["matched"] for r in rows] == [True, True, True, False]
    assert plot.settings["color_column"] == "color"
    assert plot.settings["annotation_column"] == "ion"
    assert plot.settings["title"] == "RNA match"
    assert list(plot.table["relative_intensity"]) == [25, 50, 100, 0]
    assert_frame_equal(peaks, original)


def test_unannotated_spectrum_has_no_invented_ions(peaks):
    plot = spectrum_plot(peaks.drop(columns="ion"), "Raw MS2")
    assert set(plot.data["ion"]) == {""}
    assert not any(plot.data["matched"])
    assert set(plot.data["color"].to_list()) == {None}
    assert UNKNOWN_ION_COLOR not in ION_COLORS.values()


def test_selected_ions_emphasize_exact_acquired_peaks_without_changing_signal(peaks):
    before = peaks.copy(deep=True)
    plain = spectrum_plot(peaks, "RNA match")
    selected = spectrum_plot(peaks, "RNA match", highlighted_mz=[200.5, 300.75, 200.5])

    assert selected.data["intensity"].to_list() == plain.data["intensity"].to_list()
    assert selected.data["ion"].to_list() == plain.data["ion"].to_list()
    colors = selected.data["color"].to_list()
    assert colors[1] == colors[2] == SELECTED_ION_COLOR
    assert colors[0] == plain.data["color"][0]
    assert list(selected.table["selected"]) == [False, True, True, False]
    assert_frame_equal(peaks, before)


def test_zero_height_selected_peak_remains_visible_without_inventing_intensity(peaks):
    plot = spectrum_plot(peaks, "Zero height", highlighted_mz=[400.125])
    row = plot.data.row(3, named=True)
    assert row["intensity"] == 0
    assert row["matched"] is True
    assert row["color"] == SELECTED_ION_COLOR


def test_empty_selection_keeps_the_plot_unchanged(peaks):
    assert spectrum_plot(peaks, "RNA", highlighted_mz=[]).data.equals(spectrum_plot(peaks, "RNA").data)


@pytest.mark.parametrize("selected", [
    [float("nan")], [float("inf")], [0], [-1], [True], ["bad"],
    "100.25", 100.25, [[100.25]], [100.2500001], [999.0],
])
def test_invalid_or_unobserved_highlights_are_not_snapped_to_a_peak(peaks, selected):
    with pytest.raises(ValueError, match="Highlighted|highlighted"):
        spectrum_plot(peaks, "Invalid selection", highlighted_mz=selected)


def test_ambiguous_duplicate_mass_cannot_highlight_an_arbitrary_peak(peaks):
    duplicate = pd.concat([peaks, peaks.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="exactly one acquired peak"):
        spectrum_plot(duplicate, "Ambiguous", highlighted_mz=[100.25])


@pytest.mark.parametrize("label,series", [
    ("a1-B", "a-B"), ("a-B2-", "a-B"), ("d2-", "d"),
    ("w3^2-", "w"), ("y1-H2O", "y"), ("precursor", ""),
])
def test_ion_series_uses_native_rna_labels(label, series):
    assert ion_series(label) == series


def test_chromatogram_traces_overlay_without_replacing_or_rescaling_each_other():
    points = pd.DataFrame({
        "rt": [20.0, 10.0, 20.0, 10.0, 20.0, 10.0, 10.0],
        "intensity": [200.0, 100.0, 40.0, 20.0, 10.0, 0.0, 5.0],
        "trace": ["TIC", "TIC", "BPC", "BPC", "XIC", "XIC", "M+1"],
    })
    original = points.copy(deep=True)
    plot = chromatogram_plot(points, "Acquisition chromatograms")

    assert plot.component == "SeriesPlot"
    frame = plot.data.to_pandas()
    for name, intensities in [("TIC", [100, 200]), ("BPC", [20, 40]), ("XIC", [0, 10])]:
        trace = frame[frame["series"] == name]
        assert list(trace["rt"]) == [10, 20]
        assert list(trace["intensity"]) == intensities
        assert plot.settings["series_colors"][name] == TRACE_COLORS[name]
    assert len(set(plot.settings["series_colors"].values())) == 4
    assert plot.settings["mode"] == "lines"
    assert_frame_equal(points, original)


@pytest.mark.parametrize("marginals", [False, True])
def test_peakmap_keeps_prepared_points_and_adds_summed_marginals(points, marginals):
    original = points.copy(deep=True)
    plots = peakmap_plots(points, "MS1 map", marginals)
    heatmap = plots[0]

    assert heatmap.component == "Heatmap"
    assert set(heatmap.data.iter_rows()) == set(points.itertuples(index=False, name=None))
    assert heatmap.settings["intensity_column"] == "intensity"
    assert len(plots) == (3 if marginals else 1)
    if marginals:
        rt, mz = plots[1], plots[2]
        assert rt.data["x"].to_list() == [10, 20] and rt.data["y"].to_list() == [20, 70]
        assert mz.data["x"].to_list() == [500, 500.5] and mz.data["y"].to_list() == [30, 60]
    assert_frame_equal(points, original)


def test_feature_map_preserves_coincident_rows_zero_intensities_and_negative_charge():
    features = pd.DataFrame({
        "rt": [20.0, 20.0], "mz": [500.25, 500.25],
        "intensity": [0.0, 100000.5], "label": ["[Am]CGU", "ACGU"], "charge": [-2, -3],
    })
    original = features.copy(deep=True)
    plot = feature_map_plot(features, "LFQ features")

    assert plot.component == "Heatmap"
    assert plot.data["rt"].to_list() == [20, 20]
    assert plot.data["mz"].to_list() == [500.25, 500.25]
    assert plot.data["intensity"].to_list() == [0, 100000.5]
    assert plot.settings["min_points"] == 2
    assert "charge -2" in plot.data["label"][0] and "[Am]CGU" in plot.data["label"][0]
    assert list(plot.table["charge"]) == [-2, -3]
    assert_frame_equal(features, original)


def test_feature_xic_draws_the_part_inside_its_boundaries_apart():
    trace = pd.DataFrame({"rt": [30.0, 10.0, 20.0], "intensity": [5.0, 1.0, 9.0]})
    plot = feature_xic_plot(trace, "ACGU", (15.0, 25.0))
    frame = plot.data.to_pandas()
    assert list(frame[frame["series"] == "XIC"]["rt"]) == [10, 20, 30]
    assert list(frame[frame["series"] == "Feature"]["rt"]) == [20]
    assert list(plot.table["in_feature"]) == [False, True, False]
    assert plot.settings["x_label"] == "Retention time (s)"
    unbounded = feature_xic_plot(trace, "ACGU")
    assert set(unbounded.data["series"]) == {"XIC"}
    assert not unbounded.table["in_feature"].any()


def test_abundance_bars_keep_missing_values_missing():
    plot = abundance_bars_plot(["A", "B"], [3.0, None], "ACGU", "Abundance")
    assert plot.settings["mode"] == "bars"
    assert plot.data["abundance"].to_list() == [3.0, None]
    with pytest.raises(ValueError, match="one value per acquisition"):
        abundance_bars_plot(["A"], [1.0, 2.0], "ACGU", "Abundance")
    with pytest.raises(ValueError, match="one value per acquisition"):
        abundance_bars_plot([], [], "ACGU", "Abundance")


@pytest.mark.parametrize("make", [
    lambda frame: spectrum_plot(frame, "Empty"),
    lambda frame: chromatogram_plot(frame, "Empty"),
    lambda frame: peakmap_plots(frame, "Empty"),
    lambda frame: feature_map_plot(frame, "Empty"),
    lambda frame: feature_xic_plot(frame, "Empty"),
])
def test_empty_tables_raise_an_actionable_error(make):
    columns = ["mz", "intensity", "rt", "trace", "label", "charge"]
    with pytest.raises(ValueError, match="available to plot"):
        make(pd.DataFrame(columns=columns))


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
        spectrum_plot(data, "Invalid")


@pytest.mark.parametrize("kind", ["trace_missing", "trace_null", "label_missing", "charge_fraction"])
def test_invalid_plot_metadata_is_actionable(points, kind):
    with pytest.raises(ValueError):
        if kind == "trace_missing":
            chromatogram_plot(points, "Invalid")
        elif kind == "trace_null":
            chromatogram_plot(points.assign(trace=None), "Invalid")
        elif kind == "label_missing":
            feature_map_plot(points.assign(charge=-2), "Invalid")
        else:
            feature_map_plot(points.assign(charge=-2.5, label="RNA"), "Invalid")


def test_untrusted_titles_and_labels_are_escaped(peaks, points):
    malicious = '</script><script>alert("viewer")</script><img src=x onerror=alert(1)>$x$'
    peaks["ion"] = malicious
    plots = [
        spectrum_plot(peaks, malicious),
        chromatogram_plot(points.drop_duplicates("rt").assign(trace=malicious), malicious),
        peakmap_plots(points, malicious)[0],
        feature_map_plot(points.assign(label=malicious, charge=-2), malicious),
        feature_xic_plot(points, malicious),
        abundance_bars_plot(["A"], [1.0], malicious, "Abundance"),
    ]
    for plot in plots:
        shown = repr(plot.settings) + repr(plot.data.to_dicts())
        assert "<script>" not in shown
        assert "<img src=x" not in shown
        assert "&lt;/script&gt;" in shown
        assert "&#36;x&#36;" in shown
