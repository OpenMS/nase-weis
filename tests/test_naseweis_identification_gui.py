"""Exercise native RNA assignment selection through the real results page."""

# Pytest injects the shared app fixture through the exported name.
# ruff: noqa: F811

import pytest
import streamlit as st
from test_naseweis_annotations import store_annotations
from test_naseweis_viewers_gui import (
    NATIVE_ID,
    array,
    click,
    peak_coordinates,
    plots,
    saved_result,
    widget,
)

from test_gui import _results_page, app  # noqa: F401

ASSIGNMENTS = [
    (100.0, 0.5, "d1", -1),
    (100.0, 0.5, "w1", -2),
    (200.0, 1.0, "a2-B", -1),
    (999.0, 0.2, "c4", -1),
]


def _identified_result(app, tmp_path, *, sequence="ACGU", charge=-2, assignments=None):
    results, original, snapshot = saved_result(app, tmp_path)
    native = store_annotations(tmp_path, [{
        "native_id": NATIVE_ID,
        "hits": [{
            "sequence": sequence, "charge": charge,
            "annotations": ASSIGNMENTS if assignments is None else assignments,
        }],
    }])
    destination = results / "idxml_results/search.idXML"
    destination.write_bytes(native.read_bytes())
    sources = [original, snapshot, destination, results / "mztab_results/search.mzTab"]
    _results_page(app)
    click(app, "Open spectrum viewer")
    return results, {path: path.read_bytes() for path in sources}


def _assignment_table(app):
    return next(table.value for table in app.dataframe if "observed_mz" in table.value)


def _labels(app):
    return {label["text"] for label in plots(app)[0]["layout"].get("annotations", [])}


def _selected_trace(app):
    return next((trace for trace in plots(app)[0]["data"] if trace.get("name") == "Selected ion"), None)


def test_visible_ion_table_reports_original_signal_and_unplaced_native_assignments(app, tmp_path):
    _identified_result(app, tmp_path)
    assert any(item.value == "RNA identification viewer" for item in app.subheader)
    table = _assignment_table(app).set_index("ion")
    assert set(table.index) == {"d1", "w1", "a2-B", "c4"}
    assert len(table) == 4
    assert table.loc["d1", "observed_mz"] == 100.0
    assert table.loc["d1", "observed_intensity"] == 1000.0
    assert table.loc["w1", "observed_intensity"] == 1000.0
    assert table.loc["a2-B", "observed_intensity"] == 2000.0
    assert table.loc["d1", "relative_intensity"] == 50.0
    assert table.loc["a2-B", "relative_intensity"] == 100.0
    assert table.loc["c4", "match_status"] == "no matching peak"
    assert table.loc["d1", "intensity"] == 0.5
    assert any("4 assignments · 2 located raw peaks · 1 unplaced" in item.value for item in app.caption)
    assert _selected_trace(app) is None
    assert any(item.proto.label == "Download fragment assignments (TSV)" for item in app.get("download_button"))


@pytest.mark.parametrize("series,ion,mass", [("d", "d1", 100.0), ("w", "w1", 100.0), ("a-B", "a2-B", 200.0)])
def test_ion_series_filter_changes_labels_without_changing_acquired_signal_or_saved_outputs(
    app, tmp_path, series, ion, mass
):
    _, sources = _identified_result(app, tmp_path)
    widget(app.multiselect, "Ion series").set_value([series]).run(timeout=30)
    assert not app.exception
    table = _assignment_table(app)
    assert table["ion"].tolist() == [ion]
    assert table["observed_mz"].tolist() == [mass]
    assert _labels(app) == {ion}
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]
    for path, content in sources.items():
        assert path.read_bytes() == content


def test_empty_ion_filter_keeps_raw_scan_and_can_be_corrected(app, tmp_path):
    _identified_result(app, tmp_path)
    widget(app.multiselect, "Ion series").set_value([]).run(timeout=30)
    assert not app.exception
    assert any("No assigned ions in the selected series" in item.value for item in app.info)
    assert not _labels(app)
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]
    widget(app.multiselect, "Ion series").set_value(["a-B"]).run(timeout=30)
    assert not app.exception
    assert _assignment_table(app)["ion"].tolist() == ["a2-B"]
    assert _labels(app) == {"a2-B"}


@pytest.mark.parametrize("query", ["acgu", "CONTROLLERTYPE=0", "scan=11"])
def test_identification_search_finds_rna_and_complete_native_reference(app, tmp_path, query):
    _identified_result(app, tmp_path)
    widget(app.text_input, "Find identification").set_value(query).run(timeout=30)
    assert not app.exception
    assert len(widget(app.selectbox, "Spectrum match").options) == 1
    assert NATIVE_ID in widget(app.selectbox, "Spectrum match").options[0]
    assert set(_assignment_table(app)["ion"]) == {"d1", "w1", "a2-B", "c4"}


def test_identification_search_has_recoverable_no_results_state_and_treats_query_literally(app, tmp_path):
    _identified_result(app, tmp_path)
    widget(app.text_input, "Find identification").set_value(".*").run(timeout=30)
    assert not app.exception
    assert any("No identifications match this filter" in item.value for item in app.info)
    assert not plots(app)
    widget(app.text_input, "Find identification").set_value("").run(timeout=30)
    assert not app.exception
    assert set(_assignment_table(app)["ion"]) == {"d1", "w1", "a2-B", "c4"}
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]


@pytest.mark.parametrize("sequence,charge,assignments", [
    ("UGCA", -2, ASSIGNMENTS), ("ACGU", 2, ASSIGNMENTS), ("ACGU", -2, []),
])
def test_exact_hit_without_native_assignments_explains_unannotated_scan(
    app, tmp_path, sequence, charge, assignments
):
    _identified_result(app, tmp_path, sequence=sequence, charge=charge, assignments=assignments)
    assert any("No saved fragment assignments match this exact spectrum" in item.value for item in app.info)
    assert not any("observed_mz" in table.value for table in app.dataframe)
    assert not _labels(app)
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]


@pytest.fixture
def dataframe_selection_events(monkeypatch):
    """Supply only the frontend row event AppTest's Dataframe cannot send.

    Rendering, selected-table data, native readers and plotting remain real.
    Browser acceptance separately exercises the actual row-selection control.
    """
    render = st.dataframe
    events, rendered = {}, []

    def selected_dataframe(data, *args, **kwargs):
        result = render(data, *args, **kwargs)
        key = kwargs.get("key", "")
        if kwargs.get("selection_mode") == "multi-row" and "-ions-" in key:
            rendered.append(key)
            return {"selection": {"rows": events.get(key, [])}}
        return result

    monkeypatch.setattr(st, "dataframe", selected_dataframe)
    return events, rendered


def test_selected_assignments_highlight_only_their_exact_acquired_peak(app, tmp_path, dataframe_selection_events):
    events, rendered = dataframe_selection_events
    _identified_result(app, tmp_path)
    events[rendered[-1]] = [0, 1, 3]
    app.run(timeout=30)
    assert not app.exception
    selected = _selected_trace(app)
    assert selected is not None
    assert array(selected["x"]).tolist() == [100.0]
    assert array(selected["y"]).tolist() == [50.0]
    assert any("cannot be highlighted" in item.value for item in app.info)
    selected_ions = next(item.value.removeprefix("Selected ions: ") for item in app.caption if item.value.startswith("Selected ions: "))
    assert set(selected_ions.split(", ")) == {"d1", "w1", "c4"}
    assert set(_labels(app)) == {"d1; w1", "a2-B"} or set(_labels(app)) == {"w1; d1", "a2-B"}


def test_filter_changes_do_not_apply_old_selected_row_positions_to_other_ions(app, tmp_path, dataframe_selection_events):
    events, rendered = dataframe_selection_events
    _identified_result(app, tmp_path)
    old_key = rendered[-1]
    events[old_key] = [0]
    app.run(timeout=30)
    assert array(_selected_trace(app)["x"]).tolist() == [100.0]
    widget(app.multiselect, "Ion series").set_value(["a-B"]).run(timeout=30)
    assert not app.exception
    assert rendered[-1] != old_key
    assert _selected_trace(app) is None
    assert _assignment_table(app)["observed_mz"].tolist() == [200.0]

