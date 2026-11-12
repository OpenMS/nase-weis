"""Exercise workspace and saved-run viewers through the real Streamlit UI."""

# Pytest injects the imported app fixture through parameters named app.
# ruff: noqa: F811

import base64
import json
import os
from pathlib import Path

import numpy as np
import pytest
from test_naseweis_annotations import store_annotations

from test_gui import _finished_run, _results_page, app  # noqa: F401

NATIVE_ID = "controllerType=0 controllerNumber=1 scan=11"


def acquisition(
    path,
    *,
    second=False,
    empty=False,
    ms2_only=False,
    duplicate=False,
    empty_ms1=False,
    many_scans=False,
):
    """Write a tiny real mzML with exact scan identities and known raw signal."""
    import pyopenms as oms

    path.parent.mkdir(parents=True, exist_ok=True)
    experiment = oms.MSExperiment()
    definitions = (
        []
        if empty
        else [
            (1, 10.0, "scan=10", [499.999, 500.0, 500.001, 600.0], [10, 20, 30, 5]),
            (2, 10.5, NATIVE_ID, [100.0, 200.0], [1000, 2000]),
            (1, 11.0, "scan=12", [500.0, 600.0], [40, 9]),
            (2, 12.0, "scan=13", [], []),
        ]
    )
    if second:
        definitions = [
            (1, 20.0, "scan=90", [700.0], [900]),
            (2, 20.5, NATIVE_ID, [900.0, 950.0], [30, 10]),
        ]
    if duplicate:
        definitions.append((2, 13.0, NATIVE_ID, [300.0], [500]))
    if many_scans:
        definitions.extend(
            (2, 20.0 + number, f"scan={1000 + number}", [100.0], [100])
            for number in range(1001)
        )
    for level, rt, native_id, mz, intensity in definitions:
        if ms2_only and level != 2:
            continue
        if empty_ms1 and level == 1:
            mz, intensity = [], []
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(level)
        spectrum.setRT(rt)
        spectrum.setNativeID(native_id)
        spectrum.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
        spectrum.set_peaks(
            (np.array(mz, dtype=float), np.array(intensity, dtype=float))
        )
        if level == 2:
            precursor = oms.Precursor()
            precursor.setMZ(500.0)
            precursor.setCharge(-2)
            spectrum.setPrecursors([precursor])
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)
    return path


def widget(collection, label):
    return next(element for element in collection if element.label == label)


def click(app, label):
    widget(app.button, label).click().run(timeout=30)
    assert not app.exception


def plots(app):
    return [json.loads(element.proto.spec) for element in app.get("plotly_chart")]


def array(values):
    if isinstance(values, dict) and "bdata" in values:
        return np.frombuffer(
            base64.b64decode(values["bdata"]), dtype=np.dtype(values["dtype"])
        )
    return np.asarray(values)


def trace_values(figure, name):
    trace = next(trace for trace in figure["data"] if trace.get("name") == name)
    return array(trace["x"]), array(trace["y"])


def peak_coordinates(figure):
    """Only the nonzero top of each spectrum stick represents a peak."""
    points = []
    for trace in figure["data"]:
        for mz, intensity in zip(array(trace.get("x", [])), array(trace.get("y", []))):
            if mz is not None and intensity is not None and intensity > 0:
                points.append((float(mz), float(intensity)))
    return points


def raw_page(app, path=None):
    app.run(timeout=30)
    directory = Path(app.session_state.workspace) / "naseweis"
    directory.mkdir(exist_ok=True)
    if path is not None:
        (directory / "params.json").write_text(json.dumps({"mzML-files": str(path)}))
    app.switch_page("content/naseweis_raw_viewer.py").run(timeout=30)
    assert not app.exception
    return directory


def raw_input(app, **options):
    app.run(timeout=30)
    path = acquisition(
        Path(app.session_state.workspace) / "naseweis/input-files/mzML-files/tiny.mzML",
        **options,
    )
    raw_page(app, path)
    return path


def result_content(native_id=NATIVE_ID, sequence="ACGU"):
    return (
        "NUH\taccession\topt_sequence\nNUC\trna\tACGU\n"
        "OLH\tsequence\taccession\tstart\tend\n"
        f"OLI\t{sequence}\trna\t1\t4\n"
        "OSH\tsequence\tretention_time\tcharge\texp_mass_to_charge\t"
        "search_engine_score[1]\tsearch_engine_score[2]\tspectra_ref\n"
        f"OSM\t{sequence}\t10.5\t-2\t500.0\t90.5\t0.01\tms_run[1]:{native_id}\n"
    )


def saved_result(
    app, tmp_path, *, name="Saved acquisition", native_id=NATIVE_ID, **options
):
    results = _finished_run(app, name, content=result_content(native_id))
    directory = results.parents[2]
    original = acquisition(
        directory / "input-files/mzML-files" / f"{name}.mzML", **options
    )
    snapshot = results.parent / "params.json"
    snapshot.write_text(
        json.dumps(
            {
                "mzML-files": str(original),
                "add-decoys": True,
                "enable-lfq": True,
                "FeatureFinderMetaboIdent": {"extract:mz_window": 50.0},
            }
        )
    )
    (results / "idxml_results/search.idXML").unlink()
    # A different current acquisition has the same MS2 native ID. A viewer
    # reading the live Configure selection would therefore show the wrong peaks.
    current = acquisition(tmp_path / f"current-{name}.mzML", second=True)
    (directory / "params.json").write_text(json.dumps({"mzML-files": str(current)}))
    return results, original, snapshot


def feature_table(results, *, label="ACGU", alternate=False, boundaries=True):
    quantification = results / "quantification"
    quantification.mkdir(exist_ok=True)
    extras = "\tleftWidth\trightWidth" if boundaries else ""
    values = "\t9.9\t10.2" if boundaries else ""
    rows = (
        f"10.0\t500.0\t12345\t-2\t{label}\t10.0{values}\n"
        f"11.0\t600.0\t555\t-3\t{label}-second\t11.0{values}\n"
    )
    if alternate:
        rows = f"20.0\t700.0\t777\t-2\t{label}\t20.0{values}\n"
    (quantification / "features.tsv").write_text(
        f"#rt\tmz\tintensity\tcharge\tlabel\texpected_rt{extras}\n" + rows
    )


def test_raw_viewer_empty_workspace_links_to_upload(app):
    raw_page(app)
    assert any("Upload an acquisition" in message.value for message in app.info)
    assert not app.get("plotly_chart")
    assert any(element.proto.label == "Add data" for element in app.get("page_link"))


def test_raw_loading_shows_three_real_plots_and_offline_html_exports(app):
    raw_input(app)
    assert not plots(app)
    assert any("Load the selected acquisition" in item.value for item in app.caption)
    click(app, "Load acquisition")
    assert [tab.label for tab in app.tabs] == [
        "MS1 peak map",
        "Spectra",
        "Chromatograms",
    ]
    figures = plots(app)
    assert len(figures) == 3
    assert (
        sum(
            button.proto.label == "Download interactive plot (HTML)"
            for button in app.get("download_button")
        )
        == 3
    )
    assert peak_coordinates(figures[1])
    tic_rt, tic = trace_values(figures[2], "TIC")
    np.testing.assert_allclose(tic_rt, [10.0, 11.0])
    np.testing.assert_allclose(tic, [65.0, 49.0])
    widget(app.checkbox, "Show marginal signal").check().run(timeout=30)
    assert not app.exception
    assert len(plots(app)[0]["data"]) > 1


def test_raw_spectrum_filters_use_levels_literal_native_ids_and_empty_scans(app):
    raw_input(app)
    click(app, "Load acquisition")
    widget(app.selectbox, "MS level").select("2").run(timeout=30)
    spectra = widget(app.selectbox, "Spectrum")
    assert len(spectra.options) == 2
    assert peak_coordinates(plots(app)[1]) == [(100.0, 50.0), (200.0, 100.0)]
    widget(app.text_input, "Find spectrum by native ID").set_value("scan=13").run(
        timeout=30
    )
    assert any("This spectrum has no peaks" in item.value for item in app.info)
    assert len(plots(app)) == 2
    widget(app.text_input, "Find spectrum by native ID").set_value("not[a]regex").run(
        timeout=30
    )
    assert any("No spectra match" in item.value for item in app.info)
    assert not app.exception


def test_raw_file_switch_requires_a_fresh_load_and_never_reuses_old_peaks(app):
    original = raw_input(app)
    second = acquisition(original.parent / "second.mzML", second=True)
    app.run(timeout=30)
    click(app, "Load acquisition")
    widget(app.selectbox, "Acquisition").select(str(second)).run(timeout=30)
    assert not plots(app)
    click(app, "Load acquisition")
    assert peak_coordinates(plots(app)[1]) == [(700.0, 100.0)]
    widget(app.selectbox, "Acquisition").select(str(original)).run(timeout=30)
    assert len(plots(app)) == 3
    assert all(mz < 700 for mz, _ in peak_coordinates(plots(app)[1]))


def test_raw_xic_uses_exact_peaks_and_can_replace_tic_bpc(app):
    raw_input(app)
    click(app, "Load acquisition")
    widget(app.multiselect, "Signal traces").set_value([]).run(timeout=30)
    assert any("Choose a signal trace" in item.value for item in app.info)
    widget(app.checkbox, "Add an extracted ion chromatogram").check().run(timeout=30)
    widget(app.selectbox, "Tolerance unit").select("Da").run(timeout=30)
    widget(app.number_input, "XIC m/z").set_value(500.0).run(timeout=30)
    widget(app.number_input, "Tolerance (± Da)").set_value(0.0001).run(timeout=30)
    assert not app.exception
    rt, values = trace_values(plots(app)[-1], "XIC")
    np.testing.assert_allclose(rt, [10, 11])
    np.testing.assert_allclose(values, [20, 40])


@pytest.mark.parametrize(
    "options, message",
    [
        ({"empty": True}, "contains no spectra"),
        ({"ms2_only": True}, "contains no MS1 scans"),
    ],
)
def test_raw_missing_spectrum_or_ms1_signal_is_explained(app, options, message):
    raw_input(app, **options)
    click(app, "Load acquisition")
    assert any(message in item.value for item in app.info)
    if options.get("ms2_only"):
        assert len(plots(app)) == 1  # MS2 spectrum remains usable.
    else:
        assert not plots(app)


def test_corrupt_raw_file_reports_a_read_error_without_crashing(app, tmp_path):
    path = tmp_path / "corrupt.mzML"
    path.write_text("not mzML")
    raw_page(app, path)
    click(app, "Load acquisition")
    assert any("Could not display this acquisition" in item.value for item in app.error)


def test_raw_vendor_file_explains_conversion_without_trying_to_read_it(app, tmp_path):
    path = tmp_path / "acquisition.raw"
    path.write_bytes(b"vendor input is handled by native search, not the viewer")
    raw_page(app, path)
    assert any("viewer requires mzML" in item.value for item in app.info)
    assert not any(button.label == "Load acquisition" for button in app.button)


def test_raw_ms1_scans_without_peaks_keep_scan_and_signal_metadata_available(app):
    raw_input(app, empty_ms1=True)
    click(app, "Load acquisition")
    assert any("No MS1 peaks to display" in item.value for item in app.info)
    assert any("This spectrum has no peaks" in item.value for item in app.info)
    assert len(plots(app)) == 1
    _, tic = trace_values(plots(app)[0], "TIC")
    np.testing.assert_allclose(tic, [0, 0])


def test_raw_large_scan_list_explains_cap_and_filter_can_reach_a_later_scan(app):
    raw_input(app, many_scans=True)
    click(app, "Load acquisition")
    assert any("Showing the first 1,000 matches" in item.value for item in app.caption)
    assert len(widget(app.selectbox, "Spectrum").options) == 1000
    widget(app.text_input, "Find spectrum by native ID").set_value("scan=2000").run(
        timeout=30
    )
    assert not app.exception
    assert len(widget(app.selectbox, "Spectrum").options) == 1
    assert peak_coordinates(plots(app)[1]) == [(100.0, 100.0)]


def test_raw_discovery_includes_registered_and_legacy_inputs_but_not_housekeeping(
    app, tmp_path
):
    directory = raw_page(app)
    legacy = Path(app.session_state.workspace) / "mzML-files"
    legacy.mkdir(exist_ok=True)
    uploaded = acquisition(legacy / "legacy.mzML")
    external = acquisition(tmp_path / "registered.mzML", second=True)
    note = legacy / "notes.txt"
    note.write_text("housekeeping")
    (legacy / "external_files.txt").write_text(
        f"{external}\n\n{tmp_path / 'gone.mzML'}\n{note}\n"
    )
    (directory / "params.json").write_text(
        json.dumps({"mzML-files": str(tmp_path / "missing.mzML")})
    )
    app.run(timeout=30)
    assert not app.exception
    assert widget(app.selectbox, "Acquisition").options == [
        f"{uploaded.name} · {uploaded.parent.name}",
        f"{external.name} · {external.parent.name}",
    ]


def test_results_spectrum_uses_its_saved_acquisition_and_exact_native_id(app, tmp_path):
    _, original, _ = saved_result(app, tmp_path)
    _results_page(app)
    click(app, "Open spectrum viewer")
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]
    assert any(
        original.name in item.value and NATIVE_ID in item.value for item in app.caption
    )
    assert any("No idXML annotations were saved" in item.value for item in app.caption)


@pytest.mark.parametrize("missing, duplicate", [(True, False), (False, True)])
def test_result_missing_or_ambiguous_native_ids_never_substitute_nearby_scans(
    app, tmp_path, missing, duplicate
):
    saved_result(
        app,
        tmp_path,
        native_id="scan=11" if missing else NATIVE_ID,
        duplicate=duplicate,
    )
    _results_page(app)
    click(app, "Open spectrum viewer")
    assert any(
        "does not identify exactly one scan" in item.value for item in app.warning
    )
    assert not plots(app)


@pytest.mark.parametrize("unmatched, legacy", [(True, False), (False, True)])
def test_saved_native_rna_fragment_assignments_label_the_acquired_peaks(
    app, tmp_path, unmatched, legacy
):
    results, _, _ = saved_result(app, tmp_path)
    assignments = [
        (100.0, 0.5, "d1", -1),
        (100.0, 0.5, "w1", -2),
        (200.0, 1.0, "a2-B", -1),
    ]
    if unmatched:
        assignments.append((999.0, 0.2, "c4", -1))
    native = store_annotations(
        tmp_path,
        [
            {
                "native_id": NATIVE_ID,
                "hits": [{"sequence": "ACGU", "annotations": assignments}],
            }
        ],
    )
    (results / "idxml_results" / ("id_out" if legacy else "search.idXML")).write_bytes(
        native.read_bytes()
    )
    _results_page(app)
    click(app, "Open spectrum viewer")
    labels = {item["text"] for item in plots(app)[0]["layout"]["annotations"]}
    assert {"w1; d1", "a2-B"} <= labels or {"d1; w1", "a2-B"} <= labels
    assert (
        any(
            "Some scored fragments are absent or ambiguous" in item.value
            for item in app.info
        )
        is unmatched
    )
    assert any(
        button.proto.label == "Download fragment assignments (TSV)"
        for button in app.get("download_button")
    )
    table = next(
        table.value for table in app.dataframe if "fragment_charge" in table.value
    )
    assert len(table) == len(assignments)
    assert set(table["ion"]) == (
        {"d1", "w1", "a2-B", "c4"} if unmatched else {"d1", "w1", "a2-B"}
    )
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]


def test_corrupt_native_annotations_leave_the_acquired_spectrum_available(
    app, tmp_path
):
    results, _, _ = saved_result(app, tmp_path)
    (results / "idxml_results/search.idXML").write_text("broken idXML")
    _results_page(app)
    click(app, "Open spectrum viewer")
    assert any(
        "Native fragment annotations are unavailable" in item.value
        for item in app.warning
    )
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]


@pytest.mark.parametrize("missing", ["acquisition", "snapshot"])
def test_missing_saved_acquisition_or_snapshot_explains_spectrum_unavailability(
    app, tmp_path, missing
):
    _, original, snapshot = saved_result(app, tmp_path)
    (original if missing == "acquisition" else snapshot).unlink()
    _results_page(app)
    assert any("Spectrum viewing is unavailable" in item.value for item in app.info)
    assert not any(button.label == "Open spectrum viewer" for button in app.button)
    assert not plots(app)


@pytest.mark.parametrize(
    "snapshot_content, message",
    [
        ("[]", "snapshot is invalid"),
        ("{}", "no saved acquisition path"),
        ("broken JSON", "snapshot is unavailable"),
    ],
)
def test_invalid_saved_parameters_never_fall_back_to_live_acquisition(
    app, tmp_path, snapshot_content, message
):
    _, _, snapshot = saved_result(app, tmp_path)
    snapshot.write_text(snapshot_content)
    _results_page(app)
    assert any(message in item.value for item in app.info)
    assert not plots(app)


def test_result_without_native_references_explains_why_viewer_is_unavailable(app):
    _finished_run(app)
    _results_page(app)
    assert any(
        "No spectrum references are available" in item.value for item in app.info
    )
    assert not any(button.label == "Open spectrum viewer" for button in app.button)


def test_result_empty_scan_is_explained_without_plotting_an_empty_spectrum(
    app, tmp_path
):
    saved_result(app, tmp_path, native_id="scan=13")
    _results_page(app)
    click(app, "Open spectrum viewer")
    assert any("The acquired scan contains no peaks" in item.value for item in app.info)
    assert not plots(app)


def test_result_corrupt_saved_acquisition_is_a_viewer_error_not_a_page_exception(
    app, tmp_path
):
    _, original, _ = saved_result(app, tmp_path)
    original.write_text("broken mzML")
    _results_page(app)
    click(app, "Open spectrum viewer")
    assert any("Could not display this spectrum" in item.value for item in app.error)
    assert not plots(app)


def test_lfq_feature_selection_shows_native_areas_and_exact_original_xics(
    app, tmp_path
):
    results, _, _ = saved_result(app, tmp_path)
    feature_table(results)
    _results_page(app)
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 2
    _, values = trace_values(plots(app)[1], "XIC")
    np.testing.assert_allclose(
        values, [60, 40]
    )  # 25 ppm half-window, including neighboring isotope-free peaks.
    assert len(plots(app)[1]["layout"]["shapes"]) == 2
    feature = widget(app.selectbox, "LFQ feature")
    feature.select(1).run(timeout=30)
    assert not app.exception
    assert widget(app.number_input, "XIC m/z").value == 600.0
    _, values = trace_values(plots(app)[1], "XIC")
    np.testing.assert_allclose(values, [5, 9])
    widget(app.checkbox, "Focus on the feature's elution window").uncheck().run(
        timeout=30
    )
    assert not app.exception
    assert len(plots(app)) == 2


def test_lfq_feature_map_remains_available_without_original_acquisition(app, tmp_path):
    results, original, _ = saved_result(app, tmp_path)
    feature_table(results)
    original.unlink()
    _results_page(app)
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 1
    assert any(
        "feature map is available, but its acquisition cannot be opened" in item.value
        for item in app.info
    )
    assert not any(widget.label == "XIC m/z" for widget in app.number_input)


def test_lfq_switching_saved_runs_never_reuses_another_run_feature_or_trace(
    app, tmp_path
):
    first, _, _ = saved_result(app, tmp_path, name="First")
    feature_table(first, label="FIRST")
    second, _, _ = saved_result(app, tmp_path, name="Second", second=True)
    feature_table(second, label="SECOND", alternate=True, boundaries=False)
    _results_page(app)
    click(app, "Open LFQ feature viewer")
    assert widget(app.number_input, "XIC m/z").value == 700.0
    widget(app.selectbox, "Run").select(first.parent.name).run(timeout=30)
    assert not app.exception
    assert not plots(app)
    click(app, "Open LFQ feature viewer")
    assert widget(app.number_input, "XIC m/z").value == 500.0
    widget(app.selectbox, "LFQ feature").select(1).run(timeout=30)
    assert widget(app.number_input, "XIC m/z").value == 600.0
    widget(app.selectbox, "Run").select(second.parent.name).run(timeout=30)
    assert not app.exception
    assert widget(app.number_input, "XIC m/z").value == 700.0
    _, values = trace_values(plots(app)[1], "XIC")
    np.testing.assert_allclose(values, [900])


def test_invalid_saved_lfq_window_keeps_the_feature_map_and_explains_failure(
    app, tmp_path
):
    results, _, snapshot = saved_result(app, tmp_path)
    feature_table(results)
    params = json.loads(snapshot.read_text())
    params["FeatureFinderMetaboIdent"]["extract:mz_window"] = -1
    snapshot.write_text(json.dumps(params))
    _results_page(app)
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 1
    assert any("no valid saved extraction window" in item.value for item in app.warning)


@pytest.mark.parametrize(
    "settings",
    [[], None, "bad", {"extract:mz_window": None}, {"extract:mz_window": "bad"}],
)
def test_malformed_nested_lfq_settings_preserve_feature_map_and_show_a_warning(
    app, tmp_path, settings
):
    results, _, snapshot = saved_result(app, tmp_path)
    feature_table(results)
    params = json.loads(snapshot.read_text())
    params["FeatureFinderMetaboIdent"] = settings
    snapshot.write_text(json.dumps(params))
    _results_page(app)
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 1
    assert not any(item.label == "XIC m/z" for item in app.number_input)
    assert any("LFQ viewer could not display" in item.value for item in app.warning)


def save_acquisition_identity(snapshot, path):
    from src.naseweis_provenance import acquisition_fingerprint

    params = json.loads(snapshot.read_text())
    params["_acquisition_fingerprint"] = acquisition_fingerprint(path)
    snapshot.write_text(json.dumps(params))


def unverified_identity_warnings(app):
    return [
        warning.value
        for warning in app.warning
        if "unverified" in warning.value.lower()
        or (
            "identity" in warning.value.lower()
            and any(
                term in warning.value.lower()
                for term in ["unavailable", "not recorded"]
            )
        )
    ]


def test_recorded_acquisition_identity_allows_matched_spectra_and_lfq_xic(
    app, tmp_path
):
    results, original, snapshot = saved_result(app, tmp_path)
    feature_table(results)
    save_acquisition_identity(snapshot, original)
    _results_page(app)
    click(app, "Open spectrum viewer")
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 3
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]
    _, values = trace_values(plots(app)[2], "XIC")
    np.testing.assert_allclose(values, [60, 40])
    assert not unverified_identity_warnings(app)


def test_replaced_acquisition_is_blocked_and_same_bytes_restored_with_new_mtime_work(
    app, tmp_path
):
    results, original, snapshot = saved_result(app, tmp_path)
    feature_table(results)
    saved_bytes, saved_stat = original.read_bytes(), original.stat()
    save_acquisition_identity(snapshot, original)
    _results_page(app)
    click(app, "Open spectrum viewer")
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 3

    # The replacement deliberately preserves the selected native MS2 ID while
    # changing the signal. An ID-only association would show the wrong data.
    acquisition(original, second=True)
    app.run(timeout=30)
    assert not app.exception
    assert len(plots(app)) == 1  # Saved feature map remains, raw spectrum/XIC do not.
    assert not any(button.label == "Open spectrum viewer" for button in app.button)
    assert not any(item.label == "XIC m/z" for item in app.number_input)
    assert not unverified_identity_warnings(app)
    assert any(
        any(term in item.value.lower() for term in ["fingerprint", "changed", "match"])
        for item in app.info
    )

    original.write_bytes(saved_bytes)
    os.utime(
        original,
        ns=(saved_stat.st_atime_ns, saved_stat.st_mtime_ns + 5_000_000_000),
    )
    app.run(timeout=30)
    assert not app.exception
    click(app, "Open spectrum viewer")
    assert len(plots(app)) == 3
    assert peak_coordinates(plots(app)[0]) == [(100.0, 50.0), (200.0, 100.0)]
    _, values = trace_values(plots(app)[2], "XIC")
    np.testing.assert_allclose(values, [60, 40])


def test_legacy_runs_allow_viewing_with_an_explicit_unverified_identity_warning(
    app, tmp_path
):
    results, _, _ = saved_result(app, tmp_path)
    feature_table(results)
    _results_page(app)
    click(app, "Open spectrum viewer")
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 3
    assert unverified_identity_warnings(app)


@pytest.mark.parametrize("record", [None, [], {}, "bad", {"sha256": "broken"}])
def test_malformed_identity_record_is_an_error_and_never_a_legacy_bypass(
    app, tmp_path, record
):
    results, _, snapshot = saved_result(app, tmp_path)
    feature_table(results)
    params = json.loads(snapshot.read_text())
    params["_acquisition_fingerprint"] = record
    snapshot.write_text(json.dumps(params))
    _results_page(app)
    assert not any(button.label == "Open spectrum viewer" for button in app.button)
    click(app, "Open LFQ feature viewer")
    assert len(plots(app)) == 1
    assert not any(item.label == "XIC m/z" for item in app.number_input)
    assert not unverified_identity_warnings(app)
    assert any(
        any(
            term in item.value.lower() for term in ["fingerprint", "identity", "record"]
        )
        for item in app.info
    )


def test_a_reference_to_an_ms1_scan_is_not_displayed_as_a_matched_ms2_spectrum(
    app, tmp_path
):
    saved_result(app, tmp_path, native_id="scan=10")
    _results_page(app)
    click(app, "Open spectrum viewer")
    assert not plots(app)
    assert any("does not identify exactly one" in item.value for item in app.warning)
