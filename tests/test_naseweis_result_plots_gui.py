"""Coverage and comparative plots remain tied to saved scientific evidence."""

import json

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from test_naseweis_experiment_gui import app_for, saved_experiment  # noqa: F401

from src.naseweis_experiment import sample_results
from src.naseweis_result_plots_ui import _coverage_tracks


def _coverage_file(root, sample_id, *, sequence="ACGUACGU", legacy=False):
    folder = sample_results(root, sample_id) / "mztab_results"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / ("tab_out" if legacy else "search.mzTab")
    path.write_text(
        "NUH\taccession\topt_sequence\n"
        f"NUC\tNME1\t{sequence}\n"
        "OLH\tsequence\taccession\tstart\tend\n"
        "OLI\tACGp\tNME1\t1\t3\n"
        "OLI\tUA[m5C]Gp\tNME1\t4\t7\n"
        "OSH\tsequence\tspectra_ref\n"
        "OSM\tACGp\tms_run[1]:scan=1\n"
        "OSM\tUA[m5C]Gp\tms_run[1]:scan=2\n"
    )
    return path


def _download(app, label):
    return next(item for item in app.get("download_button") if item.proto.label == label)


def _scatter(app):
    return next(
        json.loads(item.proto.spec)
        for item in app.get("plotly_chart")
        if "scatter" in item.proto.id
    )


def test_unassigned_acquisition_does_not_pool_with_a_matching_condition_name():
    available = [
        ({"file": "one.mzML", "condition": ""}, {"evidence": 1}),
        ({"file": "two.mzML", "condition": "Unassigned · one.mzML"}, {"evidence": 2}),
    ]
    tracks = _coverage_tracks(available, "Condition")
    assert tracks[0][0] != tracks[1][0]
    assert tracks[1][0] == "Unassigned · one.mzML"
    assert [track[1] for track in tracks] == [tables for _, tables in available]


def test_identification_only_experiment_has_condition_and_acquisition_coverage(saved_experiment):  # noqa: F811
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"]:
        _coverage_file(root, sample["id"], legacy=sample["id"] == "sample_004")
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = app_for(folder, mode="results").run(timeout=30)
    assert not app.exception
    assert app.selectbox[0].label == "RNA accession"
    assert app.selectbox[0].value == "NME1"
    assert next(radio for radio in app.radio if radio.label == "Coverage tracks").value == "Condition"
    _download(app, "Download coverage (HTML)")
    _download(app, "Download coverage (SVG)")
    assert not app.get("plotly_chart")
    assert any("4 of 4 acquisitions" in text.value for text in app.caption)
    next(radio for radio in app.radio if radio.label == "Coverage tracks").set_value("Acquisition").run(timeout=30)
    assert not app.exception
    app.select_slider[0].set_value(25).run(timeout=30)
    assert not app.exception
    assert app.session_state.opened_sample == "sample_001"


def test_partial_coverage_uses_only_available_saved_results_and_reports_bad_file(saved_experiment):  # noqa: F811
    folder, root, manifest = saved_experiment
    _coverage_file(root, "sample_001")
    invalid = _coverage_file(root, "sample_002")
    invalid.write_text("NUH\taccession\topt_sequence\nNUC\tbroken\n")
    manifest["status"] = "failed"
    manifest["samples"][2]["status"] = "failed"
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = app_for(folder, mode="results").run(timeout=30)
    assert not app.exception
    assert any("Could not read coverage for run2.mzML" in text.value for text in app.warning)
    assert any("1 of 4 acquisitions" in text.value for text in app.caption)
    _download(app, "Download coverage (SVG)")
    assert not any(radio.label == "Coverage tracks" for radio in app.radio)
    assert not app.get("plotly_chart")


def test_conflicting_references_do_not_share_an_axis_but_each_acquisition_remains_accessible(saved_experiment):  # noqa: F811
    folder, root, _manifest = saved_experiment
    _coverage_file(root, "sample_001")
    _coverage_file(root, "sample_002", sequence="UUUUACGU")
    app = app_for(folder, mode="results").run(timeout=30)
    assert not app.exception
    assert any("Conflicting reference sequences" in text.value for text in app.warning)
    assert not any(item.proto.label == "Download coverage (SVG)" for item in app.get("download_button"))
    assert app.session_state.opened_sample == "sample_001"
    assert _scatter(app)["layout"]["xaxis"]["type"] == "log"


def test_blank_conditions_default_to_acquisition_coverage(saved_experiment):  # noqa: F811
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"][:2]:
        _coverage_file(root, sample["id"])
        sample["condition"] = ""
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = app_for(folder, mode="results").run(timeout=30)
    assert not app.exception
    grouping = next(radio for radio in app.radio if radio.label == "Coverage tracks")
    assert grouping.value == "Acquisition"
    grouping.set_value("Condition").run(timeout=30)
    assert not app.exception
    _download(app, "Download coverage (SVG)")


def test_lfq_scatter_uses_saved_means_direction_labels_and_view_only_controls(saved_experiment):  # noqa: F811
    folder, root, _manifest = saved_experiment
    before = {path: path.read_bytes() for path in root.rglob("*") if path.is_file()}
    app = app_for(folder, mode="results").run(timeout=30)
    assert not app.exception
    spec = _scatter(app)
    assert spec["layout"]["xaxis"]["type"] == "log"
    assert "control" in spec["layout"]["xaxis"]["title"]["text"]
    assert "treated" in spec["layout"]["yaxis"]["title"]["text"]
    assert any("2 oligonucleotides measured" in text.value for text in app.caption)
    _download(app, "Download abundance comparison (HTML)")
    app.multiselect[0].set_value(["ACGU"]).run(timeout=30)
    next(radio for radio in app.radio if radio.label == "Scatter axes").set_value("Linear").run(timeout=30)
    assert not app.exception
    assert _scatter(app)["layout"]["xaxis"]["type"] == "linear"
    next(box for box in app.selectbox if box.label == "Numerator condition").set_value("control").run(timeout=30)
    assert not app.exception
    assert "treated" in _scatter(app)["layout"]["xaxis"]["title"]["text"]
    assert all(path.read_bytes() == content for path, content in before.items())


def test_single_condition_lfq_offers_direct_acquisition_scatter(saved_experiment):  # noqa: F811
    from src.naseweis_comparison import build_comparison

    folder, root, manifest = saved_experiment
    for sample in manifest["samples"]:
        sample["condition"] = "one condition"
    (root / "experiment.json").write_text(json.dumps(manifest))
    build_comparison(root, manifest["samples"])
    app = app_for(folder, mode="results").run(timeout=30)
    assert not app.exception
    assert any("not condition means" in text.value for text in app.caption)
    assert not any(box.label == "Numerator condition" for box in app.selectbox)
    assert any(box.label == "Y acquisition" for box in app.selectbox)
    assert "run1.mzML" in _scatter(app)["layout"]["xaxis"]["title"]["text"]
    assert "run2.mzML" in _scatter(app)["layout"]["yaxis"]["title"]["text"]


@pytest.mark.parametrize("means", [(float("nan"), 2.0), (0.0, 2.0)])
def test_no_paired_measurements_explains_why_no_scatter_can_be_drawn(means):
    harness = """
import streamlit as st
from src.naseweis_result_plots_ui import render_lfq_scatter
render_lfq_scatter(st.session_state.contrasts, 'treated', 'control', 'test',
                   measure='intensity', normalization='none')
"""
    contrasts = pd.DataFrame({"numerator_mean": [means[0]], "denominator_mean": [means[1]]}, index=["ACGU"])
    app = AppTest.from_string(harness)
    app.session_state.contrasts = contrasts
    app.run(timeout=30)
    assert not app.exception
    assert not app.get("plotly_chart")
    assert any("No oligonucleotides have positive measurements" in info.value for info in app.info)
    assert any("1 without positive measurements" in text.value for text in app.caption)
