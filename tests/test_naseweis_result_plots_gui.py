"""Coverage and comparative plots remain tied to saved scientific evidence."""

import json

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from test_naseweis_experiment_gui import app_for, saved_experiment  # noqa: F401

from src.naseweis_experiment import sample_results
from src.naseweis_result_plots_ui import _coverage_tracks
from tests import naseweis_drawn


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
    """Settings of the condition scatter (an Insight SeriesPlot) the last run drew."""
    return next(item.settings for item in naseweis_drawn.plots(app.drawn)
                if item.key.endswith("-scatter"))


def _coverage_notes(app):
    """The notes of the coverage map: inside the Insight SequenceCoverage, or
    under the static map an openms-insight without that component gets."""
    drawn = [note for item in naseweis_drawn.plots(app.drawn) if item.component == "SequenceCoverage"
             for note in item.settings.get("notes") or []]
    return drawn + [caption.value for caption in app.caption]


def _log(settings):
    return settings["x_label"].startswith("log10") and settings["y_label"].startswith("log10")


@pytest.fixture
def below(monkeypatch):
    """An openms-insight whose coverage map has no settings menu of its own:
    the settings and downloads are in a row below the map."""
    import src.naseweis_result_plots_ui as ui

    monkeypatch.setattr(ui, "_settings_menu", lambda: False)


@pytest.fixture
def menu(monkeypatch):
    """An openms-insight whose coverage map has its own settings and download
    menus; a stand-in for its class where the installed one has none."""
    import openms_insight

    import src.naseweis_insight as insight
    import src.naseweis_result_plots_ui as ui

    if not ui._settings_menu():
        class SequenceCoverage:  # pylint: disable=too-few-public-methods
            _supports_settings = True

        monkeypatch.setattr(openms_insight, "SequenceCoverage", SequenceCoverage, raising=False)
        monkeypatch.setattr(insight, "only_accepted", lambda _cls, **settings: settings)


def _coverage(app):
    """The coverage map (an Insight SequenceCoverage) the last run drew."""
    return next(item for item in naseweis_drawn.plots(app.drawn) if item.component == "SequenceCoverage")


def _set_tracks(app, value):
    """Pick ``value`` under Coverage tracks in the map's settings menu, as its
    browser side does: in the map's state, then a rerun."""
    key = _coverage(app).key.removesuffix("-plot") + "-state"
    state = app.session_state[key]
    state["selections"]["coverage_tracks"] = value
    app.session_state[key] = state
    return app.run(timeout=30)


def _recorded(app, monkeypatch):
    naseweis_drawn.attach(app, monkeypatch)
    return app


def test_unassigned_acquisition_does_not_pool_with_a_matching_condition_name():
    available = [
        ({"file": "one.mzML", "condition": ""}, {"evidence": 1}),
        ({"file": "two.mzML", "condition": "Unassigned · one.mzML"}, {"evidence": 2}),
    ]
    tracks = _coverage_tracks(available, "Condition")
    assert tracks[0][0] != tracks[1][0]
    assert tracks[1][0] == "Unassigned · one.mzML"
    assert [track[1] for track in tracks] == [tables for _, tables in available]


def test_identification_only_experiment_has_condition_and_acquisition_coverage(saved_experiment, monkeypatch, below):  # noqa: F811
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"]:
        _coverage_file(root, sample["id"], legacy=sample["id"] == "sample_004")
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert app.selectbox[0].label == "RNA accession"
    assert app.selectbox[0].value == "NME1"
    assert next(radio for radio in app.radio if radio.label == "Coverage tracks").value == "Condition"
    _download(app, "Download coverage (HTML)")
    _download(app, "Download coverage (SVG)")
    assert not any(item.key.endswith("-scatter") for item in app.drawn)
    assert any("4 of 4 acquisitions" in note for note in _coverage_notes(app))
    next(radio for radio in app.radio if radio.label == "Coverage tracks").set_value("Acquisition").run(timeout=30)
    assert not app.exception
    next(box for box in app.selectbox if box.label == "Bases per line").set_value(25).run(timeout=30)
    assert not app.exception
    assert app.session_state.opened_sample == "sample_001"


def test_coverage_without_insight_component_still_shows_map_and_notes(saved_experiment, monkeypatch):  # noqa: F811
    import openms_insight

    monkeypatch.delattr(openms_insight, "SequenceCoverage", raising=False)
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"]:
        _coverage_file(root, sample["id"])
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert not any(item.component == "SequenceCoverage" for item in app.drawn)
    assert any(text.value.startswith("**Sequence coverage") and "NME1" in text.value for text in app.markdown)
    assert any("4 of 4 acquisitions" in caption.value for caption in app.caption)
    _download(app, "Download coverage (SVG)")


def test_partial_coverage_uses_only_available_saved_results_and_reports_bad_file(saved_experiment, monkeypatch):  # noqa: F811
    folder, root, manifest = saved_experiment
    _coverage_file(root, "sample_001")
    invalid = _coverage_file(root, "sample_002")
    invalid.write_text("NUH\taccession\topt_sequence\nNUC\tbroken\n")
    manifest["status"] = "failed"
    manifest["samples"][2]["status"] = "failed"
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert any("Could not read coverage for run2.mzML" in text.value for text in app.warning)
    assert any("1 of 4 acquisitions" in note for note in _coverage_notes(app))
    assert naseweis_drawn.coverage_shown(app)
    assert not any(radio.label == "Coverage tracks" for radio in app.radio)
    assert not any(item.toolbar for item in app.drawn)
    assert not any(item.key.endswith("-scatter") for item in app.drawn)


def test_conflicting_references_do_not_share_an_axis_but_each_acquisition_remains_accessible(saved_experiment, monkeypatch):  # noqa: F811
    folder, root, _manifest = saved_experiment
    _coverage_file(root, "sample_001")
    _coverage_file(root, "sample_002", sequence="UUUUACGU")
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert any("Conflicting reference sequences" in text.value for text in app.warning)
    assert not any(item.proto.label == "Download coverage (SVG)" for item in app.get("download_button"))
    assert app.session_state.opened_sample == "sample_001"
    assert _log(_scatter(app))


def test_blank_conditions_default_to_acquisition_coverage(saved_experiment, monkeypatch, below):  # noqa: F811
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"][:2]:
        _coverage_file(root, sample["id"])
        sample["condition"] = ""
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    grouping = next(radio for radio in app.radio if radio.label == "Coverage tracks")
    assert grouping.value == "Acquisition"
    grouping.set_value("Condition").run(timeout=30)
    assert not app.exception
    _download(app, "Download coverage (SVG)")


def test_coverage_settings_and_downloads_are_in_the_maps_own_menus(saved_experiment, monkeypatch, menu):  # noqa: F811
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"]:
        _coverage_file(root, sample["id"], legacy=sample["id"] == "sample_004")
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert app.selectbox[0].label == "RNA accession"
    assert not any(box.label == "Bases per line" for box in app.selectbox)
    assert not any(radio.label == "Coverage tracks" for radio in app.radio)
    assert not any(item.proto.label.startswith("Download coverage") for item in app.get("download_button"))
    coverage = _coverage(app)
    assert "block_size" not in coverage.settings  # the lines fit the width
    assert [setting["label"] for setting in coverage.toolbar] == ["Coverage tracks"]
    assert coverage.toolbar[0]["options"] == ["Condition", "Acquisition"]
    assert coverage.toolbar[0]["default"] == "Condition"
    conditions = coverage.settings["track_order"]
    assert any("4 of 4 acquisitions" in note for note in _coverage_notes(app))
    _set_tracks(app, "Acquisition")
    assert not app.exception
    files = _coverage(app).settings["track_order"]
    assert files != conditions and all(name.endswith(".mzML") for name in files)
    assert len(files) == 4
    _set_tracks(app, "Condition")
    assert _coverage(app).settings["track_order"] == conditions
    assert app.session_state.opened_sample == "sample_001"


def test_blank_conditions_default_to_acquisition_tracks_in_the_menu(saved_experiment, monkeypatch, menu):  # noqa: F811
    folder, root, manifest = saved_experiment
    for sample in manifest["samples"][:2]:
        _coverage_file(root, sample["id"])
        sample["condition"] = ""
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert _coverage(app).toolbar[0]["default"] == "Acquisition"
    assert _coverage(app).settings["track_order"] == ["run1.mzML", "run2.mzML"]
    _set_tracks(app, "Condition")
    assert not app.exception
    assert all(name.startswith("Unassigned") for name in _coverage(app).settings["track_order"])


def test_one_acquisition_has_no_track_setting_in_the_menu(saved_experiment, monkeypatch, menu):  # noqa: F811
    folder, root, manifest = saved_experiment
    _coverage_file(root, "sample_001")
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert _coverage(app).toolbar is None


def test_lfq_scatter_uses_saved_means_direction_labels_and_view_only_controls(saved_experiment, monkeypatch):  # noqa: F811
    folder, root, _manifest = saved_experiment
    before = {path: path.read_bytes() for path in root.rglob("*") if path.is_file()}
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    spec = _scatter(app)
    assert _log(spec)
    assert "control" in spec["x_label"]
    assert "treated" in spec["y_label"]
    assert any("2 oligonucleotides measured" in text.value for text in app.caption)
    _download(app, "Download plotted values (TSV)")
    app.multiselect[0].set_value(["ACGU"]).run(timeout=30)
    next(radio for radio in app.radio if radio.label == "Scatter axes").set_value("Linear").run(timeout=30)
    assert not app.exception
    assert not _log(_scatter(app)) and _scatter(app)["x_label"].startswith("Mean")
    next(box for box in app.selectbox if box.label == "Numerator condition").set_value("control").run(timeout=30)
    assert not app.exception
    assert "treated" in _scatter(app)["x_label"]
    assert all(path.read_bytes() == content for path, content in before.items())


def test_single_condition_lfq_offers_direct_acquisition_scatter(saved_experiment, monkeypatch):  # noqa: F811
    from src.naseweis_comparison import build_comparison

    folder, root, manifest = saved_experiment
    for sample in manifest["samples"]:
        sample["condition"] = "one condition"
    (root / "experiment.json").write_text(json.dumps(manifest))
    build_comparison(root, manifest["samples"])
    app = _recorded(app_for(folder, mode="results"), monkeypatch).run(timeout=30)
    assert not app.exception
    assert any("not condition means" in text.value for text in app.caption)
    assert not any(box.label == "Numerator condition" for box in app.selectbox)
    assert any(box.label == "Y acquisition" for box in app.selectbox)
    assert "run1.mzML" in _scatter(app)["x_label"]
    assert "run2.mzML" in _scatter(app)["y_label"]


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
    assert not app.get("iframe") and not app.get("component_instance")
    assert any("No oligonucleotides have positive measurements" in info.value for info in app.info)
    assert any("1 without positive measurements" in text.value for text in app.caption)
