"""Full-app experiment controls and saved sample evidence stay connected."""

# Imported pytest fixtures intentionally share test parameter names.
# ruff: noqa: F811

import json
from pathlib import Path

import pytest
import streamlit as st
from test_naseweis_annotations import store_annotations
from test_naseweis_viewers_gui import (
    NATIVE_ID,
    acquisition,
    peak_coordinates,
    plots,
    result_content,
)

from src.naseweis_comparison import build_comparison
from src.naseweis_experiment import sample_results, sample_snapshot
from src.naseweis_provenance import acquisition_fingerprint
from test_gui import _finished_run, _results_page, app, topp_schemas  # noqa: F401


def configure(app):
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert not app.exception


def seed_inputs(app):
    app.run(timeout=30)
    next(
        button for button in app.button if button.label == "Load example dataset"
    ).click().run(timeout=30)
    folder = Path(app.session_state.workspace) / "naseweis"
    second = acquisition(folder / "input-files/mzML-files/second.mzML", second=True)
    configure(app)
    first = folder / "input-files/mzML-files/NME1.mzML"
    return folder, [str(first), str(second)]


def test_full_configure_selection_deselection_and_clear_keep_lfq_off(app, topp_schemas):
    folder, paths = seed_inputs(app)
    key = "naseweis-param-mzML-files"
    assert app.multiselect(key=key).value == [paths[0]]
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False
    app.multiselect(key=key).set_value(paths).run(timeout=30)
    assert not app.exception
    design = next(
        element.value
        for element in app.dataframe
        if "condition" in element.value.columns
    )
    assert design["file"].tolist() == [Path(path).name for path in paths]
    app.multiselect(key=key).set_value([paths[1]]).run(timeout=30)
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert app.multiselect(key=key).value == [paths[1]]
    app.multiselect(key=key).set_value([]).run(timeout=30)
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert app.multiselect(key=key).value == []
    assert json.loads((folder / "params.json").read_text())["mzML-files"] == []
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert not app.exception
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False
    assert not any(
        button.key == "naseweis-rna-comparison-preset" for button in app.button
    )
    assert not any(box.key == "naseweis-param-lfq-measure" for box in app.selectbox)


def test_configure_and_method_summary_work_after_experiment_design_snapshot(app, topp_schemas):
    folder, paths = seed_inputs(app)
    params_file = folder / "params.json"
    params = json.loads(params_file.read_text())
    params["mzML-files"] = paths
    params["_sample_design"] = [
        {"file": path, "condition": "control", "replicate": number, "fraction": 1}
        for number, path in enumerate(paths, 1)
    ]
    params_file.write_text(json.dumps(params))
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert app.multiselect(key="naseweis-param-mzML-files").value == paths
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False
    assert any(button.label == "📑 Method summary" for button in app.get("download_button"))
    assert json.loads(params_file.read_text())["_sample_design"] == params["_sample_design"]


def test_explicit_fitted_preset_persists_and_hiding_controls_does_not_enable_lfq(
    app, topp_schemas
):
    folder, paths = seed_inputs(app)
    app.multiselect(key="naseweis-param-mzML-files").set_value(paths).run(timeout=30)
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False
    app.checkbox(key="naseweis-param-enable-lfq").check().run(timeout=30)
    app.selectbox(key="naseweis-param-lfq-normalization").set_value("none").run(
        timeout=30
    )
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    app.button(key="naseweis-rna-comparison-preset").click().run(timeout=30)
    assert not app.exception
    prefix = "naseweis-TOPP-FeatureFinderMetaboIdent:1:"
    assert app.selectbox(key="naseweis-param-lfq-measure").value == "model_area"
    assert app.selectbox(key="naseweis-param-lfq-normalization").value == "total"
    assert app.selectbox(key=prefix + "model:type").value == "symmetric"
    expected = {
        "extract:mz_window": 10.0,
        "extract:rt_window": 0.0,
        "detect:peak_width": 60.0,
        "extract:n_isotopes": 5,
        "extract:isotope_pmin": 0.0,
    }
    for key, value in expected.items():
        assert app.number_input(key=prefix + key).value == value
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert app.selectbox(key=prefix + "model:type").value == "symmetric"
    saved = json.loads((folder / "params.json").read_text())
    assert (
        saved["FeatureFinderMetaboIdent"] | expected
        == saved["FeatureFinderMetaboIdent"]
    )
    assert saved["FeatureFinderMetaboIdent"]["EMGScoring:init_mom"] is False
    app.checkbox(key="naseweis-param-enable-lfq").uncheck().run(timeout=30)
    app.toggle(key="advanced").set_value(False).run(timeout=30)
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert not app.exception
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False
    saved = json.loads((folder / "params.json").read_text())
    assert saved["enable-lfq"] is False
    assert saved["lfq-measure"] == "model_area"
    assert saved["FeatureFinderMetaboIdent"]["model:type"] == "symmetric"


@pytest.fixture
def experiment_results(app, tmp_path):
    root = _finished_run(app, "Saved RNA experiment")
    workflow_dir = Path(app.session_state.workspace) / "naseweis"
    samples, snapshots, inputs = [], [], []
    for number, sequence in enumerate(["ACGU", "UUUU"], 1):
        identity = f"sample_{number:03d}"
        original = acquisition(
            workflow_dir / "input-files/mzML-files" / f"run-{number}.mzML",
            second=number == 2,
        )
        result = sample_results(root, identity)
        (result / "mztab_results").mkdir(parents=True)
        (result / "mztab_results/search.mzTab").write_text(
            result_content(sequence=sequence)
        )
        annotations = tmp_path / f"annotations-{number}"
        annotations.mkdir()
        native = store_annotations(
            annotations,
            [
                {
                    "native_id": NATIVE_ID,
                    "hits": [{"sequence": sequence, "annotations": []}],
                }
            ],
        )
        (result / "idxml_results").mkdir()
        (result / "idxml_results/search.idXML").write_bytes(native.read_bytes())
        quantification = result / "quantification"
        quantification.mkdir()
        (quantification / "features.tsv").write_text(
            f"#rt\tmz\tintensity\tcharge\tlabel\n10\t500\t{number * 100}\t-2\t{sequence}\n"
        )
        snapshot = sample_snapshot(root, identity)
        snapshot.write_text(
            json.dumps(
                {
                    "mzML-files": str(original),
                    "enable-lfq": True,
                    "FeatureFinderMetaboIdent": {"extract:mz_window": number * 10.0},
                    "_acquisition_fingerprint": acquisition_fingerprint(original),
                }
            )
        )
        logs = snapshot.parent / "logs"
        logs.mkdir()
        (logs / "all.log").write_text(f"sample {number} native log\n")
        samples.append(
            {
                "id": identity,
                "file": original.name,
                "path": str(original),
                "condition": "control" if number == 1 else "treated",
                "replicate": 1,
                "fraction": 1,
                "status": "finished",
            }
        )
        snapshots.append(snapshot)
        inputs.append(original)
    manifest = {
        "schema": 1,
        "lfq_enabled": True,
        "status": "finished",
        "comparison_status": "finished",
        "samples": samples,
    }
    (root / "experiment.json").write_text(json.dumps(manifest))
    build_comparison(root, samples, normalization="total")
    # Neither the master selection nor the live Configure selection identifies
    # one sample's evidence. Their contents must not reach the sample viewers.
    (root.parent / "params.json").write_text(
        json.dumps({"mzML-files": [str(path) for path in inputs], "enable-lfq": True})
    )
    (workflow_dir / "params.json").write_text(
        json.dumps({"mzML-files": [str(inputs[0])], "enable-lfq": False})
    )
    (root.parent / "logs").mkdir(exist_ok=True)
    (root.parent / "logs/all.log").write_text("master log must not be the sample log\n")
    return app, root, inputs, snapshots


@pytest.fixture
def download_bytes(monkeypatch):
    captured = {}
    native_button = st.download_button

    def record(*args, **kwargs):
        label = args[0] if args else kwargs.get("label")
        data = args[1] if len(args) > 1 else kwargs.get("data")
        if isinstance(data, (bytes, str)):
            captured[label] = data
        return native_button(*args, **kwargs)

    monkeypatch.setattr(st, "download_button", record)
    return captured


def test_real_results_sample_picker_uses_own_snapshot_spectra_and_downloads(
    experiment_results, download_bytes
):
    app, root, inputs, snapshots = experiment_results
    _results_page(app)
    assert (
        app.selectbox(
            key=next(
                box.key for box in app.selectbox if box.label == "Acquisition details"
            )
        ).value
        == "sample_001"
    )
    assert json.loads(download_bytes["Download search parameters"])[
        "mzML-files"
    ] == str(inputs[0])
    assert (
        download_bytes["Download mzTab"]
        == (
            sample_results(root, "sample_001") / "mztab_results/search.mzTab"
        ).read_bytes()
    )
    assert download_bytes["Download run log"] == b"sample 1 native log\n"
    next(
        button for button in app.button if button.label == "Open spectrum viewer"
    ).click().run(timeout=30)
    assert not app.exception
    spectrum = next(
        figure
        for figure in plots(app)
        if figure["layout"].get("xaxis", {}).get("title", {}).get("text") == "m/z"
    )
    assert peak_coordinates(spectrum) == [(100.0, 50.0), (200.0, 100.0)]
    next(box for box in app.selectbox if box.label == "Acquisition details").set_value(
        "sample_002"
    ).run(timeout=30)
    assert not app.exception
    assert json.loads(download_bytes["Download search parameters"])[
        "mzML-files"
    ] == str(inputs[1])
    assert download_bytes["Download search parameters"] == snapshots[1].read_bytes()
    assert (
        download_bytes["Download mzTab"]
        == (
            sample_results(root, "sample_002") / "mztab_results/search.mzTab"
        ).read_bytes()
    )
    assert download_bytes["Download run log"] == b"sample 2 native log\n"
    next(
        button for button in app.button if button.label == "Open spectrum viewer"
    ).click().run(timeout=30)
    assert not app.exception
    spectrum = next(
        figure
        for figure in plots(app)
        if figure["layout"].get("xaxis", {}).get("title", {}).get("text") == "m/z"
    )
    assert peak_coordinates(spectrum) == [
        (900.0, 100.0),
        (950.0, pytest.approx(100 / 3)),
    ]
    assert any(
        inputs[1].name in caption.value and NATIVE_ID in caption.value
        for caption in app.caption
    )


def test_sample_fingerprint_failure_does_not_disable_another_samples_viewer(
    experiment_results,
):
    app, _root, inputs, _snapshots = experiment_results
    inputs[0].write_bytes(b"replacement acquisition bytes")
    _results_page(app)
    assert not app.exception
    assert any("acquisition has changed" in info.value.lower() for info in app.info)
    assert not any(button.label == "Open spectrum viewer" for button in app.button)
    next(box for box in app.selectbox if box.label == "Acquisition details").set_value(
        "sample_002"
    ).run(timeout=30)
    next(
        button for button in app.button if button.label == "Open spectrum viewer"
    ).click().run(timeout=30)
    assert not app.exception
    spectrum = next(
        figure
        for figure in plots(app)
        if figure["layout"].get("xaxis", {}).get("title", {}).get("text") == "m/z"
    )
    assert peak_coordinates(spectrum)[0] == (900.0, 100.0)
