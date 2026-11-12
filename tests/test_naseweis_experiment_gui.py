"""Experiment selection, saved comparisons, and acquisition drill-down UI."""

import io
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

HARNESS = """
from pathlib import Path
from types import SimpleNamespace
import streamlit as st
from src.workflow.ParameterManager import ParameterManager
from src.workflow.StreamlitUI import StreamlitUI
from src.naseweis_experiment_ui import (
    configure_comparison, configure_experiment_inputs, render_experiment_results,
)

folder = Path(st.session_state.folder)
manager = ParameterManager(folder)
params = manager.get_parameters_from_json()
workflow = SimpleNamespace(
    workflow_dir=folder, results_dir=folder / "results", params=params,
    parameter_manager=manager, ui=StreamlitUI(folder, None, None, manager),
)
if st.session_state.mode == "configure":
    st.toggle("Advanced", key="advanced")
    configure_experiment_inputs(workflow)
    if st.session_state.get("comparison_controls", False):
        configure_comparison(workflow)
else:
    st.session_state.opened_sample = None
    def sample_details(sample):
        st.session_state.opened_sample = sample["id"]
        st.write("Sample result: " + sample["file"])
    if not render_experiment_results(workflow, sample_details):
        st.info("Legacy single-acquisition results")
"""


@pytest.fixture
def setup_inputs(tmp_path):
    folder = tmp_path / "experiment"
    inputs = folder / "input-files/mzML-files"
    inputs.mkdir(parents=True)
    files = [inputs / "control_R1.mzML", inputs / "treated_R1.mzML"]
    for path in files:
        path.touch()
    return folder, [str(path) for path in files]


def app_for(folder, *, mode="configure", params=None, comparison=False):
    folder.mkdir(parents=True, exist_ok=True)
    if params is not None:
        (folder / "params.json").write_text(json.dumps(params))
    app = AppTest.from_string(HARNESS)
    app.session_state.folder = str(folder)
    app.session_state.mode = mode
    app.session_state.location = "local"
    app.session_state.local_dir = ""
    app.session_state.previous_dir = str(folder)
    app.session_state.settings = {"data_import": {"enabled": False}}
    app.session_state.comparison_controls = comparison
    return app


@pytest.mark.parametrize("legacy_widget", [False, True])
def test_legacy_scalar_selection_becomes_a_persisted_list_without_enabling_lfq(
    setup_inputs, legacy_widget
):
    folder, files = setup_inputs
    app = app_for(
        folder,
        params={"mzML-files": files[0], "enable-lfq": False, "scientific-option": 23},
    )
    if legacy_widget:
        app.session_state["experiment-param-mzML-files"] = files[0]
    app.run(timeout=30)
    assert not app.exception
    assert app.multiselect(key="experiment-param-mzML-files").value == [files[0]]
    saved = json.loads((folder / "params.json").read_text())
    assert saved["mzML-files"] == [files[0]]
    assert saved["enable-lfq"] is False
    assert saved["scientific-option"] == 23
    assert any("Samples and conditions" in text.value for text in app.markdown)


def test_all_selected_files_have_design_rows_and_empty_selection_stays_empty(
    setup_inputs,
):
    from src.workflow.ExperimentalDesign import load_design

    folder, files = setup_inputs
    app = app_for(folder, params={"mzML-files": [files[0]], "enable-lfq": False})
    app.run(timeout=30)
    app.multiselect(key="experiment-param-mzML-files").set_value(files).run(timeout=30)
    assert not app.exception
    design = load_design(folder, files)
    assert design["file"].tolist() == [Path(path).name for path in files]
    assert list(design.columns) == ["file", "condition", "replicate", "fraction"]
    assert design["replicate"].ge(1).all()
    assert design["fraction"].eq(1).all()
    assert any("Every selected acquisition" in caption.value for caption in app.caption)
    app.multiselect(key="experiment-param-mzML-files").set_value([]).run(timeout=30)
    app.run(timeout=30)
    assert not app.exception
    assert app.multiselect(key="experiment-param-mzML-files").value == []
    assert json.loads((folder / "params.json").read_text())["mzML-files"] == []


def test_design_edits_are_saved_with_explicit_conditions_and_replicates(setup_inputs):
    from src.workflow.ExperimentalDesign import load_design

    folder, files = setup_inputs
    app = app_for(folder, params={"mzML-files": files})
    app.run(timeout=30)
    editor = next(
        element for element in app.dataframe if "condition" in element.value.columns
    )
    states = app._tree.get_widget_states()
    states.widgets.add(
        id=editor.proto.id,
        string_value=json.dumps(
            {
                "edited_rows": {
                    0: {"condition": "untreated", "replicate": 2},
                    1: {"condition": "exposed", "replicate": 3},
                },
                "added_rows": [],
                "deleted_rows": [],
            }
        ),
    )
    app._run(states, timeout=30)
    assert not app.exception
    design = load_design(folder, files)
    assert design["condition"].tolist() == ["untreated", "exposed"]
    assert design["replicate"].tolist() == [2, 3]


def test_comparison_normalization_is_visible_and_measure_requires_advanced(
    setup_inputs,
):
    folder, files = setup_inputs
    app = app_for(
        folder, params={"mzML-files": files, "enable-lfq": True}, comparison=True
    )
    app.run(timeout=30)
    assert not app.exception
    assert app.selectbox(key="experiment-param-lfq-normalization").value == "total"
    assert not any(box.key == "experiment-param-lfq-measure" for box in app.selectbox)
    app.selectbox(key="experiment-param-lfq-normalization").set_value("none").run(
        timeout=30
    )
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    app.selectbox(key="experiment-param-lfq-measure").set_value("model_area").run(
        timeout=30
    )
    assert not app.exception
    saved = json.loads((folder / "params.json").read_text())
    assert saved["lfq-normalization"] == "none"
    assert saved["lfq-measure"] == "model_area"
    assert saved["enable-lfq"] is True


@pytest.mark.parametrize("key", ["lfq-measure", "lfq-normalization"])
def test_unknown_comparison_options_are_reported_without_silent_replacement(
    setup_inputs, key
):
    folder, files = setup_inputs
    app = app_for(
        folder, params={"mzML-files": files, key: "unsupported-method"}, comparison=True
    )
    app.run(timeout=30)
    assert not app.exception
    assert any("Unsupported" in error.value for error in app.error)
    assert json.loads((folder / "params.json").read_text())[key] == "unsupported-method"


@pytest.fixture
def saved_experiment(tmp_path):
    from src.naseweis_comparison import build_comparison

    folder = tmp_path / "experiment"
    root = folder / "results"
    root.mkdir(parents=True)
    samples = []
    areas = [(10, 30), (20, 20), (40, None), (20, 20)]
    for position, (area_a, area_b) in enumerate(areas, 1):
        identity = f"sample_{position:03d}"
        samples.append(
            {
                "id": identity,
                "file": f"run{position}.mzML",
                "path": str(tmp_path / f"run{position}.mzML"),
                "condition": "control" if position <= 2 else "treated",
                "replicate": position if position <= 2 else position - 2,
                "fraction": 1,
                "status": "finished",
            }
        )
        directory = root / "samples" / identity / "results/quantification"
        directory.mkdir(parents=True)
        content = "#rt\tmz\tintensity\tcharge\tlabel\tmodel_area\n"
        content += f"10\t500\t{area_a}\t-2\tACGU\t{area_a * 2}\n"
        if area_b is not None:
            content += f"11\t600\t{area_b}\t-3\tUUUU\t{area_b * 2}\n"
        (directory / "features.tsv").write_text(content)
    manifest = {
        "schema": 1,
        "lfq_enabled": True,
        "status": "finished",
        "comparison_status": "finished",
        "samples": samples,
    }
    (root / "experiment.json").write_text(json.dumps(manifest))
    build_comparison(root, samples, measure="intensity", normalization="total")
    return folder, root, manifest


def test_saved_experiment_shows_full_matrices_missing_cells_and_condition_comparison(
    saved_experiment,
):
    folder, _root, _manifest = saved_experiment
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    assert [tab.label for tab in app.tabs][:4] == [
        "Experiment summary",
        "Sequence coverage",
        "Quantification",
        "Sample details",
    ]
    matrix = next(
        element.value
        for element in app.dataframe
        if element.value.index.name == "label"
        and any("run1.mzML" in str(column) for column in element.value.columns)
    )
    assert matrix.shape == (2, 4)
    assert matrix.loc["ACGU"].tolist() == [0.25, 0.5, 1.0, 0.5]
    assert pd.isna(matrix.loc["UUUU"].iloc[2])
    contrast = next(
        element.value
        for element in app.dataframe
        if "log2_fold_change" in element.value.columns
    )
    assert list(contrast.columns[:2]) == ["log2_fold_change", "ratio"]
    displayed = next(
        element
        for element in app.dataframe
        if "log2_fold_change" in element.value.columns
    )
    labels = json.loads(displayed.proto.columns)
    assert labels["log2_fold_change"]["label"] == "Log2 fold change"
    assert labels["ratio"]["label"] == "Ratio"
    assert contrast.loc["ACGU", "ratio"] == 2
    assert contrast.loc["ACGU", "log2_fold_change"] == 1
    assert contrast.loc["ACGU", "numerator_n"] == 2
    assert contrast.loc["ACGU", "denominator_n"] == 2
    assert not any("p-value" in column.lower() for column in contrast.columns)
    downloads = {element.proto.label for element in app.get("download_button")}
    assert downloads >= {
        "Download abundance matrix (TSV)",
        "Download raw areas (TSV)",
        "Download sample design (TSV)",
        "Download comparison settings (JSON)",
    }
    assert app.get("plotly_chart")
    assert app.session_state.opened_sample == "sample_001"
    next(box for box in app.selectbox if box.label == "Numerator condition").set_value(
        "control"
    ).run(timeout=30)
    assert not app.exception
    assert (
        next(box for box in app.selectbox if box.label == "Denominator condition").value
        == "treated"
    )
    reversed_contrast = next(
        element.value
        for element in app.dataframe
        if "log2_fold_change" in element.value.columns
    )
    assert reversed_contrast.loc["ACGU", "ratio"] == 0.5
    next(box for box in app.selectbox if box.label == "Acquisition details").set_value(
        "sample_004"
    ).run(timeout=30)
    assert not app.exception
    assert app.session_state.opened_sample == "sample_004"


def test_identification_only_experiment_has_no_quantification_controls(
    saved_experiment,
):
    folder, root, manifest = saved_experiment
    manifest.update(lfq_enabled=False, comparison_status="not_requested")
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    assert any("identification-only" in info.value for info in app.info)
    assert not any(box.label == "Numerator condition" for box in app.selectbox)
    assert not app.get("plotly_chart")


def test_failed_sample_without_results_is_a_graceful_drilldown(saved_experiment):
    folder, root, manifest = saved_experiment
    manifest["status"] = "failed"
    manifest["samples"].append(
        {
            "id": "sample_005",
            "file": "failed.mzML",
            "path": str(folder / "failed.mzML"),
            "condition": "treated",
            "replicate": 3,
            "fraction": 1,
            "status": "failed",
            "error": "NASE failed",
        }
    )
    (root / "experiment.json").write_text(json.dumps(manifest))
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    assert any("incomplete acquisitions" in warning.value for warning in app.warning)
    assert any("incomplete experiment" in warning.value for warning in app.warning)
    assert not any(
        "Download abundance matrix" in item.proto.label
        for item in app.get("download_button")
    )
    assert not any(box.label == "Numerator condition" for box in app.selectbox)
    assert not app.get("plotly_chart")
    next(box for box in app.selectbox if box.label == "Acquisition details").set_value(
        "sample_005"
    ).run(timeout=30)
    assert not app.exception
    assert app.session_state.opened_sample is None
    assert any("No sample result files" in info.value for info in app.info)


def test_missing_comparison_does_not_hide_sample_identifications(saved_experiment):
    folder, root, _manifest = saved_experiment
    (root / "comparison/comparison.json").unlink()
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    assert any("not available" in warning.value for warning in app.warning)
    assert app.session_state.opened_sample == "sample_001"


def test_legacy_results_fall_through_and_malformed_manifest_reports_error(tmp_path):
    folder = tmp_path / "experiment"
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    assert any("Legacy single" in info.value for info in app.info)
    (folder / "results").mkdir()
    (folder / "results/experiment.json").write_text("{broken")
    app.run(timeout=30)
    assert not app.exception
    assert any("Could not read this experiment" in error.value for error in app.error)


@pytest.mark.parametrize(
    "field,replacement",
    [
        ("file", "other.mzML"),
        ("path", "/other/run.mzML"),
        ("condition", "different"),
        ("replicate", 8),
        ("fraction", 2),
    ],
)
def test_saved_comparison_cannot_use_a_design_different_from_manifest(
    saved_experiment, field, replacement
):
    folder, root, _manifest = saved_experiment
    path = root / "comparison/sample-design.tsv"
    design = pd.read_csv(path, sep="\t")
    design.loc[0, field] = replacement
    design.to_csv(path, sep="\t", index=False)
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    assert any("does not match this experiment" in error.value for error in app.error)
    assert not any(
        "Download abundance matrix" in item.proto.label
        for item in app.get("download_button")
    )
    assert not app.get("plotly_chart")
    assert app.session_state.opened_sample == "sample_001"


@pytest.mark.parametrize(
    "field,value",
    [("replicate", 1.5), ("replicate", 0), ("fraction", -1), ("fraction", 1.5)],
)
def test_invalid_editor_numbers_are_not_silently_normalized(setup_inputs, field, value):
    from src.workflow.ExperimentalDesign import design_path

    folder, files = setup_inputs
    app = app_for(folder, params={"mzML-files": files, "enable-lfq": False})
    app.run(timeout=30)
    before = design_path(folder).read_bytes()
    editor = next(
        element for element in app.dataframe if "condition" in element.value.columns
    )
    states = app._tree.get_widget_states()
    states.widgets.add(
        id=editor.proto.id,
        string_value=json.dumps(
            {"edited_rows": {0: {field: value}}, "added_rows": [], "deleted_rows": []}
        ),
    )
    app._run(states, timeout=30)
    assert not app.exception
    assert any("positive integer" in error.value for error in app.error)
    assert design_path(folder).read_bytes() == before
    assert any(button.label == "Discard invalid design input" for button in app.button)
    assert json.loads((folder / "params.json").read_text())["enable-lfq"] is False


def test_invalid_stored_design_has_explicit_reset_recovery(setup_inputs):
    from src.workflow.ExperimentalDesign import design_path, load_design

    folder, files = setup_inputs
    invalid = (
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\tMSstats_Condition\tMSstats_BioReplicate\n1\t1.5\t"
        + files[0]
        + "\t1\t1\tcontrol\t1\n"
    )
    design_path(folder).write_text(invalid)
    app = app_for(folder, params={"mzML-files": files, "enable-lfq": False})
    app.run(timeout=30)
    assert not app.exception
    assert any("Fraction" in error.value for error in app.error)
    assert design_path(folder).read_text() == invalid
    next(
        button for button in app.button if button.label == "Reset sample design"
    ).click().run(timeout=30)
    assert not app.exception
    assert not app.error
    assert load_design(folder, files)["fraction"].eq(1).all()
    assert json.loads((folder / "params.json").read_text())["enable-lfq"] is False


def test_invalid_pending_design_upload_preserves_valid_saved_design(setup_inputs):
    from src.workflow.ExperimentalDesign import design_path

    folder, files = setup_inputs
    app = app_for(folder, params={"mzML-files": files})
    app.run(timeout=30)
    before = design_path(folder).read_bytes()
    app.session_state["_design-experiment-mzML-files-import"] = SimpleNamespace(
        getvalue=lambda: (
            b"source name\tcomment[data file]\tcharacteristics[biological replicate]\nsample\tcontrol_R1.mzML\t1.5\n"
        )
    )
    app.run(timeout=30)
    assert not app.exception
    assert any("biological replicate" in error.value for error in app.error)
    assert design_path(folder).read_bytes() == before
    next(
        button
        for button in app.button
        if button.label == "Discard invalid design input"
    ).click().run(timeout=30)
    assert not app.exception
    assert not app.error


@pytest.mark.parametrize("lfq", [False, "false", True])
def test_condition_and_fraction_notices_are_specific_to_enabled_lfq(setup_inputs, lfq):
    from src.workflow.ExperimentalDesign import default_design, save_design

    folder, files = setup_inputs
    design = default_design(files)
    design["fraction"] = 2
    save_design(folder, design, files)
    app = app_for(folder, params={"mzML-files": files, "enable-lfq": lfq})
    app.run(timeout=30)
    assert not app.exception
    notices = [warning.value for warning in app.warning]
    assert any("Condition comparisons require" in message for message in notices) is (
        lfq is True
    )
    assert any("unfractionated samples" in message for message in notices) is (
        lfq is True
    )


def test_downloaded_contrast_keeps_selected_direction_and_area_settings(
    saved_experiment, monkeypatch
):
    captured = {}
    native_button = st.download_button

    def download(*args, **kwargs):
        if args and args[0] == "Download condition comparison (TSV)":
            captured["tsv"] = args[1]
        return native_button(*args, **kwargs)

    monkeypatch.setattr(st, "download_button", download)
    folder, _root, _manifest = saved_experiment
    app = app_for(folder, mode="results")
    app.run(timeout=30)
    assert not app.exception
    table = pd.read_csv(io.BytesIO(captured["tsv"]), sep="\t")
    assert list(table.columns[-2:]) == ["ratio", "log2_fold_change"]
    assert table.loc[table["label"] == "ACGU", "ratio"].item() == 2
    assert table["numerator_condition"].eq("treated").all()
    assert table["denominator_condition"].eq("control").all()
    assert table["measure"].eq("intensity").all()
    assert table["normalization"].eq("total").all()
    next(box for box in app.selectbox if box.label == "Numerator condition").set_value(
        "control"
    ).run(timeout=30)
    table = pd.read_csv(io.BytesIO(captured["tsv"]), sep="\t")
    assert table["numerator_condition"].eq("control").all()
    assert table["denominator_condition"].eq("treated").all()
