"""Wrong-order input journeys remain recoverable without changing science settings."""

# Imported shared fixtures intentionally match test argument names.
# ruff: noqa: F811

import json
from pathlib import Path

import pytest

from test_gui import app, topp_schemas  # noqa: F401


def configure(app):
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert not app.exception


def load_example(app):
    app.run(timeout=30)
    next(button for button in app.button if button.label == "Load example dataset").click().run(timeout=30)
    assert not app.exception
    configure(app)
    folder = Path(app.session_state.workspace) / "naseweis"
    return folder, folder / "input-files/mzML-files/NME1.mzML"


def test_configure_before_upload_explains_how_to_add_inputs_and_keeps_lfq_off(app, topp_schemas):
    app.run(timeout=30)
    configure(app)
    assert any("No MS/MS acquisitions files yet" in item.value for item in app.caption)
    assert any("No RNA sequence database files yet" in item.value for item in app.caption)
    assert any(button.key == "_browse-mzML-files" for button in app.button)
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert app.checkbox(key="naseweis-param-enable-lfq").value is False


def test_available_example_files_are_not_reselected_after_the_user_clears_them(app, topp_schemas):
    folder, source = load_example(app)
    app.multiselect(key="naseweis-param-mzML-files").set_value([]).run(timeout=30)
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    assert source.is_file()
    configure(app)
    assert app.multiselect(key="naseweis-param-mzML-files").value == []
    assert json.loads((folder / "params.json").read_text())["mzML-files"] == []


def test_deleted_selected_acquisition_is_removed_from_configure_selection(app, topp_schemas):
    folder, source = load_example(app)
    source.unlink()
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert any("No MS/MS acquisitions files yet" in item.value for item in app.caption)
    assert json.loads((folder / "params.json").read_text())["mzML-files"] == []
    assert not any(str(source) in box.value for box in app.multiselect if isinstance(box.value, list))


def test_deleted_external_selection_does_not_fall_back_to_another_available_acquisition(app, topp_schemas, tmp_path):
    folder, source = load_example(app)
    external = tmp_path / "external.mzML"
    external.write_text("temporary external acquisition")
    (folder / "input-files/mzML-files/external_files.txt").write_text(str(external) + "\n")
    configure(app)
    app.multiselect(key="naseweis-param-mzML-files").set_value([str(external)]).run(timeout=30)
    external.unlink()
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert app.multiselect(key="naseweis-param-mzML-files").value == []
    assert source.is_file()
    assert json.loads((folder / "params.json").read_text())["mzML-files"] == []
    assert any("external.mzML" in item.value and "no longer available" in item.value for item in app.warning)


def test_deleted_selected_database_is_not_silently_replaced_by_another_database(app, topp_schemas):
    folder, _source = load_example(app)
    original = folder / "input-files/fasta-files/NME1_Target.fasta"
    other = original.with_name("other.fasta")
    other.write_text(">other\nACGU\n")
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    original.unlink()
    configure(app)
    assert app.selectbox(key="naseweis-param-fasta-files").value is None
    assert any(original.name in item.value and "no longer available" in item.value for item in app.warning)
    assert json.loads((folder / "params.json").read_text())["fasta-files"] is None


@pytest.mark.parametrize("selection", [42, {"path": "bad.mzML"}, [{"path": "bad.mzML"}], ["ok.mzML", None], [""], ["   "]])
def test_invalid_imported_acquisition_selection_has_targeted_recovery(app, topp_schemas, selection):
    folder, source = load_example(app)
    params_file = folder / "params.json"
    saved = json.loads(params_file.read_text())
    saved.update({
        "mzML-files": selection,
        "lfq-normalization": "none",
        "NucleicAcidSearchEngine": {"oligo:missed_cleavages": 3},
    })
    params_file.write_text(json.dumps(saved))
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert any("acquisition selection is invalid" in item.value for item in app.error)
    reset = next(button for button in app.button if button.label == "Reset acquisition selection")
    reset.click().run(timeout=30)
    assert not app.exception
    assert app.multiselect(key="naseweis-param-mzML-files").value == []
    after = json.loads(params_file.read_text())
    assert after["mzML-files"] == []
    assert after["lfq-normalization"] == "none"
    assert after["NucleicAcidSearchEngine"]["oligo:missed_cleavages"] == 3
    assert source.is_file()
    app.multiselect(key="naseweis-param-mzML-files").set_value([str(source)]).run(timeout=30)
    assert not app.exception
    assert json.loads(params_file.read_text())["mzML-files"] == [str(source)]


def test_lfq_settings_survive_leaving_configure_before_an_invalid_run(app, topp_schemas):
    folder, _source = load_example(app)
    app.checkbox(key="naseweis-param-enable-lfq").check().run(timeout=30)
    app.selectbox(key="naseweis-param-lfq-normalization").set_value("none").run(timeout=30)
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    key = "naseweis-TOPP-FeatureFinderMetaboIdent:1:extract:mz_window"
    app.number_input(key=key).set_value(17.0).run(timeout=30)
    app.switch_page("content/topp_workflow_execution.py").run(timeout=30)
    assert not app.exception
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    configure(app)
    assert app.checkbox(key="naseweis-param-enable-lfq").value is True
    assert app.selectbox(key="naseweis-param-lfq-normalization").value == "none"
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert app.number_input(key=key).value == 17.0
    saved = json.loads((folder / "params.json").read_text())
    assert saved["enable-lfq"] is True
    assert saved["FeatureFinderMetaboIdent"]["extract:mz_window"] == 17.0


def test_upload_submit_without_a_file_keeps_the_page_and_inputs_usable(app, monkeypatch):
    monkeypatch.setattr("src.workflow.StreamlitUI.native_dialog_available", lambda: False)
    app.run(timeout=30)
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    assert not app.exception
    button = next(button for button in app.button if button.label == "Add **MS/MS spectra**")
    button.click().run(timeout=30)
    assert not app.exception
    assert any("Nothing to add" in item.value for item in app.error)
    assert any(button.label == "Load example dataset" for button in app.button)


def test_method_summary_can_describe_mixed_nested_import_values_without_changing_them(tmp_path):
    from src.workflow.StreamlitUI import StreamlitUI

    existing = tmp_path / "existing.mzML"
    existing.touch()
    ui = object.__new__(StreamlitUI)
    ui.params = {
        "mzML-files": [str(existing), None, 3.5, {"bad": [True, "missing.mzML"]}],
        "scientific_setting": [True, 2, 3.5, ["nested"]],
    }
    before = json.dumps(ui.params)
    summary = ui.non_default_params_summary()
    assert "existing.mzML" in summary
    assert "null" in summary
    assert "3.5" in summary
    assert '"bad"' in summary
    assert "missing.mzML" in summary
    assert json.dumps(ui.params) == before


@pytest.mark.parametrize("defaults", [None, ["bad"], {"NucleicAcidSearchEngine": [None]}, {"NucleicAcidSearchEngine": "bad"}])
def test_method_summary_labels_malformed_defaults_instead_of_crashing(defaults):
    from src.workflow.StreamlitUI import StreamlitUI

    ui = object.__new__(StreamlitUI)
    ui.params = {"_defaults": defaults}
    before = json.dumps(ui.params)
    assert "invalid format" in ui.non_default_params_summary()
    assert json.dumps(ui.params) == before


def test_method_summary_tolerates_imported_strings_that_are_not_valid_paths():
    from src.workflow.StreamlitUI import StreamlitUI

    ui = object.__new__(StreamlitUI)
    too_long = "x" * 5000
    ui.params = {"scientific_text": too_long, "scientific_list": [too_long, "bad\x00path"]}
    before = json.dumps(ui.params)
    summary = ui.non_default_params_summary()
    assert too_long in summary
    assert "bad\x00path" in summary
    assert json.dumps(ui.params) == before
