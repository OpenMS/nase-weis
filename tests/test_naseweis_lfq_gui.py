"""Optional LFQ controls and saved-run feature views use the real app UI."""

import json
from pathlib import Path

import pytest

from test_gui import (  # noqa: F401
    _downloads,
    _finished_run,
    _results_page,
    app,
    topp_schemas,
)


def test_lfq_is_optional_and_tuning_requires_advanced(app, topp_schemas):  # noqa: F811
    app.run(timeout=30)
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert [tab.label for tab in app.tabs][:2] == ["Identification", "Quantification"]
    toggle = "naseweis-param-enable-lfq"
    assert app.checkbox(key=toggle).value is False
    prefix = "naseweis-TOPP-FeatureFinderMetaboIdent:1:"
    assert not any(str(item.key).startswith(prefix) for item in app.number_input)
    app.checkbox(key=toggle).check().run(timeout=30)
    assert not app.exception
    assert not any(str(item.key).startswith(prefix) for item in app.number_input)
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert not app.exception
    assert app.number_input(key=prefix + "extract:mz_window").value == 50.0
    assert app.number_input(key=prefix + "extract:rt_window").value == 240.0
    assert app.number_input(key=prefix + "extract:isotope_pmin").value == 0.05
    assert app.selectbox(key=prefix + "model:type").value == "none"
    app.number_input(key=prefix + "extract:mz_window").set_value(25.0).run(timeout=30)
    app.switch_page("content/naseweis.py").run(timeout=30)
    app.toggle(key="advanced").set_value(False).run(timeout=30)
    assert not any(str(item.key).startswith(prefix) for item in app.number_input)
    app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    app.switch_page("content/naseweis.py").run(timeout=30)
    assert app.checkbox(key=toggle).value is True
    saved = json.loads((Path(app.session_state.workspace) / "naseweis/params.json").read_text())
    assert saved["enable-lfq"] is True
    assert saved["FeatureFinderMetaboIdent"]["extract:mz_window"] == 25.0
    app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert app.number_input(key=prefix + "extract:mz_window").value == 25.0
    app.checkbox(key=toggle).uncheck().run(timeout=30)
    assert not any(str(item.key).startswith(prefix) for item in app.number_input)


def _lfq_run(app, *, table=None, enabled=True):  # noqa: F811
    results = _finished_run(app)
    snapshot = results.parent / "params.json"
    params = json.loads(snapshot.read_text()) | {"enable-lfq": enabled}
    snapshot.write_text(json.dumps(params))
    quantification = results / "quantification"
    quantification.mkdir()
    (quantification / "features.featureXML").write_text("<featureMap/>")
    if table is not None:
        (quantification / "features.tsv").write_text(table)
    return results


def test_lfq_table_and_native_downloads_come_from_saved_run(app):  # noqa: F811
    results = _lfq_run(app, table="#rt\tmz\tintensity\tcharge\tlabel\n60\t500.2\t1234.5\t-2\tACGU\n")
    _results_page(app)
    assert "Quantification" in [tab.label for tab in app.tabs]
    assert len(app.dataframe) == 2
    table = app.dataframe[1].value
    assert table.loc[0, "label"] == "ACGU"
    assert table.loc[0, "intensity"] == 1234.5
    assert _downloads(app) >= {"Download LFQ features (featureXML)", "Download LFQ feature table (TSV)"}
    # Selecting a search-only run must not retain the previous run's table.
    newer = _finished_run(app, "Search only", sequence="AAAA").parent.name
    _results_page(app)
    picker = next(box for box in app.selectbox if box.label == "Run")
    picker.select(newer).run(timeout=30)
    assert len(app.dataframe) == 1
    assert "Download LFQ feature table (TSV)" not in _downloads(app)
    assert any("not enabled" in item.value for item in app.info)
    picker = next(box for box in app.selectbox if box.label == "Run")
    picker.select(results.parent.name).run(timeout=30)
    assert app.dataframe[1].value.loc[0, "label"] == "ACGU"


@pytest.mark.parametrize("table", [None, "#rt\tmz\tintensity\tcharge\tlabel\n", "broken table\n"])
def test_empty_missing_or_corrupt_lfq_keeps_identification_results(app, table):  # noqa: F811
    _lfq_run(app, table=table)
    _results_page(app)
    assert not app.exception
    assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
    assert "Download mzTab" in _downloads(app)
    assert "Download LFQ features (featureXML)" in _downloads(app)
    if table is None:
        assert any("no completed feature table" in item.value for item in app.warning)
    elif table.startswith("broken"):
        assert any("Could not read this run's feature table" in item.value for item in app.error)
    else:
        assert any("No MS1 features" in item.value for item in app.info)


def test_lfq_without_snapshot_has_an_explanation(app):  # noqa: F811
    results = _finished_run(app)
    (results.parent / "params.json").unlink()
    _results_page(app)
    assert not app.exception
    assert any("not enabled" in item.value for item in app.info)


def test_failed_lfq_run_keeps_partial_identifications_accessible(app):  # noqa: F811
    from src.workflow.RunManager import RunManager

    _finished_run(app, "Older complete search")
    results = _lfq_run(app)
    runs = RunManager(results.parent.parent.parent)
    runs.finish(results.parent.name, "failed")
    logs = results.parent / "logs"
    logs.mkdir(exist_ok=True)
    (logs / "all.log").write_text("NASE succeeded; FeatureFinderMetaboIdent failed.\n")
    # The page offers the failed run alongside finished runs, labelled clearly.
    _results_page(app)
    picker = next(box for box in app.selectbox if box.label == "Run")
    assert any("Failed / partial results" in label for label in picker.options)
    picker.select(results.parent.name).run(timeout=30)
    assert app.dataframe[0].value.loc[0, "Sequence"] == "ACGU"
    assert "Download mzTab" in _downloads(app)
    assert "Download run log" in _downloads(app)
    assert any("This run failed" in item.value for item in app.warning)
    assert any("no completed feature table" in item.value for item in app.warning)
    # A full page rerender must preserve the failed-run selection.
    app.switch_page("content/naseweis.py").run(timeout=30)
    _results_page(app)
    assert any("This run failed" in item.value for item in app.warning)
    assert "Download mzTab" in _downloads(app)


def test_failed_run_without_artifacts_does_not_appear_as_a_result(app):  # noqa: F811
    from src.workflow.RunManager import RunManager

    results = _finished_run(app)
    (results / "mztab_results/search.mzTab").unlink()
    runs = RunManager(results.parent.parent.parent)
    runs.finish(results.parent.name, "failed")
    _results_page(app)
    assert any("No finished runs" in item.value for item in app.info)
