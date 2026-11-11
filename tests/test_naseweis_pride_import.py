"""PRIDE CSV metadata requires an explicit choice of bedRMod mapping upload."""

import io
import json
import zipfile
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from streamlit.testing.v1 import AppTest

from src.data_import import job, ui
from src.data_import.remote import RemoteFile
from test_gui import app as naseweis_app  # noqa: F401
from test_gui import topp_schemas  # noqa: F401

MAPPING = "bedrmod-mapping-files"
MAPPING_NAME = "bedRMod modification mapping (optional)"
HARNESS = """
import streamlit as st
from src.data_import.ui import link_import_section
link_import_section(st.session_state.workflow_dir, st.session_state.import_targets,
                    source=st.session_state.import_source)
"""


def targets_for(directory):
    return {
        "mzML-files": {
            "name": "MS/MS spectra",
            "file_types": ["mzML", "raw", "d"],
            "dir": str(directory / "spectra"),
            "fallback": [],
        },
        "fasta-files": {
            "name": "RNA sequences",
            "file_types": ["fasta", "fa"],
            "dir": str(directory / "fasta"),
            "fallback": [],
        },
        MAPPING: {
            "name": MAPPING_NAME,
            "file_types": ["csv"],
            "dir": str(directory / "mapping"),
            "fallback": [],
            "manual_import_sources": (ui.PRIDE,),
        },
    }


@pytest.fixture
def listing():
    names = ["run.mzML", "database.fasta", "sample-metadata.csv"]
    return SimpleNamespace(
        title="Public RNA dataset",
        url="https://www.ebi.ac.uk/pride/archive/projects/PXD000001",
        notes=[],
        files=[RemoteFile(name, 10, f"https://example.org/{name}") for name in names],
    )


@pytest.fixture
def import_app(tmp_path, monkeypatch, listing):
    resolver = Mock(return_value=listing)
    starter = Mock()
    monkeypatch.setattr(ui, "resolve", resolver)
    monkeypatch.setattr(job, "start_import", starter)

    def create(source=ui.PRIDE, *, cached=False):
        workspace = tmp_path / "workflow"
        workspace.mkdir(exist_ok=True)
        targets = targets_for(workspace)
        app = AppTest.from_string(HARNESS)
        app.session_state.workflow_dir = str(workspace)
        app.session_state.import_source = source
        app.session_state.import_targets = targets
        app.session_state.location = "local"
        app.session_state.settings = {"data_import": {"enabled": True}}
        if cached:
            # A previous app version routed this CSV before the PRIDE policy.
            app.session_state[f"{ui.LISTING_KEY}-{source}-{workspace.name}"] = {
                "title": listing.title,
                "url": listing.url,
                "notes": [],
                "items": job.build_items(listing.files, targets),
            }
        app.run(timeout=30)
        assert not app.exception
        return app, starter

    return create


def look_up(app, source):
    label = "PRIDE accession" if source == ui.PRIDE else "Link"
    next(item for item in app.text_input if item.label == label).set_value(
        "PXD000001" if source == ui.PRIDE else "https://example.org/dataset.zip"
    )
    next(button for button in app.button if button.label == "Look up").click().run(
        timeout=30
    )
    assert not app.exception


def selection_table(app):
    return next(
        element.value
        for element in app.dataframe
        if {"Import", "File", "Goes to"} <= set(element.value.columns)
    )


@pytest.mark.parametrize(
    "cached", [False, True], ids=["fresh-lookup", "cached-listing"]
)
def test_pride_csv_starts_skipped_while_spectra_and_fasta_remain_selected(
    import_app, cached
):
    app, starter = import_app(cached=cached)
    if not cached:
        look_up(app, ui.PRIDE)
    table = selection_table(app).set_index("File")
    assert not bool(table.loc["sample-metadata.csv", "Import"])
    assert table.loc["sample-metadata.csv", "Goes to"] == ui.SKIP
    assert bool(table.loc["run.mzML", "Import"])
    assert bool(table.loc["database.fasta", "Import"])
    app.button(key=f"{ui.LISTING_KEY}-{ui.PRIDE}-start").click().run(timeout=30)
    assert not app.exception
    chosen = starter.call_args.args[1]
    assert {item["name"] for item in chosen} == {"run.mzML", "database.fasta"}
    assert starter.call_args.args[2][MAPPING]["auto_route"] is False


def test_pride_csv_can_be_explicitly_selected_as_mapping(import_app):
    app, starter = import_app()
    look_up(app, ui.PRIDE)
    table = selection_table(app)
    row = int(table.index[table["File"] == "sample-metadata.csv"][0])
    edits = {
        "edited_rows": {row: {"Import": True, "Goes to": MAPPING_NAME}},
        "added_rows": [],
        "deleted_rows": [],
    }
    editor_id = next(
        element.proto.id
        for element in app.dataframe
        if "Goes to" in element.value.columns
    )
    # AppTest has no data-editor edit API. Send the same JSON widget state as
    # the browser, together with the user's explicit Import button click.
    app.button(key=f"{ui.LISTING_KEY}-{ui.PRIDE}-start").click()
    states = app._tree.get_widget_states()
    states.widgets.add(id=editor_id, string_value=json.dumps(edits))
    app._run(states, timeout=30)
    assert not app.exception
    chosen = starter.call_args.args[1]
    assert (
        next(item for item in chosen if item["name"].endswith(".csv"))["target"]
        == MAPPING
    )


def test_cloud_csv_automatic_routing_is_unchanged(import_app):
    app, starter = import_app(ui.LINK)
    look_up(app, ui.LINK)
    table = selection_table(app).set_index("File")
    assert bool(table.loc["sample-metadata.csv", "Import"])
    assert table.loc["sample-metadata.csv", "Goes to"] == MAPPING_NAME
    app.button(key=f"{ui.LISTING_KEY}-{ui.LINK}-start").click().run(timeout=30)
    assert not app.exception
    assert (
        next(
            item for item in starter.call_args.args[1] if item["name"].endswith(".csv")
        )["target"]
        == MAPPING
    )


def test_source_adjustment_does_not_mutate_metadata_or_restore_disabled_routing(
    tmp_path,
):
    targets = targets_for(tmp_path)
    targets["fasta-files"]["auto_route"] = False
    original = deepcopy(targets)
    pride = ui._targets_for_source(targets, ui.PRIDE)
    cloud = ui._targets_for_source(targets, ui.LINK)
    assert targets == original
    assert pride[MAPPING]["auto_route"] is False
    assert cloud[MAPPING].get("auto_route", True) is True
    assert pride["fasta-files"]["auto_route"] is False
    assert cloud["fasta-files"]["auto_route"] is False
    assert pride[MAPPING]["file_types"] == ["csv"]
    assert pride[MAPPING]["name"] == MAPPING_NAME


@pytest.mark.parametrize("name", ["metadata.csv", "metadata.CSV", "metadata.csv.gz"])
def test_native_job_routing_honors_manual_target_without_removing_its_file_type(
    tmp_path, name
):
    targets = targets_for(tmp_path)
    targets[MAPPING]["auto_route"] = False
    assert job.target_for(name, targets) == ""
    assert targets[MAPPING]["file_types"] == ["csv"]
    assert job.target_for("run.mzML", targets) == "mzML-files"


@pytest.fixture
def fake_fetch(monkeypatch):
    payloads = {}

    def fetch(_self, _item, remote, destination, base=0):
        del base
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payloads[remote["path"]])
        return destination

    monkeypatch.setattr(job._Job, "fetch", fetch)
    return payloads


@pytest.mark.parametrize("source", [ui.PRIDE, ui.LINK], ids=["pride", "cloud"])
def test_whole_archive_routing_skips_pride_csv_and_keeps_scientific_inputs(
    tmp_path, fake_fetch, source
):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("dataset/run.mzML", b"spectra")
        archive.writestr("dataset/database.fasta", b">rna\nACGU\n")
        archive.writestr("dataset/sample-metadata.csv", b"sample,condition\n1,A\n")
    fake_fetch["dataset.zip"] = buffer.getvalue()
    targets = ui._targets_for_source(targets_for(tmp_path), source)
    items = job.build_items(
        [
            RemoteFile(
                "dataset.zip", len(buffer.getvalue()), "https://example.org/dataset.zip"
            )
        ],
        targets,
    )
    status = job.run_import(
        str(tmp_path / "job"), items, targets, guard=False, parallel=1
    )
    assert status["state"] == "done", status
    assert not status["errors"]
    assert (Path(targets["mzML-files"]["dir"]) / "run.mzML").read_bytes() == b"spectra"
    assert (
        Path(targets["fasta-files"]["dir"]) / "database.fasta"
    ).read_bytes() == b">rna\nACGU\n"
    mapping = Path(targets[MAPPING]["dir"]) / "sample-metadata.csv"
    assert mapping.exists() is (source == ui.LINK)
    assert any("sample-metadata.csv" in name for name in status["skipped"]) is (
        source == ui.PRIDE
    )


def test_explicit_mapping_download_is_allowed_when_automatic_routing_is_disabled(
    tmp_path, fake_fetch
):
    fake_fetch["mapping.csv"] = b"id,chebi\ncustom,CHEBI:123\n"
    targets = targets_for(tmp_path)
    targets[MAPPING]["auto_route"] = False
    items = job.build_items(
        [RemoteFile("mapping.csv", 10, "https://example.org/mapping.csv")], targets
    )
    items[0]["target"] = MAPPING
    status = job.run_import(
        str(tmp_path / "job"), items, targets, guard=False, parallel=1
    )
    assert status["state"] == "done", status
    assert not status["errors"]
    assert (Path(targets[MAPPING]["dir"]) / "mapping.csv").read_bytes() == fake_fetch[
        "mapping.csv"
    ]


def test_real_naseweis_upload_registers_pride_only_manual_mapping_policy(
    naseweis_app,  # noqa: F811
    monkeypatch,
    listing,
):
    monkeypatch.setattr(ui, "resolve", Mock(return_value=listing))
    naseweis_app.session_state.settings["data_import"] = {"enabled": True}
    naseweis_app.run(timeout=30)
    naseweis_app.switch_page("content/topp_workflow_file_upload.py").run(timeout=30)
    assert not naseweis_app.exception
    next(
        item for item in naseweis_app.text_input if item.label == "PRIDE accession"
    ).set_value("PXD000001")
    naseweis_app.button(
        key=f"FormSubmitter:{ui.LISTING_KEY}-{ui.PRIDE}-form-Look up"
    ).click().run(timeout=30)
    assert not naseweis_app.exception
    table = selection_table(naseweis_app).set_index("File")
    assert not bool(table.loc["sample-metadata.csv", "Import"])
    assert table.loc["sample-metadata.csv", "Goes to"] == ui.SKIP
    assert bool(table.loc["run.mzML", "Import"])


def test_configure_mapping_browse_preserves_pride_manual_choice(
    naseweis_app,  # noqa: F811
    topp_schemas,  # noqa: F811
    monkeypatch,
    listing,
):
    monkeypatch.setattr(ui, "resolve", Mock(return_value=listing))
    naseweis_app.session_state.settings["data_import"] = {"enabled": True}
    naseweis_app.run(timeout=30)
    naseweis_app.switch_page("content/naseweis.py").run(timeout=30)
    naseweis_app.toggle(key="advanced").set_value(True).run(timeout=30)
    assert not naseweis_app.exception
    workspace = Path(naseweis_app.session_state.workspace) / "naseweis"
    # Seed a listing made before the policy change. Opening the real Configure
    # dialog must apply its mapping upload metadata to this cached CSV too.
    naseweis_app.session_state[f"{ui.LISTING_KEY}-{ui.PRIDE}-{workspace.name}"] = {
        "title": listing.title,
        "url": listing.url,
        "notes": [],
        "items": job.build_items(listing.files, targets_for(workspace)),
    }
    naseweis_app.button(key=f"_browse-{MAPPING}").click().run(timeout=30)
    assert not naseweis_app.exception
    table = selection_table(naseweis_app).set_index("File")
    assert not bool(table.loc["sample-metadata.csv", "Import"])
    assert table.loc["sample-metadata.csv", "Goes to"] == ui.SKIP
