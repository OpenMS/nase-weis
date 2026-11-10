"""Real nightly feature checks: exports, FLASH/MVH, cached digests and decoy tools."""

import json
import shutil
import sqlite3
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def nightly_runtime():
    if not shutil.which("NucleicAcidSearchEngine"):
        pytest.skip("The pinned nightly TOPP executables are not installed.")
    result = subprocess.run(
        ["NucleicAcidSearchEngine", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    revision = json.loads((ROOT / "openms-runtime.json").read_text())["revision"]
    if f"Revision: {revision}" not in result.stdout + result.stderr:
        pytest.skip("This feature test requires the pinned OpenMS nightly runtime.")


def run_search(tmp_path, name, overrides=None):
    from src.workflow.RunManager import RunManager
    from src.workflow.tasks import execute_workflow

    workflow = tmp_path / name / "naseweis"
    workflow.mkdir(parents=True)
    params = {
        "mzML-files": str(ROOT / "example-data/mzML/NME1.mzML"),
        "fasta-files": str(ROOT / "example-data/fasta/NME1_Target.fasta"),
        "add-decoys": False,
        "estimate-fdr": False,
        "max_threads": 2,
    } | (overrides or {})
    source = workflow / "params.json"
    source.write_text(json.dumps(params))
    runs = RunManager(workflow)
    run_id = runs.new_run(name, source)
    try:
        execute_workflow(str(workflow), "Workflow", "src.Workflow", run_id)
    except Exception:
        print((workflow / "logs/all.log").read_text())
        raise
    assert json.loads(source.read_text()) == params
    assert runs.get(run_id)["status"] == "finished"
    return runs.results_dir(run_id), json.loads(runs.params_file(run_id).read_text())


@pytest.mark.parametrize(
    "parameters,export",
    [
        (
            {
                "scoring:method": "mvh",
                "scoring:single_peak_credit": "false",
                "scoring:use_mass_accuracy": True,
            },
            "id_out",
        ),
        (
            {
                "preprocessing:deisotoping": "flash",
                "preprocessing:remove_precursor_peak": True,
            },
            "out_preprocessed",
        ),
        ({"precursor:correct_with_flash": True}, "out_precursor_corrections"),
        (
            {
                "scoring:method": "mvh",
                "scoring:num_intensity_classes": 4,
                "scoring:tic_fraction": 0.95,
                "scoring:use_mass_accuracy": True,
                "precursor:correct_with_flash": True,
                "preprocessing:deisotoping": "flash",
                "preprocessing:filter_window_mower": True,
                "preprocessing:filter_nlargest": True,
            },
            "out_precursor_corrections",
        ),
        (
            {
                "report:top_hits": 3,
                "report:require_coverage": 0.2,
                "report:precursor_percentage": True,
            },
            "lfq_out",
        ),
    ],
)
def test_nightly_algorithms_run_with_their_exports(
    tmp_path, nightly_runtime, parameters, export
):
    from src.naseweis_results import match_table, read_mztab
    from src.naseweis_tools import SEARCH_EXPORTS

    directory, snapshot = run_search(
        tmp_path,
        "feature",
        {
            "NucleicAcidSearchEngine": parameters,
            "search-exports": ["id_out", export] if export != "id_out" else ["id_out"],
        },
    )
    assert len(match_table(read_mztab(directory / "mztab_results/search.mzTab"))) > 0
    assert (directory / SEARCH_EXPORTS[export][1]).stat().st_size > 0
    for key, expected in parameters.items():
        assert snapshot["NucleicAcidSearchEngine"][key] == expected


def test_all_native_exports_are_generated_and_readable(tmp_path, nightly_runtime):
    from src.naseweis_tools import SEARCH_EXPORTS

    mapping = tmp_path / "mapping.csv"
    mapping.write_text("mod,chebi_id,modomics_short_name\nAm,CHEBI:28177,Am\n")
    directory, _ = run_search(
        tmp_path,
        "exports",
        {
            "search-exports": list(SEARCH_EXPORTS),
            "bedrmod-mapping-files": str(mapping),
            "NucleicAcidSearchEngine": {
                "precursor:correct_with_flash": True,
                "report:precursor_percentage": True,
                "bedrmod:organism": "Saccharomyces cerevisiae",
                "bedrmod:assembly": "test-reference",
            },
        },
    )
    for label, relative, mime in SEARCH_EXPORTS.values():
        artifact = directory / relative
        assert artifact.is_file() and artifact.stat().st_size > 0, label
    ET.parse(directory / "idxml_results/search.idXML")
    for relative in ["exports/search.oms", "exports/digest.oms"]:
        with sqlite3.connect(directory / relative) as connection:
            assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
    assert not pd.read_csv(directory / "exports/lfq-targets.tsv", sep="\t").empty
    import pyopenms

    for relative in [
        "exports/theoretical-MS2.mzML",
        "exports/experimental-MS2.mzML",
        "exports/preprocessed-MS2.mzML",
    ]:
        experiment = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(str(directory / relative), experiment)
        assert experiment.size() > 0
    assert (
        len(pd.read_csv(directory / "exports/precursor-corrections.tsv", sep="\t"))
        == 2116
    )
    assert (
        "Saccharomyces cerevisiae"
        in (directory / "exports/modifications.bedrmod").read_text()
    )


def test_cached_digest_reuses_identical_assignments_without_a_fasta(
    tmp_path, nightly_runtime
):
    from src.naseweis_results import read_mztab

    first, _ = run_search(tmp_path, "fresh", {"search-exports": ["digest_out"]})
    cached, _ = run_search(
        tmp_path,
        "cached",
        {
            "database-source": "Cached digest",
            "digest-files": str(first / "exports/digest.oms"),
            "fasta-files": None,
        },
    )
    a = read_mztab(first / "mztab_results/search.mzTab")["OSM"]
    b = read_mztab(cached / "mztab_results/search.mzTab")["OSM"]
    columns = ["sequence", "spectra_ref", "charge", "search_engine_score[1]"]
    pd.testing.assert_frame_equal(
        a[columns].sort_values(columns).reset_index(drop=True),
        b[columns].sort_values(columns).reset_index(drop=True),
    )


@pytest.mark.parametrize("suffix,only_decoy", [(False, False), (True, True)])
def test_custom_markers_and_decoy_only_exports_keep_target_sequences(
    tmp_path, nightly_runtime, suffix, only_decoy
):
    from src.naseweis_results import match_table, read_mztab

    directory, snapshot = run_search(
        tmp_path,
        "custom-decoys",
        {
            "add-decoys": True,
            "estimate-fdr": True,
            "FDR_cutoff": 1.0,
            "DecoyDatabase": {
                "decoy_string": "_CUSTOM_DECOY" if suffix else "DECOY_CUSTOM_",
                "decoy_string_position": "suffix" if suffix else "prefix",
                "only_decoy": only_decoy,
            },
        },
    )
    assert snapshot["NucleicAcidSearchEngine"]["fdr:decoy_pattern"] == (
        "_CUSTOM_DECOY" if suffix else "DECOY_CUSTOM_"
    )
    assert snapshot["DecoyDatabase"]["only_decoy"] is only_decoy
    matches = match_table(read_mztab(directory / "mztab_results/search.mzTab"))
    assert len(matches) > 0
    assert not matches["Accession"].str.contains("CUSTOM_DECOY|DECOY_CUSTOM").any()
    if only_decoy:
        generated = (directory / "database/decoys.fasta").read_text()
        assert ">NME1\n" not in generated
        assert ">NME1\n" in (directory / "database/target-and-decoy.fasta").read_text()


def test_provided_target_decoy_database_does_not_generate_more_decoys(
    tmp_path, nightly_runtime
):
    from src.naseweis_results import match_table, read_mztab

    first, _ = run_search(
        tmp_path, "generated", {"add-decoys": True, "estimate-fdr": True}
    )
    second, _ = run_search(
        tmp_path,
        "provided",
        {
            "fasta-files": str(first / "database/target-decoy.fasta"),
            "estimate-fdr": True,
        },
    )
    a = match_table(read_mztab(first / "mztab_results/search.mzTab"))
    b = match_table(read_mztab(second / "mztab_results/search.mzTab"))
    pd.testing.assert_frame_equal(a, b)
    assert not (second / "database").exists()


@pytest.mark.parametrize(
    "tool_params,neighbor",
    [
        ({"type": "RNA", "method": "reverse", "only_decoy": True}, False),
        ({"type": "RNA", "method": "shuffle", "shuffle_decoy_ratio": 2}, False),
        ({"type": "protein", "method": "shuffle", "entrapment_ratio": 1}, False),
        ({"type": "protein", "method": "reverse"}, True),
    ],
)
def test_standalone_database_builder_modes(
    tmp_path, nightly_runtime, tool_params, neighbor
):
    from src.workflow.RunManager import RunManager
    from src.workflow.tasks import execute_workflow

    workflow = tmp_path / "decoy-database"
    workflow.mkdir()
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">test\n"
        + (
            "ACGUAUGCAUCGACUGAUGCUAGCUAGC"
            if tool_params["type"] == "RNA"
            else "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
        )
        + "\n"
    )
    params = {
        "database-files": [str(fasta)],
        "DecoyDatabase": tool_params,
        "max_threads": 1,
    }
    if neighbor:
        params["neighbor-files"] = str(fasta)
    source = workflow / "params.json"
    source.write_text(json.dumps(params))
    runs = RunManager(workflow)
    run_id = runs.new_run("database mode", source)
    try:
        execute_workflow(str(workflow), "DecoyWorkflow", "src.DecoyWorkflow", run_id)
    except Exception:
        print((workflow / "logs/all.log").read_text())
        raise
    directory = runs.results_dir(run_id) / "database"
    text = (directory / "database.fasta").read_text()
    expected = (
        1
        if tool_params.get("only_decoy")
        else 3
        if tool_params.get("shuffle_decoy_ratio") == 2
        else 4
        if tool_params.get("entrapment_ratio") or neighbor
        else 2
    )
    assert sum(line.startswith(">") for line in text.splitlines()) == expected
    if neighbor:
        assert (directory / "neighbors.fasta").is_file()
        assert (directory / "relevant-targets-and-decoys.fasta").is_file()
