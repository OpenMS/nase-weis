"""Application contracts for inputs, parameters, failures and example loading.

TOPP execution is replaced here so every failure branch is deterministic. The
integration suite separately checks the real tools and reference results.
"""

import json
from copy import deepcopy
from pathlib import Path
from unittest.mock import Mock

import pytest


@pytest.fixture
def workflow(tmp_path, monkeypatch):
    # Import after collection: inherited tests temporarily replace Streamlit
    # while importing their own isolated ParameterManager instances.
    from src.Workflow import Workflow
    from src.workflow.FileManager import FileManager
    from src.workflow.Logger import Logger
    from src.workflow.ParameterManager import ParameterManager

    directory = tmp_path / "naseweis"
    directory.mkdir()
    spectra, database = tmp_path / "spectra.mzML", tmp_path / "rna.fasta"
    spectra.touch()
    database.write_text(">rna\nACGU\n")
    wf = object.__new__(Workflow)
    wf.workflow_dir = directory
    wf.results_dir = directory / "results"
    wf.results_dir.mkdir()
    wf.logger = Logger(directory)
    wf.file_manager = FileManager(directory, wf.results_dir)
    wf.parameter_manager = ParameterManager(directory)
    wf.params = {
        "mzML-files": str(spectra),
        "fasta-files": str(database),
        "add-decoys": True,
        "FDR_cutoff": 0.05,
    }
    wf.parameter_manager.write_parameters(wf.params)
    # A real ParamXML flag makes this fixture detect regressions in flag
    # discovery too, rather than relying only on the app's hardcoded list.
    for tool in ["DecoyDatabase", "NucleicAcidSearchEngine"]:
        original = (Path("tests/data/naseweis/nightly") / f"{tool}.ini").read_text(
            encoding="latin-1"
        )
        marker = '<ITEM name="future_flag" type="bool" value="false"/>'
        # Place an extra bool in instance 1, without altering its real IO fields.
        position = original.index(">", original.index('<NODE name="1"')) + 1
        (wf.parameter_manager.ini_dir / f"{tool}.ini").write_text(
            original[:position] + marker + original[position:], encoding="latin-1"
        )
    monkeypatch.setattr(wf.parameter_manager, "create_ini", Mock(return_value=True))
    monkeypatch.setattr("src.Workflow.shutil.which", lambda name: f"/tools/{name}")
    wf.executor = Mock()
    wf.executor.run_topp.return_value = True
    return wf


@pytest.mark.parametrize("key", ["mzML-files", "fasta-files"])
@pytest.mark.parametrize("selection", [None, "", [], ["a", "b"], "missing", "."])
def test_invalid_input_never_runs_a_tool_or_rewrites_parameters(
    workflow, key, selection
):
    before = workflow.parameter_manager.params_file.read_bytes()
    workflow.params[key] = selection
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert workflow.parameter_manager.params_file.read_bytes() == before
    assert "ERROR" in (workflow.workflow_dir / "logs/minimal.log").read_text()


@pytest.mark.parametrize("cutoff", [-0.01, 1.01, float("nan"), float("inf")])
def test_invalid_fdr_does_not_start_a_search(workflow, cutoff):
    before = workflow.parameter_manager.params_file.read_bytes()
    workflow.params["FDR_cutoff"] = cutoff
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert workflow.parameter_manager.params_file.read_bytes() == before


@pytest.mark.parametrize("missing", ["DecoyDatabase", "NucleicAcidSearchEngine"])
def test_missing_tool_is_an_actionable_failure(workflow, monkeypatch, missing):
    monkeypatch.setattr(
        "src.Workflow.shutil.which", lambda name: None if name == missing else name
    )
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert (
        "add them to PATH" in (workflow.workflow_dir / "logs/minimal.log").read_text()
    )


@pytest.mark.parametrize("tool", ["DecoyDatabase", "NucleicAcidSearchEngine"])
def test_unreadable_tool_parameters_do_not_run(workflow, tool):
    workflow.parameter_manager.create_ini.side_effect = lambda name: name != tool
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert (
        f"Cannot read {tool}"
        in (workflow.workflow_dir / "logs/minimal.log").read_text()
    )


@pytest.mark.parametrize("use_decoys", [True, False])
def test_explicit_scientific_settings_survive_and_fdr_is_reset(workflow, use_decoys):
    workflow.params.update(
        {
            "add-decoys": use_decoys,
            "FDR_cutoff": 0.01,
            "ms1_resolution": 100,
            "ms2_resolution": 100,
            "_defaults": {
                "NucleicAcidSearchEngine": {
                    "oligo:enzyme": "RNase_A",
                    "fragment:mass_tolerance": 50.0,
                    "fdr:decoy_pattern": "OLD_",
                    "fdr:remove_decoys": True,
                }
            },
            "NucleicAcidSearchEngine": {
                "oligo:enzyme": "no cleavage",
                "fragment:mass_tolerance": 0.02,
                "fragment:mass_tolerance_unit": "Da",
                "precursor:use_avg_mass": False,
                "precursor:include_unknown_charge": True,
                "modifications:variable": ["m5C", "Am"],
            },
            "DecoyDatabase": {
                "type": "protein",
                "only_decoy": False,
                "decoy_string": "WRONG_",
            },
        }
    )
    original = workflow.params
    untouched = deepcopy(original)
    assert workflow.execution() is True
    assert original == untouched, "execution mutated the caller's parameter object"
    saved = json.loads(workflow.parameter_manager.params_file.read_text())
    engine = saved["NucleicAcidSearchEngine"]
    assert engine["oligo:enzyme"] == "no cleavage"
    assert engine["fragment:mass_tolerance"] == 0.02
    assert engine["fragment:mass_tolerance_unit"] == "Da"
    assert engine["precursor:use_avg_mass"] is False
    assert engine["precursor:include_unknown_charge"] is True
    assert engine["modifications:variable"] == ["m5C", "Am"]
    assert engine["fdr:decoy_pattern"] == ("WRONG_" if use_decoys else "")
    assert engine["fdr:cutoff"] == (0.01 if use_decoys else 1.0)
    assert engine["fdr:remove_decoys"] is use_decoys
    assert saved["DecoyDatabase"]["type"] == "RNA"
    assert saved["DecoyDatabase"]["only_decoy"] is False
    assert saved["DecoyDatabase"]["decoy_string"] == "WRONG_"
    for tool in (["DecoyDatabase"] if use_decoys else []) + ["NucleicAcidSearchEngine"]:
        assert "future_flag" in saved["_flag_params"][tool]
    calls = workflow.executor.run_topp.call_args_list
    assert [call.args[0] for call in calls] == (
        ["DecoyDatabase"] if use_decoys else []
    ) + ["NucleicAcidSearchEngine"]
    search_io = calls[-1].args[1]
    assert search_io["in"] == [original["mzML-files"]]
    expected_database = (
        str(workflow.results_dir / "database/target-decoy.fasta")
        if use_decoys
        else original["fasta-files"]
    )
    assert search_io["database"] == [expected_database]


@pytest.mark.parametrize("outcomes,expected_calls", [([False], 1), ([True, False], 2)])
def test_failed_tool_stops_the_pipeline(workflow, outcomes, expected_calls):
    workflow.executor.run_topp.side_effect = outcomes
    assert workflow.execution() is False
    assert workflow.executor.run_topp.call_count == expected_calls


def test_example_loading_is_repeatable_and_preserves_user_files(workflow):
    # Existing files with the example names belong to the user too.
    for key, filename in [
        ("mzML-files", "NME1.mzML"),
        ("fasta-files", "NME1_Target.fasta"),
    ]:
        target = workflow.workflow_dir / "input-files" / key / filename
        target.parent.mkdir(parents=True)
        target.write_bytes(b"existing user data")
    stored = workflow.params | {"FDR_cutoff": 0.01, "custom-note": "keep me"}
    workflow.parameter_manager.write_parameters(stored)
    workflow.load_example_data()
    workflow.load_example_data()
    saved = workflow.parameter_manager.read_parameters_strict()
    assert saved["FDR_cutoff"] == 0.01 and saved["custom-note"] == "keep me"
    for key in ["mzML-files", "fasta-files"]:
        assert Path(saved[key]).read_bytes() == b"existing user data"


def test_example_loading_copies_both_inputs(workflow):
    workflow.load_example_data()
    saved = workflow.parameter_manager.read_parameters_strict()
    assert (
        Path(saved["mzML-files"]).read_bytes()
        == Path("example-data/mzML/NME1.mzML").read_bytes()
    )
    assert (
        Path(saved["fasta-files"]).read_bytes()
        == Path("example-data/fasta/NME1_Target.fasta").read_bytes()
    )


@pytest.mark.parametrize(
    "changes,reason",
    [
        ({"search-exports": ["imaginary"]}, "Unknown output"),
        (
            {
                "NucleicAcidSearchEngine": {
                    "precursor:correct_with_flash": True,
                    "precursor:use_avg_mass": True,
                }
            },
            "monoisotopic",
        ),
        ({"DecoyDatabase": {"decoy_string": " "}}, "marker must not be empty"),
        (
            {
                "add-decoys": False,
                "estimate-fdr": True,
                "NucleicAcidSearchEngine": {"fdr:decoy_pattern": ""},
            },
            "nonempty decoy",
        ),
        ({"DecoyDatabase": {"shuffle_decoy_ratio": 2}}, "1:1"),
        (
            {"NucleicAcidSearchEngine": {"old_removed_setting": True}},
            "Unsupported parameters",
        ),
        ({"DecoyDatabase": {"unknown_setting": True}}, "Unsupported parameters"),
    ],
)
def test_nightly_invalid_settings_fail_before_any_tool(workflow, changes, reason):
    before = workflow.parameter_manager.params_file.read_bytes()
    workflow.params.update(changes)
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert workflow.parameter_manager.params_file.read_bytes() == before
    assert reason in (workflow.workflow_dir / "logs/minimal.log").read_text()


def test_missing_export_capability_is_reported(workflow):
    path = workflow.parameter_manager.ini_dir / "NucleicAcidSearchEngine.ini"
    path.write_text(
        path.read_text(encoding="latin-1").replace(
            'name="db_out"', 'name="old_output"'
        ),
        encoding="latin-1",
    )
    workflow.params["search-exports"] = ["db_out"]
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert "pinned nightly" in (workflow.workflow_dir / "logs/minimal.log").read_text()


def test_cached_digest_ignores_fasta_and_does_not_require_decoy_tool(
    workflow, monkeypatch
):
    cached = workflow.workflow_dir / "cached.oms"
    cached.touch()
    workflow.params.update(
        {
            "database-source": "Cached digest",
            "digest-files": str(cached),
            "fasta-files": None,
        }
    )
    monkeypatch.setattr(
        "src.Workflow.shutil.which",
        lambda name: None if name == "DecoyDatabase" else name,
    )
    assert workflow.execution() is True
    assert workflow.executor.run_topp.call_count == 1
    io = workflow.executor.run_topp.call_args.args[1]
    assert io["digest"] == [str(cached)]
    assert "database" not in io


def test_vendor_directory_and_custom_mapping_reach_native_inputs(workflow):
    folder = workflow.workflow_dir / "sample.d"
    folder.mkdir()
    mapping = workflow.workflow_dir / "modifications.csv"
    mapping.write_text("mod,chebi_id\nAm,28177\n")
    workflow.params.update(
        {
            "mzML-files": str(folder),
            "add-decoys": False,
            "bedrmod-mapping-files": str(mapping),
        }
    )
    assert workflow.execution() is True
    io = workflow.executor.run_topp.call_args.args[1]
    assert io["in"] == [str(folder)]
    assert io["bedrmod_chebi_mapping"] == [str(mapping)]


def test_missing_modification_mapping_fails(workflow):
    workflow.params.update(
        {"add-decoys": False, "bedrmod-mapping-files": "missing.csv"}
    )
    assert workflow.execution() is False
    workflow.executor.run_topp.assert_not_called()
    assert (
        "mapping no longer exists"
        in (workflow.workflow_dir / "logs/minimal.log").read_text()
    )


def test_imported_file_parameters_cannot_replace_workspace_io(workflow):
    workflow.params.update(
        {
            "NucleicAcidSearchEngine": {
                "out": "/elsewhere/result",
                "in": "other.mzML",
                "digest": "other.oms",
            },
            "DecoyDatabase": {"in": ["other.fasta"], "out": "/elsewhere/db"},
        }
    )
    assert workflow.execution() is True
    saved = workflow.parameter_manager.read_parameters_strict()
    assert not {"in", "out", "digest"} & saved["NucleicAcidSearchEngine"].keys()
    assert not {"in", "out"} & saved["DecoyDatabase"].keys()
