"""Imported defaults must not redirect the files owned by a workflow run."""

from copy import deepcopy
from unittest.mock import Mock

import pytest
from test_naseweis_decoy_workflow import builder  # noqa: F401
from test_naseweis_workflow import workflow  # noqa: F401

from src.naseweis_tools import (
    DECOY_INPUTS,
    DECOY_OUTPUTS,
    SEARCH_EXPORTS,
    SEARCH_INPUTS,
)

SEARCH_FILES = SEARCH_INPUTS | set(SEARCH_EXPORTS) | {"out"}
DECOY_FILES = DECOY_INPUTS | DECOY_OUTPUTS


def recording_executor(subject):
    """Run the real command builder, stopping only at subprocess execution."""
    from src.workflow.CommandExecutor import CommandExecutor

    executor = CommandExecutor(
        subject.workflow_dir, subject.logger, subject.parameter_manager
    )
    executor._max_threads = 1
    executor.run_command = Mock(return_value=True)
    subject.executor = executor
    return executor


def imported_files(keys, layer):
    return {key: f"/imported-{layer}/{key.replace(':', '-')}.wrong" for key in keys}


def assert_workspace_files(command, file_keys, expected):
    assert not any(argument.startswith("/imported-") for argument in command)
    for key in file_keys:
        flag = f"-{key}"
        assert command.count(flag) == int(key in expected), (key, command)
        if key in expected:
            assert command[command.index(flag) + 1] == expected[key]


@pytest.mark.parametrize("cached_digest", [False, True])
def test_search_commands_ignore_file_defaults_and_keep_scientific_options(
    workflow,  # noqa: F811
    cached_digest,
):
    # TOPP declares log as a string, but it can still write outside the run.
    search_paths = SEARCH_FILES | {"log"}
    decoy_paths = DECOY_FILES | {"log"}
    workflow.params.update(
        {
            "estimate-fdr": False,
            "search-exports": ["id_out", "out_preprocessed"],
            "_defaults": {
                "NucleicAcidSearchEngine": imported_files(search_paths, "defaults")
                | {"scoring:method": "mvh"},
                "DecoyDatabase": imported_files(decoy_paths, "defaults") | {"seed": 42},
            },
            "NucleicAcidSearchEngine": imported_files(search_paths, "overrides")
            | {"report:top_hits": 2},
            "DecoyDatabase": imported_files(decoy_paths, "overrides"),
        }
    )
    if cached_digest:
        digest = workflow.workflow_dir / "cached.oms"
        digest.touch()
        workflow.params.update(
            {"database-source": "Cached digest", "digest-files": str(digest)}
        )
    original = workflow.params
    before = deepcopy(original)
    executor = recording_executor(workflow)
    assert workflow.execution() is True
    assert original == before
    commands = [call.args[0] for call in executor.run_command.call_args_list]
    assert [command[0] for command in commands] == (
        [] if cached_digest else ["DecoyDatabase"]
    ) + ["NucleicAcidSearchEngine"]
    database = (
        workflow.params["digest-files"]
        if cached_digest
        else str(workflow.results_dir / "database/target-decoy.fasta")
    )
    assert_workspace_files(
        commands[-1],
        search_paths,
        {
            "in": workflow.params["mzML-files"],
            "digest" if cached_digest else "database": database,
            "out": str(workflow.results_dir / "mztab_results/search.mzTab"),
            "id_out": str(workflow.results_dir / "idxml_results/search.idXML"),
            "out_preprocessed": str(
                workflow.results_dir / "exports/preprocessed-MS2.mzML"
            ),
        },
    )
    assert commands[-1][commands[-1].index("-scoring:method") + 1] == "mvh"
    assert commands[-1][commands[-1].index("-report:top_hits") + 1] == "2"
    if not cached_digest:
        assert_workspace_files(
            commands[0],
            decoy_paths,
            {"in": original["fasta-files"], "out": database},
        )
        assert commands[0][commands[0].index("-seed") + 1] == "42"
    snapshot = workflow.parameter_manager.read_parameters_strict()
    for tool, keys in [
        ("NucleicAcidSearchEngine", search_paths),
        ("DecoyDatabase", decoy_paths),
    ]:
        assert not keys & snapshot[tool].keys()
        assert not keys & snapshot["_defaults"][tool].keys()


def test_builder_commands_ignore_file_defaults_and_preserve_neighbor_outputs(builder):  # noqa: F811
    # Native log paths are workspace-owned even though their schema type is string.
    owned_paths = DECOY_FILES | {"log"}
    builder.params.update(
        {
            "neighbor-files": builder.params["database-files"][0],
            "_defaults": {
                "DecoyDatabase": imported_files(owned_paths, "defaults") | {"seed": 23},
                "unrelated": {
                    "in": "keep-this-setting",
                    "log": "keep-this-log-setting",
                },
            },
            "DecoyDatabase": imported_files(owned_paths, "overrides")
            | {"type": "protein", "method": "reverse"},
        }
    )
    original = deepcopy(builder.params)
    executor = recording_executor(builder)
    assert builder.execution() is True
    command = executor.run_command.call_args.args[0]
    assert_workspace_files(
        command,
        owned_paths,
        {
            "in": builder.params["database-files"][0],
            "out": str(builder.results_dir / "database/database.fasta"),
            "NeighborSearch:in_relevant_proteins": builder.params["neighbor-files"],
            "NeighborSearch:out_neighbor": str(
                builder.results_dir / "database/neighbors.fasta"
            ),
            "NeighborSearch:out_relevant": str(
                builder.results_dir / "database/relevant-targets-and-decoys.fasta"
            ),
        },
    )
    assert command[command.index("-seed") + 1] == "23"
    assert command[command.index("-method") + 1] == "reverse"
    assert builder.params == original
    snapshot = builder.parameter_manager.read_parameters_strict()
    assert not owned_paths & snapshot["DecoyDatabase"].keys()
    assert not owned_paths & snapshot["_defaults"]["DecoyDatabase"].keys()
    assert snapshot["_defaults"]["unrelated"] == {
        "in": "keep-this-setting",
        "log": "keep-this-log-setting",
    }
