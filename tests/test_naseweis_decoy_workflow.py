"""Database preparation contracts independent of the NASE search workflow."""

import json
from copy import deepcopy

import pytest
from test_naseweis_workflow import workflow  # noqa: F401


@pytest.fixture
def builder(workflow):  # noqa: F811
    from src.DecoyWorkflow import DecoyWorkflow

    builder = object.__new__(DecoyWorkflow)
    builder.__dict__.update(workflow.__dict__)
    builder.params = {
        "database-files": [workflow.params["fasta-files"]],
        "DecoyDatabase": {"type": "RNA", "method": "shuffle"},
    }
    return builder


@pytest.mark.parametrize(
    "selection", [None, [], "rna.fasta", [None], ["missing"], ["."]]
)
def test_invalid_input_cannot_start_builder(builder, selection):
    builder.params["database-files"] = selection
    assert builder.execution() is False
    builder.executor.run_topp.assert_not_called()


@pytest.mark.parametrize(
    "changes,reason",
    [
        ({"DecoyDatabase": {"type": "RNA", "entrapment_ratio": 1}}, "require protein"),
        ({"neighbor-files": "missing.fasta"}, "require protein"),
        (
            {"DecoyDatabase": {"type": "protein"}, "neighbor-files": "missing.fasta"},
            "no longer exists",
        ),
        (
            {"DecoyDatabase": {"method": "reverse", "shuffle_decoy_ratio": 2}},
            "shuffle method",
        ),
        ({"DecoyDatabase": {"missing_option": True}}, "Unsupported DecoyDatabase"),
    ],
)
def test_invalid_builder_modes_fail_without_rewriting_settings(
    builder, changes, reason
):
    original = builder.parameter_manager.params_file.read_bytes()
    builder.params.update(changes)
    assert builder.execution() is False
    builder.executor.run_topp.assert_not_called()
    assert builder.parameter_manager.params_file.read_bytes() == original
    assert reason in (builder.workflow_dir / "logs/minimal.log").read_text()


def test_missing_builder_parameters_is_actionable(builder):
    builder.parameter_manager.create_ini.return_value = False
    assert builder.execution() is False
    builder.executor.run_topp.assert_not_called()
    assert "not installed" in (builder.workflow_dir / "logs/minimal.log").read_text()


@pytest.mark.parametrize("success", [True, False])
def test_builder_uses_collected_inputs_and_propagates_failure(builder, success):
    builder.executor.run_topp.return_value = success
    original = deepcopy(builder.params)
    assert builder.execution() is success
    assert builder.params == original
    name, io = builder.executor.run_topp.call_args.args
    assert name == "DecoyDatabase"
    assert io["in"] == [builder.params["database-files"]]
    assert io["out"] == [str(builder.results_dir / "database/database.fasta")]
    snapshot = json.loads(builder.parameter_manager.params_file.read_text())
    assert "future_flag" in snapshot["_flag_params"]["DecoyDatabase"]


def test_neighbor_outputs_and_workspace_io_override_imported_paths(builder):
    relevant = builder.params["database-files"][0]
    builder.params.update(
        {
            "neighbor-files": relevant,
            "DecoyDatabase": {
                "type": "protein",
                "in": ["elsewhere"],
                "out": "/elsewhere/db",
            },
        }
    )
    assert builder.execution() is True
    io = builder.executor.run_topp.call_args.args[1]
    assert io["NeighborSearch:in_relevant_proteins"] == [relevant]
    assert io["NeighborSearch:out_neighbor"] == [
        str(builder.results_dir / "database/neighbors.fasta")
    ]
    assert io["NeighborSearch:out_relevant"] == [
        str(builder.results_dir / "database/relevant-targets-and-decoys.fasta")
    ]
    snapshot = builder.parameter_manager.read_parameters_strict()
    assert not {"in", "out"} & snapshot["DecoyDatabase"].keys()
