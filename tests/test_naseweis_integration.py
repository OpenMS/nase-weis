"""Execute real TOPP searches on the bundled example, including the worker path."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from src.naseweis_results import coverage_html, match_table, read_mztab
from src.workflow.RunManager import RunManager

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(
    not all(
        shutil.which(tool) for tool in ["DecoyDatabase", "NucleicAcidSearchEngine"]
    ),
    reason="Real OpenMS TOPP executables are required for this integration test.",
)


@pytest.mark.parametrize("use_decoys", [True, False])
def test_example_search_in_fresh_worker(tmp_path, use_decoys, request):
    workflow = tmp_path / "naseweis"
    workflow.mkdir()
    # Include stale FDR and average-mass settings from an older workspace.
    params = {
        "mzML-files": str(ROOT / "example-data/mzML/NME1.mzML"),
        "fasta-files": str(ROOT / "example-data/fasta/NME1_Target.fasta"),
        "add-decoys": use_decoys,
        "FDR_cutoff": 0.05,
        "max_threads": 1 if request.config.getoption("--naseweis-baseline") else 2,
        "_defaults": {"NucleicAcidSearchEngine": {"fdr:decoy_pattern": "STALE_"}},
        "NucleicAcidSearchEngine": {"precursor:use_avg_mass": False},
    }
    param_file = workflow / "params.json"
    param_file.write_text(json.dumps(params))
    runs = RunManager(workflow)
    run_id = runs.new_run("NME1 integration", param_file)
    script = (
        "from src.workflow.tasks import execute_workflow; "
        "import sys; "
        "execute_workflow(sys.argv[1], 'Workflow', 'src.Workflow', sys.argv[2])"
    )
    completed = subprocess.run(
        [sys.executable, "-c", script, str(workflow), run_id],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, (
        completed.stdout + completed.stderr + (workflow / "logs/all.log").read_text()
    )
    results_dir = runs.results_dir(run_id)
    tables = read_mztab(results_dir / "mztab_results/search.mzTab")
    if request.config.getoption("--naseweis-baseline"):
        from naseweis_reference import REFERENCE, assert_nme1_reference

        assert_nme1_reference(tables, use_decoys, reference=REFERENCE / "nightly")
    matches = match_table(tables)
    assert len(matches) > 0
    assert "NME1" in tables["NUC"]["accession"].values
    assert "sequence coverage" in coverage_html(tables, "NME1")
    assert (results_dir / "idxml_results/search.idXML").stat().st_size > 100
    assert "WORKFLOW FINISHED" in (workflow / "logs/minimal.log").read_text()
    # A worker only updates its run snapshot, not the user's current configuration.
    assert json.loads(param_file.read_text()) == params
    snapshot = json.loads(runs.params_file(run_id).read_text())
    assert snapshot["NucleicAcidSearchEngine"]["fdr:decoy_pattern"] == (
        "DECOY_" if use_decoys else ""
    )
    assert snapshot["NucleicAcidSearchEngine"]["precursor:use_avg_mass"] is False
    if use_decoys:
        assert matches["q-value (%)"].between(0, 5).all()
        assert (results_dir / "database/target-decoy.fasta").is_file()
        assert not matches["Accession"].str.contains("DECOY_").any()
    else:
        assert "q-value (%)" not in matches
        assert not (results_dir / "database").exists()


def test_missing_input_is_a_failed_worker_job(tmp_path):
    workflow = tmp_path / "naseweis"
    workflow.mkdir()
    (workflow / "params.json").write_text('{"mzML-files": "missing.mzML"}')
    from src.workflow.tasks import WorkflowExecutionError, execute_workflow

    with pytest.raises(WorkflowExecutionError):
        execute_workflow(str(workflow), "Workflow", "src.Workflow")
    log = (workflow / "logs/minimal.log").read_text()
    assert "ERROR" in log
    assert "WORKFLOW FINISHED" not in log
