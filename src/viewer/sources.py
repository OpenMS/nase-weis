"""What a Viewer dataset holds.

A dataset is a saved run (``<run id>``), shown whole as its experiment, or one
acquisition of it (``<run id>/<sample id>``). ``resolve`` points a copy of the
workflow at it, the way the result sections expect.
"""
from __future__ import annotations

from copy import copy
from dataclasses import dataclass, field
from pathlib import Path

SEP = "/"
EXPERIMENT, ACQUISITION = "experiment", "acquisition"


def member_id(run_id: str, sample_id: str) -> str:
    return f"{run_id}{SEP}{sample_id}"


def split(dataset: str) -> tuple[str, str | None]:
    run_id, _, sample_id = dataset.partition(SEP)
    return run_id, sample_id or None


def manifest(wf, run_id: str) -> dict | None:
    """A run's experiment; None for a run of one file without one (or an
    unreadable manifest, which the result sections report)."""
    from src.naseweis_experiment import load_manifest

    try:
        return load_manifest(wf.runs.results_dir(run_id))
    except (OSError, ValueError):
        return None


def kind_of(wf, dataset: str) -> str:
    run_id, sample_id = split(dataset)
    if sample_id is None and manifest(wf, run_id) is not None:
        return EXPERIMENT
    return ACQUISITION


@dataclass
class Source:
    """A dataset with the workflow pointed at it."""

    wf: object
    run_id: str
    sample: dict | None
    status: str
    cache_dir: Path
    notes: list = field(default_factory=list)

    @property
    def dataset(self) -> str:
        return member_id(self.run_id, self.sample["id"]) if self.sample else self.run_id

    def mzml(self) -> Path | None:
        """The searched acquisition; None (with the reason in ``notes``) when
        it is not there as it was searched."""
        import streamlit as st

        from src.naseweis_viewers import read_run_snapshot, saved_acquisition_revision

        online = st.session_state.get("settings", {}).get("online_deployment", False)
        try:
            params = read_run_snapshot(self.wf._result_snapshot())
            return Path(saved_acquisition_revision(
                params, self.wf.workflow_dir, online, verify_content=False)[0])
        except (OSError, ValueError) as error:
            if str(error) not in self.notes:
                self.notes.append(str(error))
            return None

    def lfq_table(self) -> Path | None:
        from src.naseweis_lfq import LFQ_TABLE

        path = Path(self.wf.results_dir, LFQ_TABLE)
        return path if path.is_file() else None

    def lfq_requested(self) -> bool:
        """The run was started with quantification, whether or not it finished:
        its section then says why it has no features instead of offering to add it."""
        import json

        from src.naseweis_tools import enabled

        try:
            saved = json.loads(Path(self.wf._result_snapshot()).read_text())
        except (OSError, ValueError):
            return False
        return isinstance(saved, dict) and enabled(saved.get("enable-lfq", False))


def resolve(wf, dataset: str) -> Source | None:
    """``dataset`` as a Source; None when the run is gone."""
    from src.naseweis_experiment import sample_results, sample_snapshot

    run_id, sample_id = split(dataset)
    _, statuses = wf.run_statuses()
    if run_id not in statuses:
        return None
    child = copy(wf)
    child.file_manager = copy(wf.file_manager)
    child._use_run(run_id)
    child._result_run_status = statuses[run_id]
    sample = None
    if sample_id is not None:
        found = manifest(wf, run_id)
        sample = next((s for s in (found or {}).get("samples", []) if s["id"] == sample_id), None)
        if sample is None:
            return None
        root = wf.runs.results_dir(run_id)
        try:
            child.results_dir = sample_results(root, sample_id)
            child._sample_snapshot = sample_snapshot(root, sample_id)
        except (OSError, ValueError):
            return None
        child._sample_log_dir = child._sample_snapshot.parent / "logs"
    cache = Path(wf.workflow_dir, "insight_cache", run_id, sample_id or "run")
    return Source(child, run_id, sample, statuses[run_id], cache)
