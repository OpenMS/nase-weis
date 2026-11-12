"""Execute an immutable, explicitly designed cohort as isolated sample searches.

Malformed scientific metadata raises ValueError consistently for the app UI.
"""

import csv
import io
import json
import math
import os
import re
import uuid
from copy import copy, deepcopy
from decimal import Decimal, InvalidOperation
from pathlib import Path

from src.naseweis_tools import enabled

MANIFEST = "experiment.json"
_STATUSES = {"pending", "running", "finished", "failed"}
_DESIGN_NUMBERS = {
    "fraction",
    "replicate",
    "msstats_bioreplicate",
    "biological replicate",
    "biological_replicate",
    "characteristics[biological replicate]",
    "comment[biological replicate]",
    "comment[fraction identifier]",
}


def validate_design_source(content: str | bytes) -> None:
    """Reject invalid raw design numbers before the framework normalizes them.

    Only the app's biological-replicate and fraction columns are inspected.
    The existing OpenMS/SDRF parser remains responsible for table structure;
    omitted numeric columns retain its defaults, but present blank cells are
    invalid. Blank lines separate the two OpenMS design tables.
    """
    if isinstance(content, bytes):
        try:
            content = content.decode("utf-8-sig")
        except UnicodeError as error:
            raise ValueError("The experimental design must be UTF-8 text.") from error
    if not isinstance(content, str):
        raise ValueError("The experimental design must be text.")  # noqa: TRY004
    reader = csv.reader(
        io.StringIO(content.lstrip("\ufeff")), delimiter="\t", strict=True
    )
    columns = None
    try:
        for row in reader:
            if not row or (len(row) == 1 and not row[0].strip()):
                columns = None
                continue
            if columns is None:
                columns = [
                    (number, name.strip())
                    for number, name in enumerate(row)
                    if " ".join(name.strip().casefold().split()) in _DESIGN_NUMBERS
                ]
                continue
            for number, name in columns:
                value = row[number].strip() if number < len(row) else ""
                try:
                    numeric = Decimal(value)
                except InvalidOperation:
                    numeric = Decimal("NaN")
                if (
                    not numeric.is_finite()
                    or numeric <= 0
                    or numeric != numeric.to_integral_value()
                ):
                    raise ValueError(
                        f"Experimental design {name} on line {reader.line_num} must be a positive integer; received {value!r}."
                    )
    except csv.Error as error:
        raise ValueError(
            "The experimental design is not valid tab-separated text."
        ) from error


def _sample_directory(results_root: str | Path, sample_id: str) -> Path:
    if (
        not isinstance(sample_id, str)
        or not re.fullmatch(r"sample_[0-9]{3,}", sample_id)
        or not set(sample_id.removeprefix("sample_")) - {"0"}
    ):
        raise ValueError("Invalid controlled experiment sample ID.")
    try:
        directory = Path(results_root).resolve() / "samples" / sample_id
        resolved = directory.resolve()
    except (OSError, RuntimeError) as error:
        raise ValueError("Unsafe experiment sample path.") from error
    if resolved != directory:
        raise ValueError(
            "Unsafe experiment sample path: symbolic links are not allowed."
        )
    return directory


def _sample_artifact(results_root: str | Path, sample_id: str, name: str) -> Path:
    path = _sample_directory(results_root, sample_id) / name
    try:
        resolved = path.resolve()
    except (OSError, RuntimeError) as error:
        raise ValueError("Unsafe experiment artifact path.") from error
    if resolved != path:
        raise ValueError(
            "Unsafe experiment artifact path: symbolic links are not allowed."
        )
    return path


def sample_results(results_root: str | Path, sample_id: str) -> Path:
    """Derive a controlled sample result directory; never use manifest paths."""
    return _sample_artifact(results_root, sample_id, "results")


def sample_snapshot(results_root: str | Path, sample_id: str) -> Path:
    """Derive a controlled sample parameter snapshot; reject linked artifacts."""
    return _sample_artifact(results_root, sample_id, "params.json")


def _positive_integer(value, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive integer.")  # noqa: TRY004
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{label} must be a positive integer.") from error
    if not math.isfinite(number) or not number.is_integer() or number < 1:
        raise ValueError(f"{label} must be a positive integer.")
    return int(number)


def _selected_paths(params: dict) -> list[Path]:
    selected = params.get("mzML-files")
    if not isinstance(selected, list) or not selected:
        raise ValueError("Select one or more acquisitions for the experiment.")
    paths = []
    for value in selected:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                "Each selected acquisition must be one existing file path."
            )
        try:
            path = Path(value).resolve(strict=True)
        except OSError as error:
            raise ValueError(
                f"Selected acquisition is unavailable: {Path(value).name}. Select an existing acquisition on Configure."
            ) from error
        suffix = path.suffix.lower()
        if suffix not in {".mzml", ".raw", ".d"} or not (
            path.is_file() if suffix != ".d" else path.is_dir()
        ):
            raise ValueError("Select existing mzML/RAW files or Bruker .d directories.")
        paths.append(path)
    if len(set(paths)) != len(paths):
        raise ValueError("The experiment selects the same acquisition more than once.")
    if len({path.name.casefold() for path in paths}) != len(paths):
        raise ValueError("Selected acquisitions must have distinct file names.")
    return paths


def _design_samples(params: dict, paths: list[Path], lfq: bool) -> list[dict]:
    design = params.get("_sample_design")
    if design is None and "_sample_design" not in params and not lfq:
        design = [
            {"file": str(path), "condition": "", "replicate": number, "fraction": 1}
            for number, path in enumerate(paths, 1)
        ]
    if not isinstance(design, list) or len(design) != len(paths):
        raise ValueError(
            "Save an experimental design with exactly one row per selected acquisition."
        )
    by_path, by_name = (
        {str(path): path for path in paths},
        {path.name: path for path in paths},
    )

    def selected(reference):
        if not isinstance(reference, str) or not reference:
            raise ValueError("Each design row must identify a selected acquisition.")
        path = by_name.get(reference) or by_path.get(str(Path(reference).resolve()))
        if path is None:
            raise ValueError(
                "Experimental design rows must match the selected acquisitions exactly."
            )
        return path

    rows = {}
    for row in design:
        if not isinstance(row, dict):
            raise ValueError("Each experimental design row must be an object.")  # noqa: TRY004
        path = selected(row.get("path", row.get("file")))
        if "path" in row and "file" in row and selected(row["file"]) != path:
            raise ValueError("A design row has conflicting acquisition paths.")
        if path in rows:
            raise ValueError(
                "The experimental design lists an acquisition more than once."
            )
        condition = row.get("condition", "")
        if not isinstance(condition, str) or (lfq and not condition.strip()):
            raise ValueError("Assign a nonblank condition to every LFQ acquisition.")
        replicate = _positive_integer(
            row.get("replicate", len(rows) + 1) if not lfq else row.get("replicate"),
            "Biological replicate",
        )
        fraction = _positive_integer(
            row.get("fraction", 1) if not lfq else row.get("fraction"), "Fraction"
        )
        if lfq and fraction != 1:
            raise ValueError(
                "RNA LFQ currently requires unfractionated acquisitions (fraction 1)."
            )
        rows[path] = {
            "condition": condition,
            "replicate": replicate,
            "fraction": fraction,
        }
    return [
        {
            "id": f"sample_{number:03d}",
            "file": path.name,
            "path": str(path),
            **rows[path],
            "status": "pending",
        }
        for number, path in enumerate(paths, 1)
    ]


def preflight_experiment(params: dict) -> tuple[list[dict], bool]:
    """Validate selection, immutable design and comparison choices before work."""
    paths = _selected_paths(params)
    lfq = enabled(params.get("enable-lfq", False))
    if lfq and len(paths) < 2:
        raise ValueError("Relative LFQ requires at least two selected acquisitions.")
    samples = _design_samples(params, paths, lfq)
    if lfq:
        if params.get("lfq-measure", "intensity") not in ("intensity", "model_area"):
            raise ValueError("Choose LFQ measure intensity or model_area.")
        if params.get("lfq-normalization", "total") not in ("total", "none"):
            raise ValueError("Choose LFQ normalization total or none.")
        if params.get("lfq-measure", "intensity") == "model_area":
            defaults = params.get("_defaults", {})
            if not isinstance(defaults, dict):
                raise ValueError("Invalid saved feature extraction defaults.")
            model_defaults = defaults.get("FeatureFinderMetaboIdent", {})
            model_overrides = params.get("FeatureFinderMetaboIdent", {})
            if not isinstance(model_defaults, dict) or not isinstance(
                model_overrides, dict
            ):
                raise ValueError("Invalid saved feature extraction settings.")
            model = (model_defaults | model_overrides).get("model:type", "none")
            if model not in ("symmetric", "asymmetric"):
                raise ValueError(
                    "Fitted model_area requires a symmetric or asymmetric elution model."
                )
    return samples, lfq


def _atomic_manifest(root: Path, manifest: dict) -> None:
    temporary = root / f".{MANIFEST}.{uuid.uuid4().hex}.tmp"
    try:
        temporary.write_text(
            json.dumps(manifest, indent=2, allow_nan=False), encoding="utf-8"
        )
        os.replace(temporary, root / MANIFEST)
    finally:
        temporary.unlink(missing_ok=True)


def load_manifest(results_root: str | Path) -> dict | None:
    """Read a cohort manifest while validating its controlled artifact layout."""
    root = Path(results_root)
    path = root / MANIFEST
    if path.is_symlink():
        raise ValueError("This experiment's manifest cannot be a symbolic link.")
    if not path.exists():
        return None
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ValueError(
            "This experiment's manifest is unavailable or invalid."
        ) from error
    if (
        not isinstance(manifest, dict)
        or type(manifest.get("schema")) is not int
        or manifest["schema"] != 1
        or not isinstance(manifest.get("samples"), list)
        or not manifest["samples"]
        or type(manifest.get("lfq_enabled")) is not bool
        or not isinstance(manifest.get("status"), str)
        or manifest["status"] not in _STATUSES
        or not isinstance(manifest.get("comparison_status"), str)
        or manifest["comparison_status"]
        not in {"not_requested", "pending", "finished", "failed"}
    ):
        raise ValueError("This experiment's manifest has an invalid schema or status.")
    if manifest["lfq_enabled"] != (manifest["comparison_status"] != "not_requested"):
        raise ValueError(
            "This experiment's LFQ and comparison statuses are inconsistent."
        )
    seen_ids, seen_files = set(), set()
    required = {"id", "file", "path", "condition", "replicate", "fraction", "status"}
    for sample in manifest["samples"]:
        if (
            not isinstance(sample, dict)
            or not required <= sample.keys()
            or sample.keys() - required - {"error"}
        ):
            raise ValueError("This experiment's sample manifest is invalid.")
        sample_results(root, sample["id"])
        sample_snapshot(root, sample["id"])
        filename, acquisition = sample["file"], sample["path"]
        if (
            not isinstance(filename, str)
            or not filename
            or re.search(r"[\\/]", filename)
            or not isinstance(acquisition, str)
            or not acquisition
            or re.split(r"[\\/]", acquisition)[-1] != filename
            or not isinstance(sample["condition"], str)
            or not isinstance(sample["status"], str)
            or sample["status"] not in _STATUSES
            or sample["id"] in seen_ids
            or filename.casefold() in seen_files
        ):
            raise ValueError("This experiment's sample identity or status is invalid.")
        _positive_integer(sample["replicate"], "Saved biological replicate")
        fraction = _positive_integer(sample["fraction"], "Saved fraction")
        if manifest["lfq_enabled"] and (
            not sample["condition"].strip() or fraction != 1
        ):
            raise ValueError("This experiment's saved LFQ design is invalid.")
        seen_ids.add(sample["id"])
        seen_files.add(filename.casefold())
    if manifest["status"] == "finished" and (
        any(sample["status"] != "finished" for sample in manifest["samples"])
        or manifest["comparison_status"]
        != ("finished" if manifest["lfq_enabled"] else "not_requested")
    ):
        raise ValueError("This experiment's completion status is inconsistent.")
    return manifest


class _SampleLogger:
    """Keep sample logs locally while forwarding the complete experiment log."""

    def __init__(self, parent, local, sample):
        self.parent, self.local, self.sample = parent, local, sample

    def log(self, message, level=0):
        self.parent.log(
            f"[{self.sample['id']}: {self.sample['file']}] {message}", level
        )
        self.local.log(message, level)


def _shareable_database(child, root: Path) -> str | None:
    """Only reuse a native-generated database inside this controlled experiment."""
    generated = getattr(child, "_generated_database", None)
    if not isinstance(generated, (str, Path)):
        return None
    try:
        path = Path(generated).resolve()
        if path.is_file() and path.is_relative_to(root.resolve()):
            return str(path)
    except (OSError, RuntimeError):
        pass
    return None


def run_experiment(workflow, execute_sample) -> bool:
    """Run all selected samples sequentially, then compare only a complete cohort.

    The master parameter object and its persisted snapshot remain untouched.
    Each sample has a separate snapshot, results and logs. Only the global PID
    directory is shared so the existing cancellation controls still apply.
    """
    try:
        samples, lfq = preflight_experiment(workflow.params)
    except (OSError, ValueError, TypeError) as error:
        workflow.logger.log(f"ERROR: Cannot start experiment: {error}")
        return False
    from src.workflow.CommandExecutor import CommandExecutor
    from src.workflow.FileManager import FileManager
    from src.workflow.Logger import Logger
    from src.workflow.ParameterManager import ParameterManager

    root = Path(workflow.results_dir)
    root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": 1,
        "lfq_enabled": lfq,
        "status": "running",
        "comparison_status": "pending" if lfq else "not_requested",
        "samples": samples,
    }
    _atomic_manifest(root, manifest)
    shared_database = None
    for sample in samples:
        sample["status"] = "running"
        _atomic_manifest(root, manifest)
        logger = workflow.logger
        child = None
        try:
            snapshot = sample_snapshot(root, sample["id"])
            results = sample_results(root, sample["id"])
            results.mkdir(parents=True, exist_ok=True)
            child = copy(workflow)
            for private in ("_shared_database", "_generated_database"):
                if hasattr(child, private):
                    delattr(child, private)
            child._experiment_sample_id = sample["id"]
            if shared_database is not None:
                child._shared_database = shared_database
            child.results_dir = results
            child.params = deepcopy(workflow.params)
            child.params["mzML-files"] = sample["path"]
            child.params.pop("_acquisition_fingerprint", None)
            child.params["_sample_design"] = [
                {
                    key: sample["path"] if key == "file" else sample[key]
                    for key in ("file", "condition", "replicate", "fraction")
                }
            ]
            child.parameter_manager = ParameterManager(workflow.workflow_dir)
            child.parameter_manager.params_file = snapshot
            child.parameter_manager.write_parameters(child.params)
            child.file_manager = FileManager(workflow.workflow_dir, results)
            child.logger = _SampleLogger(
                workflow.logger, Logger(snapshot.parent), sample
            )
            logger = child.logger
            child.executor = CommandExecutor(
                workflow.workflow_dir, child.logger, child.parameter_manager
            )
            child.executor.pid_dir.mkdir(parents=True, exist_ok=True)
            logger.log("Starting sample workflow.")
            if not execute_sample(child):
                raise ValueError(
                    "The sample workflow reported failure; inspect this sample's saved log."
                )
            sample["status"] = "finished"
        # Record arbitrary tool/callback failures without discarding other samples.
        except Exception as error:  # noqa: BLE001
            sample["status"], sample["error"] = "failed", str(error)
            logger.log(
                f"ERROR: Sample {sample['id']} ({sample['file']}) failed: {error}"
            )
        finally:
            if child is not None:
                shared_database = _shareable_database(child, root) or shared_database
        _atomic_manifest(root, manifest)
    if any(sample["status"] != "finished" for sample in samples):
        manifest["status"] = "failed"
        if lfq:
            manifest["comparison_status"] = "failed"
            manifest["error"] = (
                "Comparison was not generated because at least one sample failed."
            )
        _atomic_manifest(root, manifest)
        return False
    if lfq:
        try:
            from src.naseweis_comparison import build_comparison

            build_comparison(
                root,
                samples,
                measure=workflow.params.get("lfq-measure", "intensity"),
                normalization=workflow.params.get("lfq-normalization", "total"),
            )
            manifest["comparison_status"] = "finished"
        except (OSError, ValueError) as error:
            manifest["status"], manifest["comparison_status"] = "failed", "failed"
            manifest["error"] = str(error)
            workflow.logger.log(f"ERROR: Experiment comparison failed: {error}")
            _atomic_manifest(root, manifest)
            return False
    manifest["status"] = "finished"
    _atomic_manifest(root, manifest)
    return True
