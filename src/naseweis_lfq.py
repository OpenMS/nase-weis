"""Per-acquisition RNA LFQ following Herbert's TOPPAS extraction workflow.

The workflow extracts MS1 features from NASE targets and exports all native
feature annotations. It does not link acquisitions or normalize intensities.
"""

import csv
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from src.naseweis_tools import tool_flags, tool_schema

LFQ_TOOLS = ("FeatureFinderMetaboIdent", "TextExporter")
LFQ_DESCRIPTION = (
    "Extract MS1 features from this acquisition using the identified RNA "
    "oligonucleotides, then export their intensities and native annotations. "
    "This follows the TOPPAS workflow without cross-run linking or normalization."
)
TOPPAS_SOURCE = "Herbert/pipeline2_fixed.toppas"
TOPPAS_SHA256 = "1cecf741a05e22123d3471c841910097329552b3f662064c1b8a5077f2141000"

# Preserve the scientific settings in the runnable Herbert TOPPAS workflow.
# isotope_pmin takes precedence over n_isotopes, as in the native tool.
FEATURE_FINDER_DEFAULTS = {
    "extract:mz_window": 50.0,
    "extract:rt_window": 240.0,
    "extract:n_isotopes": 2,
    "extract:isotope_pmin": 0.05,
    "detect:peak_width": 0.0,
    "detect:min_peak_width": 0.2,
    "detect:signal_to_noise": 0.8,
    "model:type": "none",
    "model:add_zeros": 0.2,
    "model:unweighted_fit": False,
    "model:no_imputation": False,
    "model:each_trace": False,
    "model:check:min_area": 1.0,
    "model:check:boundaries": 0.5,
    "model:check:width": 10.0,
    "model:check:asymmetry": 10.0,
    "EMGScoring:max_iteration": 100,
    "EMGScoring:init_mom": False,
}
# These layout options are workflow-owned: metadata includes the RNA labels.
# Explicit TSV fixes the output format while retaining the TOPPAS export layout.
TEXT_EXPORTER_DEFAULTS = {
    "no_ids": True,
    "feature:add_metavalues": 0,
    "feature:minimal": False,
    "out_type": "tsv",
    "quoting": "none",
    "replacement": "_",
}
FEATURE_FINDER_INPUTS = {"in", "id"}
FEATURE_FINDER_OUTPUTS = {
    "out", "lib_out", "chrom_out", "candidates_out", "trafo_out"
}
TEXT_EXPORTER_INPUTS = {"in"}
TEXT_EXPORTER_OUTPUTS = {
    "out", "consensus:centroids", "consensus:elements", "consensus:features"
}
LFQ_FEATURES = "quantification/features.featureXML"
LFQ_TABLE = "quantification/features.tsv"
LFQ_TARGET_COLUMNS = (
    "CompoundName", "SumFormula", "Mass", "Charge", "RetentionTime",
    "RetentionTimeRange", "IsoDistribution",
)


def _native_value(value, entry: dict):
    """Validate schema types and restrictions before an expensive search starts."""
    native_type = entry["type"]
    if native_type == "bool" or (
        native_type == "string" and entry.get("restrictions") == "true,false"
    ):
        if isinstance(value, bool):
            boolean = value
        elif isinstance(value, str) and value.lower() in ("true", "false"):
            boolean = value.lower() == "true"
        else:
            raise ValueError("must be true or false")
        return boolean if native_type == "bool" else str(boolean).lower()
    if native_type in ("int", "double", "float"):
        if isinstance(value, bool) or not isinstance(value, (int, float, str)):
            raise ValueError("must be a finite number")
        try:
            number = float(value)
        except (ValueError, OverflowError):
            raise ValueError("must be a finite number") from None
        if not math.isfinite(number):
            raise ValueError("must be a finite number")
        if native_type == "int" and not number.is_integer():
            raise ValueError("must be an integer")
        limits = entry.get("restrictions", "")
        if ":" in limits:
            lower, upper = limits.split(":", 1)
            if (lower and number < float(lower)) or (upper and number > float(upper)):
                raise ValueError(f"must be within {limits}")
        return int(number) if native_type == "int" else number
    if not isinstance(value, str):
        raise TypeError("must be text")
    choices = entry.get("restrictions", "")
    if choices and value not in choices.split(",") and not (
        value == "" and entry.get("value") == ""
    ):
        raise ValueError(f"must be one of {choices}")
    return value


def prepare_lfq(params: dict, parameter_manager, spectra: str) -> str | None:
    """Prepare enabled LFQ settings in a copied run snapshot, without writing it.

    Workspace paths and the feature-table layout cannot be overridden by
    imported native settings. Scientific overrides retain normal precedence.
    """
    if not isinstance(spectra, str) or Path(spectra).suffix.lower() not in {
        ".mzml", ".raw", ".d"
    }:
        return "Label-free quantification requires an mzML, RAW or .d acquisition."
    defaults = params.get("_defaults", {})
    if not isinstance(defaults, dict) or not isinstance(params.get("_flag_params", {}), dict):
        return "Invalid stored LFQ parameter defaults or flags; reset or repair the parameters."
    prepared = {}
    for tool, preset, inputs, outputs in (
        (LFQ_TOOLS[0], FEATURE_FINDER_DEFAULTS, FEATURE_FINDER_INPUTS, FEATURE_FINDER_OUTPUTS),
        (LFQ_TOOLS[1], TEXT_EXPORTER_DEFAULTS, TEXT_EXPORTER_INPUTS, TEXT_EXPORTER_OUTPUTS),
    ):
        stored_defaults, overrides = defaults.get(tool, {}), params.get(tool, {})
        if not isinstance(stored_defaults, dict) or not isinstance(overrides, dict):
            return f"Invalid stored {tool} parameters; reset or repair its settings."
        if not parameter_manager.create_ini(tool):
            return f"Cannot read {tool} parameters. Install the pinned OpenMS nightly tools and add them to PATH."
        path = parameter_manager.ini_dir / f"{tool}.ini"
        try:
            schema = tool_schema(path)
        except (OSError, ET.ParseError, KeyError):
            return f"Cannot read {tool} parameter schema; reinstall the pinned OpenMS nightly runtime."
        required = inputs | {"out"}
        if not required <= schema.keys():
            return f"This {tool} installation lacks required LFQ inputs or outputs; use the pinned nightly runtime."
        # The executor merges defaults again, so remove file paths from both
        # stored layers, including newly added native file options in a nightly.
        files = inputs | outputs | {
            key for key, entry in schema.items()
            if entry["type"] in {"input-file", "output-file"}
        }
        clean_defaults = {key: value for key, value in stored_defaults.items() if key not in files}
        settings = preset | clean_defaults | {key: value for key, value in overrides.items() if key not in files}
        if tool == "TextExporter":
            settings.update(TEXT_EXPORTER_DEFAULTS)
            clean_defaults.update(TEXT_EXPORTER_DEFAULTS)
        unknown = sorted(settings.keys() - schema.keys())
        if unknown:
            return f"Unsupported parameters for {tool}: {', '.join(unknown)}. Use the pinned nightly runtime or remove stale settings."
        for key, value in settings.items():
            entry = schema[key]
            try:
                if entry["is_list"]:
                    if not isinstance(value, (list, tuple)):
                        raise ValueError("must be a list")
                    settings[key] = [_native_value(item, entry) for item in value]
                else:
                    settings[key] = _native_value(value, entry)
            except (ValueError, TypeError) as error:
                return f"Invalid {tool} parameter {key}: {error}."
        if tool == "FeatureFinderMetaboIdent" and (
            settings["extract:rt_window"] == 0 and settings["detect:peak_width"] == 0
        ):
            return "FeatureFinderMetaboIdent requires a positive extract:rt_window or detect:peak_width to define the extraction range."
        prepared[tool] = (settings, clean_defaults, tool_flags(path))
    for tool, (settings, clean_defaults, flags) in prepared.items():
        params[tool] = settings
        params.setdefault("_defaults", {})[tool] = clean_defaults
        params.setdefault("_flag_params", {})[tool] = flags
    return None


def run_lfq(workflow, targets_path: str | Path) -> bool:
    """Extract this run's features, then export them; propagate either failure."""
    if not Path(targets_path).is_file():
        workflow.logger.log("ERROR: NASE did not produce the RNA targets required for label-free quantification.")
        return False
    features = workflow.results_dir / LFQ_FEATURES
    table = workflow.results_dir / LFQ_TABLE
    try:
        with Path(targets_path).open(encoding="utf-8-sig", newline="") as source:
            records = csv.reader(source, delimiter="\t", strict=True)
            if tuple(next(records, ())) != LFQ_TARGET_COLUMNS:
                raise ValueError("the native NASE target header is missing or invalid")
            count = 0
            for number, row in enumerate(records, 2):
                if not row or not any(row):
                    continue
                if len(row) != len(LFQ_TARGET_COLUMNS):
                    raise ValueError(f"target row {number} has an invalid column count")
                count += 1
    except (OSError, UnicodeError, csv.Error, ValueError) as error:
        workflow.logger.log(f"ERROR: Cannot read RNA LFQ targets: {error}.")
        return False
    features.parent.mkdir(parents=True, exist_ok=True)
    if not count:
        # FFMetaboIdent rejects an empty target list. An empty native feature
        # map represents this successful zero-identification case accurately.
        import pyopenms as poms

        poms.FeatureXMLFile().store(str(features), poms.FeatureMap())
        workflow.logger.log("No identified RNA targets; label-free quantification has no features to extract.")
    else:
        workflow.logger.log("Extracting RNA MS1 features for label-free quantification...")
        if not workflow.executor.run_topp(
            "FeatureFinderMetaboIdent",
            {"in": [workflow.params["mzML-files"]], "id": [str(targets_path)], "out": [str(features)]},
        ):
            return False
    workflow.logger.log("Exporting RNA feature intensities and annotations...")
    return bool(workflow.executor.run_topp("TextExporter", {"in": [str(features)], "out": [str(table)]}))


_NUMERIC_COLUMNS = {
    "rt", "mz", "intensity", "charge", "width", "quality", "rt_quality",
    "mz_quality", "rt_start", "rt_end", "expected_rt", "rt_deviation",
    "expected_im", "num_of_masstraces", "isotope_probability",
}


def read_lfq_table(path: str | Path) -> pd.DataFrame:
    """Read native TextExporter features, preserving labels and all meta columns.

    A header with no rows is a valid empty export. Truncated or malformed files
    raise ValueError rather than silently reporting corrupt intensities.
    """
    header, rows = None, []
    with Path(path).open(encoding="utf-8-sig", newline="") as source:
        try:
            for number, row in enumerate(csv.reader(source, delimiter="\t", strict=True), 1):
                if not row or not any(row):
                    continue
                marker = row[0].lstrip("#")
                if marker == "rt" or (marker == "FEATURE" and row[0].startswith("#")):
                    if not {"mz", "intensity"} <= set(row):
                        raise ValueError(f"Invalid feature header on line {number}.")
                    candidate = row[1:] if marker == "FEATURE" else row[:]
                    candidate[0] = candidate[0].lstrip("#")
                    if header is not None or len(set(candidate)) != len(candidate):
                        raise ValueError(f"Invalid or duplicate feature header on line {number}.")
                    header = candidate
                    continue
                if row[0].startswith("#"):
                    continue
                if header is None:
                    raise ValueError(f"Missing TextExporter feature header before line {number}.")
                if row[0] in {"RUN", "PROTEIN", "PEPTIDE", "UNASSIGNEDPEPTIDE"}:
                    continue
                if row[0] == "FEATURE":
                    row = row[1:]
                if len(row) != len(header):
                    raise ValueError(f"Invalid feature row on line {number}: expected {len(header)} columns, got {len(row)}.")
                rows.append(row)
        except csv.Error as error:
            raise ValueError(f"Invalid TextExporter feature table: {error}") from error
    if header is None:
        raise ValueError("Missing TextExporter feature header; the export may be truncated.")
    table = pd.DataFrame(rows, columns=header)
    for column in table.columns.intersection(list(_NUMERIC_COLUMNS)):
        values = table[column].replace({"": None, "NA": None, "nan": None})
        try:
            table[column] = pd.to_numeric(values, errors="raise")
        except (ValueError, TypeError) as error:
            raise ValueError(f"Invalid numeric values in feature column {column}.") from error
    return table
