"""Runtime schemas and output contracts for NASEWEIS's OpenMS tools."""

import xml.etree.ElementTree as ET
from pathlib import Path

# Native output option -> UI label, relative output path, MIME type.
SEARCH_EXPORTS = {
    "id_out": ("idXML (TOPPView)", "idxml_results/search.idXML", "application/xml"),
    "db_out": (
        "Identification database (OMS)",
        "exports/search.oms",
        "application/octet-stream",
    ),
    "digest_out": (
        "Reusable database digest (OMS)",
        "exports/digest.oms",
        "application/octet-stream",
    ),
    "lfq_out": (
        "Label-free quantification targets",
        "exports/lfq-targets.tsv",
        "text/tab-separated-values",
    ),
    "bedrmod_out": (
        "RNA modifications (bedRMod)",
        "exports/modifications.bedrmod",
        "text/plain",
    ),
    "theo_ms2_out": (
        "Theoretical MS2 spectra",
        "exports/theoretical-MS2.mzML",
        "application/xml",
    ),
    "exp_ms2_out": (
        "Matched experimental MS2 spectra",
        "exports/experimental-MS2.mzML",
        "application/xml",
    ),
    "out_preprocessed": (
        "Preprocessed MS2 spectra",
        "exports/preprocessed-MS2.mzML",
        "application/xml",
    ),
    "out_precursor_corrections": (
        "FLASH precursor correction audit",
        "exports/precursor-corrections.tsv",
        "text/tab-separated-values",
    ),
}
SEARCH_INPUTS = {"in", "database", "digest", "bedrmod_chebi_mapping"}
DECOY_INPUTS = {"in", "NeighborSearch:in_relevant_proteins"}
DECOY_OUTPUTS = {"out", "NeighborSearch:out_neighbor", "NeighborSearch:out_relevant"}


def tool_schema(path: str | Path) -> dict[str, dict]:
    """Read the executable's ParamXML without a pyOpenMS version dependency."""
    root = ET.parse(path).getroot()
    result = {}

    def walk(node, parts):
        for child in node:
            if child.tag == "NODE":
                walk(child, parts + [child.attrib["name"]])
            elif child.tag in ("ITEM", "ITEMLIST") and "1" in parts:
                key = ":".join(parts[parts.index("1") + 1 :] + [child.attrib["name"]])
                result[key] = dict(child.attrib, is_list=child.tag == "ITEMLIST")

    walk(root, [])
    return result


def tool_flags(path: str | Path, fallback=()) -> list[str]:
    """Discover every actual CLI flag; string-valued true/false options stay values."""
    try:
        return sorted(
            {key for key, entry in tool_schema(path).items() if entry["type"] == "bool"}
            | set(fallback)
        )
    except (OSError, ET.ParseError):
        return list(fallback)


def enabled(value) -> bool:
    """Imported JSON and TOPP parameter files may store booleans as strings."""
    return value.lower() == "true" if isinstance(value, str) else bool(value)


def validate_search_options(
    search: dict, exports: list[str], *, cached_digest=False
) -> str | None:
    """Reject incompatible scientific settings before starting a potentially long run."""
    if enabled(search.get("precursor:correct_with_flash", False)):
        if enabled(search.get("precursor:use_avg_mass", False)):
            return "FLASH precursor correction requires monoisotopic precursor masses (disable use_avg_mass)."
        if search.get("precursor:mass_tolerance_unit", "ppm") != "ppm":
            return "FLASH precursor correction requires precursor tolerance in ppm."
    elif "out_precursor_corrections" in exports:
        return "Enable FLASH precursor correction before requesting its audit export."
    if (
        search.get("preprocessing:deisotoping") == "flash"
        and search.get("fragment:mass_tolerance_unit", "ppm") != "ppm"
    ):
        return "FLASH RNA deisotoping requires fragment tolerance in ppm."
    if enabled(search.get("modifications:resolve_ambiguities", False)):
        ions = search.get("fragment:ions", ["a-B"])
        if "a-B" not in ions:
            return "Modification ambiguity resolution requires a-B fragment ions."
    if cached_digest and "digest_out" in exports:
        return "A cached digest cannot be exported again; choose a FASTA database to generate a new digest."
    unknown = set(exports) - SEARCH_EXPORTS.keys()
    return (
        f"Unknown output selection: {', '.join(sorted(unknown))}" if unknown else None
    )
