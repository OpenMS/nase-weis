"""Keep every nightly parameter and file IO reachable, with explicit compatibility guards."""

from pathlib import Path

import pytest

from src.naseweis_tools import (
    DECOY_INPUTS,
    DECOY_OUTPUTS,
    SEARCH_EXPORTS,
    SEARCH_INPUTS,
    enabled,
    tool_flags,
    tool_schema,
    validate_search_options,
)

SCHEMAS = Path(__file__).parent / "data/naseweis/nightly"


@pytest.mark.parametrize(
    "tool,inputs,outputs",
    [
        ("NucleicAcidSearchEngine", SEARCH_INPUTS, set(SEARCH_EXPORTS) | {"out"}),
        ("DecoyDatabase", DECOY_INPUTS, DECOY_OUTPUTS),
    ],
)
def test_every_nightly_file_option_has_a_workflow_handler(tool, inputs, outputs):
    schema = tool_schema(SCHEMAS / f"{tool}.ini")
    assert {
        key for key, entry in schema.items() if entry["type"] == "input-file"
    } == inputs
    assert {
        key for key, entry in schema.items() if entry["type"] == "output-file"
    } == outputs


def test_new_flags_are_discovered_but_string_booleans_remain_values():
    flags = tool_flags(SCHEMAS / "NucleicAcidSearchEngine.ini")
    assert {
        "precursor:correct_with_flash",
        "scoring:use_mass_accuracy",
        "preprocessing:filter_window_mower",
        "report:precursor_percentage",
    } <= set(flags)
    assert "scoring:single_peak_credit" not in flags


def test_missing_and_broken_schema_retains_legacy_flags(tmp_path):
    path = tmp_path / "missing.ini"
    assert tool_flags(path, ["legacy"]) == ["legacy"]
    path.write_text("<broken")
    assert tool_flags(path, ["legacy"]) == ["legacy"]


@pytest.mark.parametrize(
    "value,result",
    [(True, True), (False, False), ("true", True), ("false", False), (None, False)],
)
def test_imported_booleans(value, result):
    assert enabled(value) is result


@pytest.mark.parametrize(
    "params,exports,cached,reason",
    [
        (
            {"precursor:correct_with_flash": True, "precursor:use_avg_mass": True},
            [],
            False,
            "monoisotopic",
        ),
        (
            {
                "precursor:correct_with_flash": True,
                "precursor:mass_tolerance_unit": "Da",
            },
            [],
            False,
            "precursor tolerance",
        ),
        ({}, ["out_precursor_corrections"], False, "Enable FLASH"),
        (
            {
                "preprocessing:deisotoping": "flash",
                "fragment:mass_tolerance_unit": "Da",
            },
            [],
            False,
            "fragment tolerance",
        ),
        (
            {"modifications:resolve_ambiguities": True, "fragment:ions": ["b", "y"]},
            [],
            False,
            "a-B",
        ),
        ({}, ["digest_out"], True, "cached digest"),
        ({}, ["unknown"], False, "Unknown output"),
    ],
)
def test_incompatible_settings_have_actionable_errors(params, exports, cached, reason):
    assert (
        reason.lower()
        in validate_search_options(params, exports, cached_digest=cached).lower()
    )


def test_compatible_settings_pass():
    assert (
        validate_search_options(
            {
                "precursor:correct_with_flash": True,
                "preprocessing:deisotoping": "flash",
            },
            list(set(SEARCH_EXPORTS) - {"digest_out"}),
            cached_digest=True,
        )
        is None
    )
