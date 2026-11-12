"""Scientific contracts for RNA abundance pooling and relative experiments."""

import json
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from src.naseweis_comparison import (
    build_comparison,
    compare_conditions,
    load_comparison,
)


def samples(count=2):
    return [
        {
            "id": f"sample_{index:03d}",
            "file": f"sample-{index}.mzML",
            "path": f"/inputs/sample-{index}.mzML",
            "condition": "control" if index == 1 else "treated",
            "replicate": 1,
            "fraction": 1,
            "status": "finished",
        }
        for index in range(1, count + 1)
    ]


def write_features(root, identity, rows):
    path = root / "samples" / identity / "results" / "quantification" / "features.tsv"
    path.parent.mkdir(parents=True)
    table = pd.DataFrame(rows)
    if table.empty:
        table = pd.DataFrame(columns=["rt", "mz", "intensity", "charge"])
    else:
        table["rt"], table["mz"] = 60.0, 500.25
        if "charge" not in table:
            table["charge"] = -2
    table = table[
        ["rt", "mz", "intensity"]
        + [column for column in table if column not in {"rt", "mz", "intensity"}]
    ]
    table.to_csv(path, sep="\t", index=False)
    path.write_text("#" + path.read_text())
    return path


@pytest.fixture
def experiment(tmp_path):
    design = samples()
    write_features(
        tmp_path,
        "sample_001",
        [
            {
                "label": "pA[Am]CUc",
                "intensity": 10,
                "charge": -2,
                "sum_formula": "C10H20O5",
            },
            {
                "label": "pA[Am]CUc",
                "intensity": 20,
                "charge": -3,
                "sum_formula": "C10H20O5",
            },
            {
                "label": "pA[Am]CUc+[Na-H]",
                "intensity": 30,
                "charge": -2,
                "sum_formula": "C10H19NaO5",
            },
            {"label": "ACGU", "intensity": 40, "sum_formula": "C20H40O10"},
            {"label": "AAAA", "intensity": 0, "sum_formula": "C5H10O2"},
        ],
    )
    write_features(
        tmp_path,
        "sample_002",
        [
            {"label": "pA[Am]CUc", "intensity": 120, "sum_formula": "C10H20O5"},
            {"label": "UGCA", "intensity": 120, "sum_formula": "C20H40O10"},
        ],
    )
    return tmp_path, design


def test_charge_adduct_pooling_preserves_modified_rna_chemistry_and_missing_values(
    experiment,
):
    root, design = experiment
    original = deepcopy(design)
    result = build_comparison(root, design)
    raw, relative = result["raw"], result["abundances"]

    assert list(raw.index) == ["ACGU", "UGCA", "pA[Am]CUc"]
    assert list(raw.columns) == ["sample_001", "sample_002"]
    assert raw.loc["pA[Am]CUc"].tolist() == [60, 120]
    assert raw.loc["ACGU", "sample_001"] == 40
    assert pd.isna(raw.loc["ACGU", "sample_002"])
    assert pd.isna(raw.loc["UGCA", "sample_001"])
    assert relative.loc["pA[Am]CUc"].tolist() == [0.6, 0.5]
    assert relative.sum().tolist() == [1, 1]
    assert result["metadata"]["sample_totals"] == {"sample_001": 100, "sample_002": 240}
    assert result["metadata"]["factors"] == {"sample_001": 100, "sample_002": 240}
    assert design == original
    assert {path.name for path in (root / "comparison").iterdir()} == {
        "raw-abundances.tsv",
        "abundances.tsv",
        "sample-design.tsv",
        "comparison.json",
    }
    restored = load_comparison(root)
    for key in ["raw", "abundances", "design"]:
        assert_frame_equal(result[key], restored[key])
    assert restored["metadata"] == result["metadata"]


def test_none_normalization_retains_raw_native_area(experiment):
    root, design = experiment
    result = build_comparison(root, design, normalization="none")
    assert_frame_equal(result["raw"], result["abundances"])
    assert result["metadata"]["factors"] == {"sample_001": 1.0, "sample_002": 1.0}
    summary = compare_conditions(
        result["abundances"], result["design"], "treated", "control"
    )
    assert summary.loc["pA[Am]CUc", "ratio"] == 2
    assert summary.loc["pA[Am]CUc", "log2_fold_change"] == 1
    assert pd.isna(summary.loc["ACGU", "ratio"])


def test_fitted_area_is_explicit_and_not_substituted_with_native_intensity(tmp_path):
    write_features(
        tmp_path,
        "sample_001",
        [
            {
                "label": "ACGU",
                "intensity": 999,
                "model_area": "42.25",
                "sum_formula": "C10H20",
            },
        ],
    )
    result = build_comparison(
        tmp_path, samples(1), measure="model_area", normalization="none"
    )
    assert result["raw"].loc["ACGU", "sample_001"] == 42.25
    assert result["metadata"]["measure"] == "model_area"
    with pytest.raises(ValueError, match="model_area"):
        root = tmp_path / "missing"
        write_features(root, "sample_001", [{"label": "ACGU", "intensity": 999}])
        build_comparison(root, samples(1), measure="model_area")


def test_only_trailing_adduct_suffix_is_removed(tmp_path):
    write_features(
        tmp_path,
        "sample_001",
        [
            {"label": "pA[m5C]Gp+[Na-H]", "intensity": 2},
            {"label": "pA[m5C]Gp", "intensity": 3},
            {"label": "pA[m5C]Gc", "intensity": 4},
            {"label": "A+[native]G", "intensity": 5},
        ],
    )
    raw = build_comparison(tmp_path, samples(1), normalization="none")["raw"]
    assert set(raw.index) == {"pA[m5C]Gp", "pA[m5C]Gc", "A+[native]G"}
    assert raw.loc["pA[m5C]Gp", "sample_001"] == 5


def test_formula_conflicts_for_exact_native_identity_are_rejected(experiment):
    root, design = experiment
    path = root / "samples/sample_002/results/quantification/features.tsv"
    path.write_text(path.read_text().replace("C10H20O5", "C10H99O5"))
    with pytest.raises(ValueError, match="Conflicting sum formulas"):
        build_comparison(root, design)
    assert not (root / "comparison").exists()


def test_unknown_formulas_are_flagged_without_merging_different_isobar_labels(tmp_path):
    write_features(
        tmp_path,
        "sample_001",
        [
            {"label": "ACGU", "intensity": 1},
            {"label": "UGCA", "intensity": 2},
        ],
    )
    result = build_comparison(tmp_path, samples(1))
    assert set(result["raw"].index) == {"ACGU", "UGCA"}
    assert any(
        "incomplete formula" in warning for warning in result["metadata"]["warnings"]
    )


@pytest.mark.parametrize(
    "measure,value",
    [
        ("intensity", -1),
        ("intensity", float("inf")),
        ("model_area", -1),
        ("model_area", float("nan")),
        ("model_area", float("inf")),
        ("model_area", "invalid"),
    ],
)
def test_invalid_native_areas_fail_without_creating_comparison(
    tmp_path, measure, value
):
    row = {"label": "ACGU", "intensity": 1, measure: value}
    write_features(tmp_path, "sample_001", [row])
    with pytest.raises(ValueError, match="numeric|finite|negative"):
        build_comparison(tmp_path, samples(1), measure=measure)
    assert not (tmp_path / "comparison").exists()


@pytest.mark.parametrize("label", ["", " ", "+[Na-H]"])
def test_missing_or_empty_cleaned_rna_labels_are_rejected(tmp_path, label):
    write_features(tmp_path, "sample_001", [{"label": label, "intensity": 1}])
    with pytest.raises(ValueError, match="label"):
        build_comparison(tmp_path, samples(1))


def test_valid_empty_or_zero_area_samples_retain_missing_matrix_columns(tmp_path):
    design = samples(3)
    design[2]["status"] = "failed"
    write_features(tmp_path, "sample_001", [])
    write_features(tmp_path, "sample_002", [{"label": "ACGU", "intensity": 0}])
    result = build_comparison(tmp_path, design)
    assert result["raw"].empty and result["abundances"].empty
    assert list(result["raw"].columns) == [sample["id"] for sample in design]
    assert result["metadata"]["factors"] == dict.fromkeys(result["raw"].columns)
    assert any(
        "did not finish" in warning for warning in result["metadata"]["warnings"]
    )
    assert_frame_equal(load_comparison(tmp_path)["raw"], result["raw"])


def test_empty_fitted_area_sample_is_valid_without_an_unavailable_metadata_column(
    tmp_path,
):
    write_features(tmp_path, "sample_001", [])
    result = build_comparison(tmp_path, samples(1), measure="model_area")
    assert result["abundances"].empty


def test_a_failed_or_empty_sample_does_not_receive_invented_zero_measurements(tmp_path):
    design = samples(2)
    design[1]["status"] = "failed"
    write_features(tmp_path, "sample_001", [{"label": "ACGU", "intensity": 5}])
    result = build_comparison(tmp_path, design)
    assert result["raw"].loc["ACGU", "sample_001"] == 5
    assert pd.isna(result["raw"].loc["ACGU", "sample_002"])
    assert pd.isna(result["abundances"].loc["ACGU", "sample_002"])


@pytest.mark.parametrize("labels", [["ACGU", "ACGU"], ["ACGU", "UGCA"]])
def test_area_pooling_or_total_overflow_is_rejected_before_publishing_a_matrix(
    tmp_path, labels
):
    write_features(
        tmp_path,
        "sample_001",
        [{"label": label, "intensity": 1e308} for label in labels],
    )
    with pytest.raises(ValueError, match="overflowed"):
        build_comparison(tmp_path, samples(1))
    assert not (tmp_path / "comparison").exists()


def test_sample_output_cannot_escape_the_experiment_through_a_symlink(tmp_path):
    root = tmp_path / "experiment"
    outside = tmp_path / "outside.tsv"
    outside.write_text("#rt\tmz\tintensity\n")
    link = root / "samples/sample_001/results/quantification/features.tsv"
    link.parent.mkdir(parents=True)
    link.symlink_to(outside)
    with pytest.raises(ValueError, match="outside this experiment"):
        build_comparison(root, samples(1))


def test_biological_replicates_are_weighted_equally_after_technical_means():
    matrix = pd.DataFrame(
        {
            "sample_001": [20.0, np.nan, 0],
            "sample_002": [40.0, 40, 0],
            "sample_003": [90.0, np.nan, 0],
            "sample_004": [10.0, np.nan, 0],
            "sample_005": [30.0, 20, 0],
        },
        index=pd.Index(["ACGU", "[m5C]GGp", "AAAA"], name="label"),
    )
    design = pd.DataFrame(
        {
            "condition": ["treated", "treated", "treated", "control", "control"],
            "replicate": [1, 1, 2, 1, 2],
        },
        index=matrix.columns,
    )
    before_matrix, before_design = matrix.copy(deep=True), design.copy(deep=True)
    result = compare_conditions(matrix, design, "treated", "control")
    assert result.loc["ACGU", "numerator_mean"] == 60
    assert result.loc["ACGU", "numerator_std"] == pytest.approx(
        np.std([30, 90], ddof=1)
    )
    assert result.loc["ACGU", "denominator_mean"] == 20
    assert result.loc["ACGU", "denominator_std"] == pytest.approx(
        np.std([10, 30], ddof=1)
    )
    assert result.loc["ACGU", "numerator_n"] == result.loc["ACGU", "denominator_n"] == 2
    assert result.loc["ACGU", "ratio"] == 3
    assert result.loc["ACGU", "log2_fold_change"] == pytest.approx(np.log2(3))
    assert result.loc["[m5C]GGp", "numerator_n"] == 1
    assert result.loc["[m5C]GGp", "ratio"] == 2
    assert pd.isna(result.loc["[m5C]GGp", "numerator_std"])
    assert pd.isna(result.loc["AAAA", "ratio"])
    assert "p_value" not in result
    assert_frame_equal(matrix, before_matrix)
    assert_frame_equal(design, before_design)


@pytest.mark.parametrize(
    "changes",
    [
        {"id": "../escape"},
        {"replicate": 0},
        {"replicate": 1.5},
        {"replicate": True},
        {"replicate": "bad"},
        {"fraction": 0},
        {"condition": None},
        {"path": ""},
    ],
)
def test_invalid_sample_design_is_rejected_before_reading_files(tmp_path, changes):
    design = samples(1)
    design[0].update(changes)
    with pytest.raises(ValueError):
        build_comparison(tmp_path, design)


@pytest.mark.parametrize("design", [[], ["bad"], samples(1) * 2])
def test_empty_invalid_or_duplicate_sample_id_designs_are_rejected(tmp_path, design):
    with pytest.raises(ValueError):
        build_comparison(tmp_path, design)


@pytest.mark.parametrize(
    "setting",
    [
        {"measure": "unknown"},
        {"normalization": "unknown"},
        {"measure": {}},
        {"normalization": []},
    ],
)
def test_unknown_analysis_methods_are_rejected(tmp_path, setting):
    with pytest.raises(ValueError, match="measure"):
        build_comparison(tmp_path, samples(1), **setting)


@pytest.mark.parametrize(
    "case",
    [
        "same",
        "unknown",
        "missing_design",
        "duplicate_design",
        "bad_replicate",
        "negative",
        "infinite",
        "not_table",
    ],
)
def test_invalid_condition_comparisons_have_actionable_errors(case):
    matrix = pd.DataFrame({"a": [1.0], "b": [2.0]}, index=["ACGU"])
    design = pd.DataFrame(
        {"condition": ["A", "B"], "replicate": [1, 1]}, index=["a", "b"]
    )
    numerator, denominator = "B", "A"
    if case == "same":
        numerator = denominator
    elif case == "unknown":
        numerator = "C"
    elif case == "missing_design":
        design = design.iloc[:1]
    elif case == "duplicate_design":
        design.index = ["a", "a"]
    elif case == "bad_replicate":
        design["replicate"] = 0
    elif case == "negative":
        matrix.loc["ACGU", "a"] = -1
    elif case == "infinite":
        matrix.loc["ACGU", "a"] = np.inf
    else:
        matrix = None
    with pytest.raises(ValueError):
        compare_conditions(matrix, design, numerator, denominator)


@pytest.mark.parametrize(
    "change",
    ["json", "schema", "design", "columns", "rows", "area", "empty", "identities"],
)
def test_corrupt_saved_comparisons_are_not_presented_as_valid(experiment, change):
    root, design = experiment
    build_comparison(root, design)
    output = root / "comparison"
    if change == "json":
        (output / "comparison.json").write_text("{")
    elif change == "schema":
        (output / "comparison.json").write_text(json.dumps({"schema_version": 99}))
    elif change == "design":
        (output / "sample-design.tsv").write_text("wrong\tcolumns\n1\t2\n")
    elif change == "columns":
        path = output / "abundances.tsv"
        path.write_text(path.read_text().replace("sample_001", "sample_999"))
    elif change == "rows":
        path = output / "abundances.tsv"
        path.write_text(path.read_text().replace("ACGU", "AAAA"))
    elif change == "empty":
        (output / "raw-abundances.tsv").write_text("")
    elif change == "identities":
        path = output / "raw-abundances.tsv"
        path.write_text(path.read_text().replace("label\t", "wrong_identity\t"))
    else:
        path = output / "raw-abundances.tsv"
        path.write_text(path.read_text().replace("40\t", "-1\t"))
    with pytest.raises(ValueError):
        load_comparison(root)
