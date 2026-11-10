"""Exact native spectrum joins and faithful RNA fragment labels for viewers."""

from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from src.naseweis_annotations import (
    ANNOTATION_COLUMNS,
    annotate_peaks,
    match_annotations,
    native_spectrum_id,
    read_annotations,
)

NATIVE_ID = "controllerType=0 controllerNumber=1 scan=2593"


def store_annotations(tmp_path, records):
    """Write a real minimal NASE-style idXML with the installed pyOpenMS API."""
    import pyopenms as oms

    protein = oms.ProteinIdentification()
    protein.setIdentifier("native-rna-search")
    protein.setSearchEngine("NucleicAcidSearchEngine")
    identifications = oms.PeptideIdentificationList()
    for record in records:
        identification = oms.PeptideIdentification()
        identification.setIdentifier("native-rna-search")
        identification.setSpectrumReference(record.get("native_id", NATIVE_ID))
        identification.setRT(1750.0)
        identification.setMZ(815.6)
        identification.setScoreType("hyperscore")
        identification.setHigherScoreBetter(True)
        hits = []
        for definition in record.get("hits", [{}]):
            hit = oms.PeptideHit()
            hit.setCharge(definition.get("charge", -2))
            hit.setScore(90.5)
            if definition.get("sequence", "UAAUGp") is not None:
                hit.setMetaValue("label", definition.get("sequence", "UAAUGp"))
            hit.setMetaValue("molecule_type", definition.get("molecule_type", "RNA"))
            annotations = []
            for values in definition.get(
                "annotations", [(305.01806640625, 0.5, "c1", -1)]
            ):
                annotation = oms.PeptideHit_PeakAnnotation()
                (
                    annotation.mz,
                    annotation.intensity,
                    annotation.annotation,
                    annotation.charge,
                ) = values
                annotations.append(annotation)
            hit.setPeakAnnotations(annotations)
            hits.append(hit)
        identification.setHits(hits)
        identifications.push_back(identification)
    destination = tmp_path / "rna.idXML"
    oms.IdXMLFile().store(str(destination), [protein], identifications)
    return destination


@pytest.mark.parametrize(
    "reference, expected",
    [
        (f"ms_run[1]:{NATIVE_ID}", NATIVE_ID),
        ("ms_run[12]:index=153", "index=153"),
        (
            "ms_run[2]:sample=5 period=1 cycle=77 experiment=2",
            "sample=5 period=1 cycle=77 experiment=2",
        ),
        (NATIVE_ID, NATIVE_ID),
        ("index=153", "index=153"),
        ("scan=12", "scan=12"),
        ("ms_run[x]:scan=12", "ms_run[x]:scan=12"),
        ("prefix ms_run[1]:scan=12", "prefix ms_run[1]:scan=12"),
        (" ms_run[1]:scan=12", " ms_run[1]:scan=12"),
        ("ms_run[1]:ms_run[2]:scan=12", "ms_run[2]:scan=12"),
        ("", ""),
    ],
)
def test_native_reference_preserves_the_complete_native_identifier(reference, expected):
    assert native_spectrum_id(reference) == expected


@pytest.mark.parametrize("reference", [None, 123, [], float("nan")])
def test_nontext_native_references_are_actionable_errors(reference):
    with pytest.raises(ValueError, match="spectrum reference must be text"):
        native_spectrum_id(reference)


def test_real_idxml_reads_rna_label_signed_charges_and_all_peak_assignments(tmp_path):
    definitions = [
        {
            "hits": [
                {
                    "annotations": [
                        (305.01806640625, 0.5, "c1", -1),
                        (305.01806640625, 0.5, "w1", -2),
                        (436.1012, 0.25, "a2-B", -1),
                    ]
                },
                {"sequence": "A[m1A]UGp", "charge": -3},
            ]
        },
        {"native_id": "index=7", "hits": [{"sequence": "ACGU", "charge": 2}]},
    ]
    path = store_annotations(tmp_path, definitions)
    table = read_annotations(path)
    assert list(table.columns) == ANNOTATION_COLUMNS
    assert table["native_id"].tolist() == [NATIVE_ID] * 4 + ["index=7"]
    assert table["sequence"].tolist() == ["UAAUGp"] * 3 + ["A[m1A]UGp", "ACGU"]
    assert table["charge"].tolist() == [-2, -2, -2, -3, 2]
    assert sorted(zip(table["ion"], table["fragment_charge"])) == [
        ("a2-B", -1),
        ("c1", -1),
        ("c1", -1),
        ("c1", -1),
        ("w1", -2),
    ]
    assert table["identification_index"].tolist() == [0, 0, 0, 0, 1]
    assert table["hit_index"].tolist() == [0, 0, 0, 1, 0]
    assert table.loc[0, "mz"] == pytest.approx(305.01806640625, abs=1e-10)


@pytest.mark.parametrize("sequence", [None, "", " "])
def test_missing_rna_labels_never_fall_back_to_peptide_sequences(tmp_path, sequence):
    path = store_annotations(tmp_path, [{"hits": [{"sequence": sequence}]}])
    with pytest.raises(ValueError, match="RNA label"):
        read_annotations(path)


@pytest.mark.parametrize(
    "record, message",
    [
        ({"native_id": ""}, "Missing native spectrum ID"),
        ({"hits": [{"charge": 0}]}, "nonzero RNA precursor charge"),
        ({"hits": [{"molecule_type": "protein"}]}, "not RNA"),
        ({"hits": [{"annotations": [(100.0, 1.0, "", -1)]}]}, "fragment ion label"),
    ],
)
def test_invalid_native_annotation_identity_is_rejected(tmp_path, record, message):
    with pytest.raises(ValueError, match=message):
        read_annotations(store_annotations(tmp_path, [record]))


@pytest.mark.parametrize(
    "records", [[], [{"hits": []}], [{"hits": [{"annotations": []}]}]]
)
def test_valid_idxml_without_annotations_returns_a_typed_empty_table(tmp_path, records):
    table = read_annotations(store_annotations(tmp_path, records))
    assert table.empty
    assert list(table.columns) == ANNOTATION_COLUMNS
    assert table["charge"].dtype == np.dtype("int64")
    assert table["mz"].dtype == np.dtype("float64")


def test_missing_unreadable_or_invalid_idxml_reports_the_file_problem(tmp_path):
    with pytest.raises(
        FileNotFoundError, match="NASE idXML annotations are unavailable"
    ):
        read_annotations(tmp_path / "missing.idXML")
    with pytest.raises(OSError, match="Select a NASE idXML file"):
        read_annotations(tmp_path)
    invalid = tmp_path / "invalid.idXML"
    invalid.write_text("not XML")
    with pytest.raises(ValueError, match="use an intact NASE idXML output"):
        read_annotations(invalid)


def test_exact_hit_selection_uses_native_id_rna_label_and_signed_charge(tmp_path):
    table = read_annotations(
        store_annotations(
            tmp_path,
            [
                {"hits": [{}, {"charge": 2}, {"sequence": "ACGU"}]},
                {"native_id": "controllerType=0 controllerNumber=1 scan=2594"},
            ],
        )
    )
    original = table.copy(deep=True)
    selected = match_annotations(table, f"ms_run[1]:{NATIVE_ID}", "UAAUGp", -2)
    assert len(selected) == 1
    assert selected.iloc[0]["charge"] == -2
    assert selected.iloc[0]["hit_index"] == 0
    assert match_annotations(table, "scan=2593", "UAAUGp", -2).empty
    assert match_annotations(table, NATIVE_ID, "UAAUGp", -3).empty
    assert match_annotations(table, NATIVE_ID, "UAAUG", -2).empty
    pd.testing.assert_frame_equal(table, original)


def test_conflicting_duplicate_native_hits_are_not_guessed(tmp_path):
    table = read_annotations(
        store_annotations(
            tmp_path,
            [{"hits": [{}, {"annotations": [(400.0, 0.4, "d3", -1)]}]}],
        )
    )
    with pytest.raises(ValueError, match="Conflicting native annotation hits"):
        match_annotations(table, NATIVE_ID, "UAAUGp", -2)


def test_identical_repeated_hits_and_multiple_assignments_remain_available(tmp_path):
    hit = {"annotations": [(100.0, 0.5, "c1", -1), (100.0, 0.5, "w1", -2)]}
    table = read_annotations(
        store_annotations(tmp_path, [{"hits": [hit, deepcopy(hit)]}])
    )
    selected = match_annotations(table, NATIVE_ID, "UAAUGp", -2)
    assert len(selected) == 4
    assert sorted(selected["ion"]) == ["c1", "c1", "w1", "w1"]


def test_peak_labels_keep_raw_heights_indices_order_and_all_native_input_rows():
    peaks = pd.DataFrame(
        {
            "mz": [400.0, 100.0, 200.0],
            "intensity": [1e7, 2500.0, 4500.0],
            "extra": [3, 1, 2],
        },
        index=["last", "first", "middle"],
    )
    assignments = pd.DataFrame(
        {
            "mz": [100.0, 100.000001, 200.0, 900.0],
            "intensity": [0.5, 0.5, 1.0, 0.1],
            "ion": ["c1", "w1", "a2-B", "d4"],
        }
    )
    original_peaks, original_assignments = (
        peaks.copy(deep=True),
        assignments.copy(deep=True),
    )
    annotated, unmatched = annotate_peaks(peaks, assignments)
    pd.testing.assert_frame_equal(annotated.drop(columns="ion"), original_peaks)
    pd.testing.assert_frame_equal(assignments, original_assignments)
    assert annotated["ion"].tolist() == ["", "c1; w1", "a2-B"]
    assert unmatched["ion"].tolist() == ["d4"]
    assert unmatched["match_status"].tolist() == ["no matching peak"]
    assert annotated["intensity"].tolist() == [1e7, 2500.0, 4500.0]


def test_nearby_or_duplicate_raw_peaks_remain_ambiguous_instead_of_nearest_matching():
    peaks = pd.DataFrame(
        {"mz": [100.0, 100.000015, 200.0, 200.0], "intensity": [1, 9, 2, 3]}
    )
    assignments = pd.DataFrame({"mz": [100.000008, 200.0], "ion": ["d1", "w3"]})
    annotated, unmatched = annotate_peaks(peaks, assignments)
    assert annotated["ion"].tolist() == ["", "", "", ""]
    assert unmatched["ion"].tolist() == ["d1", "w3"]
    assert unmatched["match_status"].tolist() == ["ambiguous peak match"] * 2


def test_serialization_tolerance_and_exact_only_mode_are_explicit():
    peaks = pd.DataFrame({"mz": [100.0], "intensity": [1000.0]})
    assignments = pd.DataFrame({"mz": [100.000009, 100.0], "ion": ["c2", "w2"]})
    annotated, unmatched = annotate_peaks(peaks, assignments)
    assert annotated.loc[0, "ion"] == "c2; w2"
    assert unmatched.empty
    exact, unmatched = annotate_peaks(peaks, assignments, tolerance_da=0)
    assert exact.loc[0, "ion"] == "w2"
    assert unmatched["ion"].tolist() == ["c2"]


def test_repeated_dataframe_indices_do_not_label_the_wrong_peak():
    peaks = pd.DataFrame(
        {"mz": [100.0, 200.0], "intensity": [10, 20]}, index=["peak", "peak"]
    )
    annotated, unmatched = annotate_peaks(
        peaks, pd.DataFrame({"mz": [200.0], "ion": ["d1"]})
    )
    assert annotated["ion"].tolist() == ["", "d1"]
    assert list(annotated.index) == ["peak", "peak"]
    assert unmatched.empty


def test_empty_spectra_or_assignments_still_have_valid_viewer_tables():
    empty_peaks = pd.DataFrame(columns=["mz", "intensity"])
    assignments = pd.DataFrame({"mz": [100.0], "ion": ["w1"]})
    annotated, unmatched = annotate_peaks(empty_peaks, assignments)
    assert annotated.empty and "ion" in annotated
    assert unmatched["ion"].tolist() == ["w1"]
    peaks = pd.DataFrame({"mz": [100.0], "intensity": [1.0]})
    annotated, unmatched = annotate_peaks(peaks, assignments.iloc[0:0])
    assert annotated["ion"].tolist() == [""]
    assert unmatched.empty


@pytest.mark.parametrize("tolerance", [-1, float("inf"), float("nan"), True, "invalid"])
def test_invalid_matching_tolerances_never_guess_peak_labels(tolerance):
    with pytest.raises(ValueError, match="Annotation tolerance"):
        annotate_peaks(
            pd.DataFrame({"mz": [100.0], "intensity": [1.0]}),
            pd.DataFrame({"mz": [100.0], "ion": ["w1"]}),
            tolerance_da=tolerance,
        )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), 0.0, -100.0])
def test_invalid_peak_coordinates_are_reported(value):
    with pytest.raises(ValueError, match="Spectrum m/z"):
        annotate_peaks(
            pd.DataFrame({"mz": [value], "intensity": [1.0]}),
            pd.DataFrame({"mz": [100.0], "ion": ["w1"]}),
        )


def test_missing_columns_and_nonfinite_raw_intensities_are_actionable():
    with pytest.raises(ValueError, match="missing columns intensity"):
        annotate_peaks(
            pd.DataFrame({"mz": [100.0]}), pd.DataFrame({"mz": [100.0], "ion": ["w1"]})
        )
    with pytest.raises(ValueError, match="missing columns ion"):
        annotate_peaks(
            pd.DataFrame({"mz": [100.0], "intensity": [1.0]}),
            pd.DataFrame({"mz": [100.0]}),
        )
    with pytest.raises(ValueError, match="Spectrum intensity"):
        annotate_peaks(
            pd.DataFrame({"mz": [100.0], "intensity": [float("nan")]}),
            pd.DataFrame({"mz": [100.0], "ion": ["w1"]}),
        )
