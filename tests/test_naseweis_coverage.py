"""Native-coordinate, count, modification and vector-export coverage checks."""

import copy
import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from src.naseweis_coverage import (
    build_coverage,
    render_coverage_html,
    render_coverage_svg,
)


def tables(sequence="ACGUACGU", hits=None, matches=None, accession="rna"):
    if hits is None:
        hits = [("ACGp", accession, 1, 3)]
    if matches is None:
        matches = [("ACGp", "scan=1"), ("ACGp", "scan=2")]
    return {
        "NUC": pd.DataFrame({"accession": [accession], "opt_sequence": [sequence]}),
        "OLI": pd.DataFrame(hits, columns=["sequence", "accession", "start", "end"]),
        "OSM": pd.DataFrame(matches, columns=["sequence", "spectra_ref"]),
    }


def test_pool_acquisitions_without_double_counting_native_mappings_or_scans():
    first = tables(
        hits=[("ACGp", "rna", 1, 3)] * 2 + [("ACGp", "rna", 5, 7)],
        matches=[("ACGp", "scan=1")] * 3 + [("ACGp", "scan=2")],
    )
    second = tables(matches=[("ACGp", "scan=1")])
    model = build_coverage([("control", first), ("control", second)], "rna")
    assert model["reference"] == list("ACGUACGU")
    assert model["covered_bases"] == 6
    assert model["coverage"] == 0.75
    group = model["groups"][0]
    assert group["acquisitions"] == 2
    assert group["assignments"] == 3  # Same scan ID in another acquisition is distinct.
    assert group["oligonucleotides"] == 1
    assert [(hit["start"], hit["count"]) for hit in group["hits"]] == [(1, 3), (5, 2)]


def test_ambiguous_scan_supports_exact_sequences_without_inventing_identity():
    data = tables(
        hits=[("ACGp", "rna", 1, 3), ("A[m5C]Gp", "rna", 1, 3)],
        matches=[("ACGp", "scan=1"), ("A[m5C]Gp", "scan=1")],
    )
    model = build_coverage([("sample", data)], "rna")
    assert [hit["sequence"] for hit in model["groups"][0]["hits"]] == ["ACGp", "A[m5C]Gp"]
    assert [hit["count"] for hit in model["groups"][0]["hits"]] == [1, 1]
    assert model["groups"][0]["assignments"] == 2
    assert "not independent totals" in render_coverage_svg(model)


def test_modifications_are_native_bracket_tokens_at_reference_positions():
    data = tables(
        hits=[("pA[m5C]Gp", "rna", 5, 7), ("[5'-p]ACGU[3'-c]", "rna", 1, 4)],
        matches=[("pA[m5C]Gp", "scan=1"), ("[5'-p]ACGU[3'-c]", "scan=2")],
    )
    model = build_coverage([("sample", data)], "rna")
    first, modified = model["groups"][0]["hits"]
    assert first["modifications"] == []
    assert modified["modifications"] == [{"position": 6, "token": "[m5C]"}]
    svg = render_coverage_svg(model)
    assert svg.count('class="modification"') == 1
    assert "Position 6: [m5C]" in svg
    assert "terminal p/c groups are not modification calls" in svg


def test_different_modified_labels_keep_their_own_spectral_counts():
    data = tables(
        hits=[("A[m5C]Gp", "rna", 1, 3), ("ACGp", "rna", 1, 3)],
        matches=[("ACGp", "scan=1"), ("ACGp", "scan=2"), ("A[m5C]Gp", "scan=3")],
    )
    model = build_coverage([("sample", data)], "rna")
    assert {hit["sequence"]: hit["count"] for hit in model["groups"][0]["hits"]} == {
        "ACGp": 2, "A[m5C]Gp": 1,
    }


@pytest.mark.parametrize("sequence,token", [
    ("AYGp", "Y"), ("AIGp", "I"), ("ADGp", "D"),
    ("*ACGp", None), ("[5'-p*]ACGp", None),
])
def test_native_single_letter_modifications_and_terminal_phosphorothioates(sequence, token):
    data = tables(sequence="ACG", hits=[(sequence, "rna", 1, 3)], matches=[(sequence, "scan=1")])
    model = build_coverage([("sample", data)], "rna")
    hit = model["groups"][0]["hits"][0]
    assert model["covered_bases"] == 3
    assert hit["modifications"] == ([] if token is None else [{"position": 2, "token": token}])
    svg = render_coverage_svg(model)
    assert svg.count('class="modification"') == (0 if token is None else 1)


def test_unique_scan_count_falls_back_to_rows_for_missing_references():
    data = tables(matches=[("ACGp", "scan=1"), ("ACGp", None), ("ACGp", None)])
    assert build_coverage([("sample", data)], "rna")["groups"][0]["assignments"] == 3
    data["OSM"] = data["OSM"].drop(columns="spectra_ref")
    assert build_coverage([("sample", data)], "rna")["groups"][0]["assignments"] == 3


def test_same_sequence_on_multiple_accessions_does_not_inflate_selected_accession():
    data = tables(hits=[("ACGp", "rna", 1, 3), ("ACGp", "other", 1, 3)])
    model = build_coverage([("sample", data)], "rna")
    assert len(model["groups"][0]["hits"]) == 1
    assert model["groups"][0]["assignments"] == 2


def test_missing_osm_is_unknown_and_does_not_turn_into_zero_count():
    data = tables()
    del data["OSM"]
    model = build_coverage([("sample", data), ("sample", tables())], "rna")
    assert model["groups"][0]["hits"][0]["count"] is None
    assert model["groups"][0]["assignments"] is None
    assert len(model["warnings"]) == 1
    assert 'data-count="unknown"' in render_coverage_svg(model)


def test_empty_matches_have_zero_count_and_still_show_exported_mapping():
    model = build_coverage([("sample", tables(matches=[]))], "rna")
    assert model["groups"][0]["hits"][0]["count"] == 0
    assert model["groups"][0]["assignments"] == 0
    assert model["covered_bases"] == 3


def test_empty_group_stays_in_order_and_does_not_require_this_accession():
    other = tables(accession="other", hits=[("ACGp", "other", 1, 3)])
    model = build_coverage([("control", tables()), ("treatment", other)], "rna")
    assert [group["name"] for group in model["groups"]] == ["control", "treatment"]
    assert model["groups"][1]["coverage"] == 0
    assert model["groups"][1]["hits"] == []
    assert "No mapped oligonucleotides in this block" in render_coverage_svg(model)


def test_reference_without_hits_renders_zero_coverage():
    data = tables(hits=[], matches=[])
    model = build_coverage([("sample", data)], "rna")
    assert model["coverage"] == 0
    assert "0 / 8 bases mapped (0.0%)" in render_coverage_svg(model)


@pytest.mark.parametrize("changed", ["ACGUACGA", "A[m5C]GUACGU"])
def test_conflicting_reference_sequences_cannot_share_a_plot(changed):
    with pytest.raises(ValueError, match="Conflicting reference sequences"):
        build_coverage([("control", tables()), ("treatment", tables(sequence=changed))], "rna")


def test_multiple_reference_rows_with_different_sequence_are_rejected():
    data = tables()
    data["NUC"] = pd.concat([data["NUC"], pd.DataFrame({"accession": ["rna"], "opt_sequence": ["CCCC"]})])
    with pytest.raises(ValueError, match="Conflicting reference sequences"):
        build_coverage([("sample", data)], "rna")


def test_terminal_reference_aliases_share_an_axis():
    first = tables(sequence="pACGUACGUp")
    second = tables(sequence="[5'-p]ACGUACGU[3'-p]")
    model = build_coverage([("control", first), ("treatment", second)], "rna")
    assert model["length"] == 8


@pytest.mark.parametrize("sequence", [None, "", pd.NA])
def test_no_reference_sequence_is_an_actionable_error(sequence):
    with pytest.raises(ValueError, match="No reference sequence was exported"):
        build_coverage([("sample", tables(sequence=sequence))], "rna")


def test_missing_reference_for_mapped_sample_is_not_borrowed_silently():
    data = tables()
    data["NUC"] = pd.DataFrame(columns=["accession", "opt_sequence"])
    with pytest.raises(ValueError, match="No reference sequence was exported for mapped accession"):
        build_coverage([("control", tables()), ("treatment", data)], "rna")


@pytest.mark.parametrize("field, value", [
    ("start", 0), ("start", -1), ("start", 1.5), ("start", float("nan")),
    ("end", float("inf")), ("end", True), ("end", "bad"), ("end", None),
])
def test_invalid_coordinates_fail_instead_of_being_truncated(field, value):
    data = tables()
    data["OLI"][field] = value
    with pytest.raises(ValueError, match="positive integer"):
        build_coverage([("sample", data)], "rna")


@pytest.mark.parametrize("start,end", [(4, 3), (7, 9)])
def test_out_of_range_coordinates_are_rejected(start, end):
    data = tables(hits=[("ACGp", "rna", start, end)])
    with pytest.raises(ValueError, match="exceed the reference"):
        build_coverage([("sample", data)], "rna")


def test_coordinate_length_must_match_native_modified_sequence_length():
    data = tables(hits=[("A[m5C]Gp", "rna", 1, 4)])
    with pytest.raises(ValueError, match="coordinate span"):
        build_coverage([("sample", data)], "rna")


@pytest.mark.parametrize("section,column,message", [
    ("NUC", "accession", "nucleotide accession"),
    ("OLI", "accession", "oligonucleotide accession"),
    ("OLI", "start", "incomplete oligonucleotide coordinates"),
])
def test_missing_native_columns_fail_clearly(section, column, message):
    data = tables()
    data[section] = data[section].drop(columns=column)
    with pytest.raises(ValueError, match=message):
        build_coverage([("sample", data)], "rna")


@pytest.mark.parametrize("section", ["OLI", "OSM"])
@pytest.mark.parametrize("sequence", [None, ""])
def test_missing_native_oligo_identity_is_rejected(section, sequence):
    data = tables()
    data[section]["sequence"] = sequence
    with pytest.raises(ValueError, match="native RNA sequence"):
        build_coverage([("sample", data)], "rna")


def test_empty_or_invalid_group_selection_is_actionable():
    with pytest.raises(ValueError, match="No acquisition results"):
        build_coverage([], "rna")
    with pytest.raises(ValueError, match="Choose an RNA accession"):
        build_coverage([("sample", tables())], "")
    with pytest.raises(ValueError, match="text label and native mzTab tables"):
        build_coverage([("", tables())], "rna")


def test_shared_color_scale_is_independent_of_other_groups():
    model = build_coverage([("control", tables()), ("treatment", tables(matches=[("ACGp", "scan=1")]))], "rna")
    svg = ET.fromstring(render_coverage_svg(model))
    hits = [element for element in svg.iter() if element.get("class") == "mapped-oligo"]
    assert [element.get("fill") for element in hits] == ["#fdc78d", "#fcfdbf"]
    assert "129+" in render_coverage_svg(model)


def test_overlapping_hits_pack_to_different_lanes_but_nonoverlap_reuses_lane():
    data = tables(hits=[("ACGp", "rna", 1, 3), ("CGUp", "rna", 2, 4), ("ACGp", "rna", 5, 7)])
    model = build_coverage([("sample", data)], "rna")
    svg = ET.fromstring(render_coverage_svg(model))
    hits = [element for element in svg.iter() if element.get("class") == "mapped-oligo"]
    assert hits[0].get("y") != hits[1].get("y")
    assert hits[0].get("y") == hits[2].get("y")


def test_wrap_continues_bar_at_native_positions_and_marks_modification_once():
    data = tables(sequence="A" * 60, hits=[("A[m5C]AAp", "rna", 49, 52)], matches=[("A[m5C]AAp", "scan=1")])
    model = build_coverage([("sample", data)], "rna")
    svg = render_coverage_svg(model)
    assert "Positions 1–50" in svg
    assert "Positions 51–60" in svg
    assert svg.count('class="mapped-oligo"') == 2
    assert svg.count('class="modification"') == 1
    assert "Position 50: [m5C]" in svg
    assert "4 / 60 bases mapped (6.7%)" in svg


@pytest.mark.parametrize("size", [0, 9, 101, 10.0, True])
def test_invalid_block_size_does_not_produce_broken_vector_layout(size):
    with pytest.raises(ValueError, match="between 10 and 100"):
        render_coverage_svg(build_coverage([("sample", tables())], "rna"), block_size=size)


def test_export_is_standalone_parseable_vector_and_escapes_user_text():
    accession = 'RNA<script>alert("x")</script>'
    group = '<img src=x onerror="x">'
    data = tables(accession=accession, hits=[("A[m5C]Gp", accession, 1, 3)], matches=[])
    model = build_coverage([(group, data)], accession)
    svg, html = render_coverage_svg(model), render_coverage_html(model)
    assert ET.fromstring(svg).tag == "{http://www.w3.org/2000/svg}svg"
    assert "<script>" not in html
    assert "<img" not in html
    assert "RNA&lt;script&gt;" in html
    assert "&lt;img" in html
    assert "<svg " in html
    assert "<title>" in svg
    assert "javascript:" not in html
    assert "<script src" not in html


def test_bracket_modified_reference_glyph_has_exact_hover_identity():
    data = tables(sequence="A[m5C]GUACGU")
    model = build_coverage([("sample", data)], "rna")
    svg = render_coverage_svg(model)
    assert "Position 2: [m5C]" in svg
    assert "Reference * = bracket-annotated nucleotide" in svg


def test_build_and_render_do_not_mutate_native_dataframes():
    data = tables()
    before = copy.deepcopy(data)
    model = build_coverage([("sample", data)], "rna")
    render_coverage_html(model)
    for section in before:
        pd.testing.assert_frame_equal(data[section], before[section])
