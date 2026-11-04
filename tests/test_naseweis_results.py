"""Regression tests for NASE result parsing, modified bases and coverage."""

import pandas as pd
import pytest

from src.naseweis_results import coverage_html, match_table, read_mztab, split_sequence


@pytest.fixture
def result_file(tmp_path):
    path = tmp_path / "result.mzTab"
    path.write_text(
        "MTD\tmzTab-version\t1.0.0\n"
        "NUH\taccession\tdescription\topt_sequence\n"
        "NUC\tRNA<script>\tnull\tA[Am]CGU\n\n"
        "COM\tA comment between sections\n"
        "OLH\tsequence\taccession\tstart\tend\n"
        "OLI\tA[Am]Cp\tRNA<script>\t1\t3\n"
        "OLI\tA[Am]Cp\tsecond-RNA\t1\t3\n"
        "OSH\tsequence\tretention_time\tcharge\texp_mass_to_charge\tsearch_engine_score[1]\tsearch_engine_score[2]\n"
        "OSM\tA[Am]Cp\t25.5\t-2\t500.2\t120\t0.01\n",
        encoding="utf-8",
    )
    return path


def test_mztab_sections_and_match_table(result_file):
    tables = read_mztab(result_file)
    assert list(tables) == ["NUC", "OLI", "OSM"]
    matches = match_table(tables)
    assert matches.loc[0, "Accession"] == "RNA<script>; second-RNA"
    assert matches.loc[0, "q-value (%)"] == 1.0
    assert matches.loc[0, "Charge"] == -2
    assert matches.loc[0, "Retention time (s)"] == 25.5


def test_unfiltered_output_has_no_q_value_column(result_file):
    tables = read_mztab(result_file)
    tables["OSM"]["search_engine_score[2]"] = -1.0
    assert "q-value (%)" not in match_table(tables)


def test_empty_sections(tmp_path):
    path = tmp_path / "empty.mzTab"
    path.write_text(
        "NUH\taccession\topt_sequence\nOLH\tsequence\taccession\tstart\tend\nOSH\tsequence\tcharge\n"
    )
    tables = read_mztab(path)
    assert tables["NUC"].empty
    assert match_table(tables).empty


def test_malformed_row_reports_line(tmp_path):
    path = tmp_path / "invalid.mzTab"
    path.write_text("OSH\tsequence\tcharge\nOSM\tACG\n")
    with pytest.raises(ValueError, match="line 2"):
        read_mztab(path)


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("", []),
        ("ACGU", list("ACGU")),
        ("pA[Am]C[m5C]Gp", ["pA", "[Am]", "C", "[m5C]", "Gp"]),
    ],
)
def test_modified_sequence(sequence, expected):
    assert split_sequence(sequence) == expected


@pytest.mark.parametrize("sequence", ["A[Am", "A[[Am]]", "p", "<script>"])
def test_invalid_sequence(sequence):
    with pytest.raises(ValueError):
        split_sequence(sequence)


def test_coverage_escapes_accessions_and_handles_modifications(result_file):
    html = coverage_html(read_mztab(result_file), "RNA<script>")
    assert "60.0% sequence coverage" in html
    assert "RNA&lt;script&gt;" in html
    assert "<script>" not in html
    assert "A[Am]Cp" in html
    assert "positions 1–3" in html
    assert "http" not in html


def test_coverage_rejects_out_of_range_coordinates(result_file):
    tables = read_mztab(result_file)
    tables["OLI"].loc[0, "end"] = 100
    with pytest.raises(ValueError, match="exceed the reference"):
        coverage_html(tables, "RNA<script>")


def test_cache_is_scoped_to_full_path_and_revision(result_file, tmp_path):
    from src.Workflow import load_results

    stat = result_file.stat()
    first = load_results(str(result_file), stat.st_mtime_ns, stat.st_size)
    second_path = tmp_path / "other-workspace.mzTab"
    second_path.write_text(result_file.read_text().replace("500.2", "600.2"))
    # Same file size and deliberately identical revision: the path still differs.
    second = load_results(str(second_path), stat.st_mtime_ns, stat.st_size)
    assert first["OSM"].loc[0, "exp_mass_to_charge"] == 500.2
    assert second["OSM"].loc[0, "exp_mass_to_charge"] == 600.2


@pytest.mark.parametrize(
    "sequence,expected",
    [
        ("A[Am]Cc", ["A", "[Am]", "Cc"]),
        ("[5'-p]ACGU[3'-p]", ["pA", "C", "G", "Up"]),
    ],
)
def test_terminal_modifications(sequence, expected):
    assert split_sequence(sequence) == expected


def test_long_coverage_wraps_into_blocks():
    tables = {
        "NUC": pd.DataFrame({"accession": ["rna"], "opt_sequence": ["A" * 60]}),
        "OLI": pd.DataFrame(
            {"sequence": ["AAAA"], "accession": ["rna"], "start": [49], "end": [52]}
        ),
        "OSM": pd.DataFrame({"sequence": ["AAAA"]}),
    }
    plot = coverage_html(tables, "rna")
    assert "Positions 1–50" in plot
    assert "Positions 51–60" in plot
    assert plot.count('class="hit"') == 2
    assert "6.7% sequence coverage" in plot


def test_cached_results_refresh_after_an_in_place_update(result_file):
    from src.Workflow import load_results

    stat = result_file.stat()
    old = load_results(str(result_file), stat.st_mtime_ns, stat.st_size)
    result_file.write_text(result_file.read_text().replace("500.2", "700.2"))
    new = load_results(str(result_file), stat.st_mtime_ns + 1, stat.st_size)
    assert old["OSM"].loc[0, "exp_mass_to_charge"] == 500.2
    assert new["OSM"].loc[0, "exp_mass_to_charge"] == 700.2


def test_mutating_a_cached_table_does_not_change_other_views(result_file):
    from src.Workflow import load_results

    stat = result_file.stat()
    args = (str(result_file), stat.st_mtime_ns, stat.st_size)
    view = load_results(*args)
    view["OSM"].loc[0, "sequence"] = "CHANGED"
    assert load_results(*args)["OSM"].loc[0, "sequence"] == "A[Am]Cp"


@pytest.mark.parametrize("sequence", [None, ""])
def test_missing_reference_sequence_has_an_explanation(sequence):
    tables = {"NUC": pd.DataFrame({"accession": ["rna"], "opt_sequence": [sequence]})}
    assert "No reference sequence was exported" in coverage_html(tables, "rna")


def test_mixed_missing_q_values_are_not_reported_as_negative_percentages():
    tables = {
        "OSM": pd.DataFrame(
            {
                "sequence": ["ACG", "ACG", "ACG"],
                "search_engine_score[2]": [0.0, -1, float("nan")],
            }
        )
    }
    result = match_table(tables)
    assert result.loc[0, "q-value (%)"] == 0.0
    assert result["q-value (%)"].iloc[1:].isna().all()
    assert len(result) == 3


@pytest.mark.parametrize("sequence", ["c", "AcG", "ApCG"])
def test_invalid_terminal_markers_are_rejected(sequence):
    with pytest.raises(ValueError):
        split_sequence(sequence)


def test_legacy_coverage_imports_still_work(result_file):
    from src import plot_full_coverage as legacy

    tables = legacy.read_mzTab(result_file)
    assert legacy.split_sequence("A[Am]Cp") == ["A", "[Am]", "Cp"]
    # Both accessions need a reference sequence for a complete export.
    tables["NUC"] = tables["NUC"].iloc[:1]
    assert "60.0% sequence coverage" in legacy.make_coverage_html(tables)


def test_reference_without_any_hits_has_zero_coverage():
    tables = {
        "NUC": pd.DataFrame({"accession": ["rna"], "opt_sequence": ["ACGU"]}),
        "OLI": pd.DataFrame(columns=["accession", "sequence", "start", "end"]),
    }
    assert "0.0% sequence coverage" in coverage_html(tables, "rna")
