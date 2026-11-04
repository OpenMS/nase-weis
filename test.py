"""Validate the bundled RNA example without downloads or synthetic delays."""

from pathlib import Path

import pyopenms as poms


def test_example_contains_msms_spectra():
    experiment = poms.MSExperiment()
    poms.MzMLFile().load("example-data/mzML/NME1.mzML", experiment)
    assert sum(spectrum.getMSLevel() == 2 for spectrum in experiment) == 2116


def test_example_rna_database():
    lines = Path("example-data/fasta/NME1_Target.fasta").read_text().splitlines()
    assert lines[0] == ">NME1"
    sequence = poms.NASequence.fromString("".join(lines[1:]))
    assert sequence.size() == 340
