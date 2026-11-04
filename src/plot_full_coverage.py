"""Compatibility imports for NASEWEIS's original coverage module.

Authors: Hendrik Weisser, Samuel Wein.
License: BSD 3-clause.
"""

from src.naseweis_results import coverage_html, read_mztab, split_sequence

read_mzTab = read_mztab


def make_coverage_html(mztab1):
    """Render the accessions in one NASE mzTab file, entirely offline."""
    return "\n".join(
        coverage_html(mztab1, accession) for accession in mztab1["NUC"]["accession"]
    )
