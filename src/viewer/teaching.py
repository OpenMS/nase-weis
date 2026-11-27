"""What each view and component of Explore results shows, for its info icon.

One short paragraph per view (``VIEWS``, the preset's help) and per component
(``COMPONENTS``), written for someone who runs the search rather than for
someone who wrote it.
"""
from __future__ import annotations

VIEWS = {
    "summary": "How many spectra, oligonucleotides and RNAs this search found, and how the "
               "scores and precursor mass errors of its matches are distributed.",
    "matches": "Every spectrum match of the search. A click on a match shows its spectrum with "
               "the fragment ions NASE assigned, and those fragments as a table.",
    "spectra": "The acquisition itself: its scans, the picked spectrum, the MS1 chromatogram and "
               "peak map, with the identified spectra linked to them.",
    "coverage": "The RNAs of the database the search identified. Pick an RNA to see its "
                "oligonucleotides and which residues they cover, an oligonucleotide to see its "
                "spectrum matches, and a match to see its fragment ions.",
    "quantification": "MS1 feature areas of the identified oligonucleotides and their extracted "
                      "ion chromatograms.",
    "export": "The result files of the run: mzTab, idXML, the search database, the parameters "
              "and the run's log.",
    "identification": "Whether the search behaved: the number of matches, the distribution of "
                      "their scores and precursor mass errors, and the identified spectra.",
    "experiment_summary": "The acquisitions of the experiment, their conditions and status, and "
                          "what the search found across all of them.",
    "experiment_coverage": "The RNAs identified in the experiment. Pick an RNA in the table to see "
                           "its oligonucleotides and its coverage per condition or per acquisition.",
    "experiment_quantification": "The comparison of the conditions: the abundance of each "
                                 "oligonucleotide in both, their ratio and how it is spread. "
                                 "Change the conditions and the normalization with Edit comparison "
                                 "above.",
    "abundances": "The abundance of every oligonucleotide in every sample, the picked one by "
                  "condition, and the most variable ones as z-scores.",
    "samples": "The results of one acquisition of the experiment at a time.",
}

COMPONENTS = {
    "summary": "Spectrum matches that passed the FDR cutoff, the distinct oligonucleotides they "
               "name and the RNAs of the database these map to.",
    "matches": "The spectrum matches as a table, with a CSV download. Each row is one spectrum "
               "assigned to one oligonucleotide.",
    "spectrum_viewer": "The spectrum of the picked match with its theoretical fragments; matched "
                       "peaks are coloured by ion series.",
    "coverage": "The reference sequence of the picked RNA with every identified oligonucleotide "
                "drawn under the residues it covers. Colours count the spectra supporting it; "
                "modified residues are marked.",
    "rna_table": "The RNAs of the database the search identified, with their length, the "
                 "oligonucleotides and spectra found for them and the share of residues covered. "
                 "A click shows that RNA's oligonucleotides and coverage.",
    "oligo_table": "The oligonucleotides of the picked RNA, with where they lie in it, how many "
                   "spectra match them and the best score. A click lists their spectrum matches.",
    "oligo_matches": "The spectrum matches of the picked oligonucleotide, sorted by score. A click "
                     "shows the match's fragment ions.",
    "match_spectrum": "The spectrum of the picked match. Peaks NASE assigned to a fragment ion are "
                      "labelled and coloured by ion series; a click on a peak picks it in the "
                      "fragment table.",
    "fragment_table": "The fragment ions NASE assigned in the picked match, with their observed and "
                      "expected m/z and error. A click marks the peak in the fragment spectrum.",
    "quantification": "MS1 feature areas of the identified oligonucleotides, with the extracted "
                      "ion chromatogram of the picked one.",
    "downloads": "The result files of this run. mzTab holds every result table; idXML opens in "
                 "TOPPView.",
    "file_scans": "Every scan of the acquisition with its MS level, retention time and precursor. "
                  "A click shows that spectrum.",
    "file_spectrum": "The peaks of the picked spectrum. Zoom with the mouse; double-click resets.",
    "file_chromatogram": "Total ion and base peak chromatogram of the MS1 scans. A click picks the "
                         "nearest spectrum.",
    "file_peak_map": "Retention time against m/z of the MS1 peaks, coloured by intensity. Drag to "
                     "zoom.",
    "match_table": "The identified spectra, sorted by score. A click picks the match and its "
                   "spectrum in the other components.",
    "score_distribution": "How the search scores of the matches are spread. A long tail of low "
                          "scores points to random matches.",
    "precursor_error": "The difference between measured and calculated precursor m/z in ppm. A "
                       "narrow peak near zero means the mass tolerance fits the instrument.",
    "experiment_summary": "What the search found across the experiment's acquisitions, and the "
                          "status of each acquisition.",
    "experiment_coverage": "The reference sequence of the picked RNA with the oligonucleotides "
                           "identified in each condition (or acquisition) as tracks below it.",
    "experiment_quantification": "The comparison of two conditions: their mean abundance per "
                                 "oligonucleotide plotted against each other, and the ratios as "
                                 "a table.",
    "sample_details": "All results of one acquisition: summary, matches, spectra, coverage and "
                      "quantification.",
    "abundance_table": "The abundance of every oligonucleotide in every sample, normalized as set "
                       "in Edit comparison. A click picks it in the profile.",
    "abundance_profile": "The picked oligonucleotide's abundance in every sample, grouped by "
                         "condition, with each condition's mean.",
    "zscore_heatmap": "The most variable oligonucleotides, each scaled to its mean and spread "
                      "across the samples, so conditions that differ show as blocks of colour.",
}


def component(name: str) -> str | None:
    return COMPONENTS.get(name)


def view(name: str, kind: str = "") -> str | None:
    """The help of view ``name``; an experiment's view is looked up as ``experiment_<name>`` first."""
    if kind == "experiment" and f"experiment_{name}" in VIEWS:
        return VIEWS[f"experiment_{name}"]
    return VIEWS.get(name)
