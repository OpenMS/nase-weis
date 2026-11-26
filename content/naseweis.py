"""Describe experiment: the acquisitions, the RNA database and the search settings."""

import streamlit as st

from src.common.common import page_setup
from src.naseweis_data_folder import data_folder_notice
from src.Workflow import Workflow
from src.workflow import guide

page_setup()
wf = Workflow()

st.title("Describe experiment")
guide.page_intro(
    """
Choose your acquisitions and your RNA database below. The search settings start at values that suit a typical experiment.

- **MS/MS acquisitions**: one file per acquisition, mzML, Thermo *.raw* or a Bruker timsTOF *.d* folder
  (vendor files are converted for you). **📂 Browse** opens your files by folder and adds new ones.
- **RNA database**: the expected sequences as FASTA. Decoys are added for the false discovery rate
  unless the file already has them. A cached digest from an earlier search can be used instead.

No data at hand? **Load example dataset** adds the NME1 spectra with their database.
""",
    background="""
- **mzML** is the open file format for raw MS data. Other vendor formats convert to it with
  [ProteoWizard msConvert](https://proteowizard.sourceforge.io/) (keep *peak picking* on).
- **NucleicAcidSearchEngine** (NASE) digests the RNA sequences with the chosen RNase, adds the
  modifications and matches the oligonucleotides to the MS/MS spectra.
- **Decoys** are shuffled sequences that cannot be present; how many of them match estimates the
  false discovery rate, reported per match as a **q-value**.
- **Label-free quantification** extracts the MS1 signal of every identified oligonucleotide in each
  acquisition and compares the abundances of modified RNAs between conditions.
""",
    background_title="Terms explained",
)

wf.analysis_section()
data_folder_notice()
wf.show_parameter_section()

# Opens Run analysis and starts the analysis there at once
wf.show_run_button(
    label="Run analysis",
    summary=wf.run_bar_summary(),
    name_label="Experiment name",
    help="Save these settings, open Run analysis and start the search of the selected acquisitions.",
)
