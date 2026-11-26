"""Run analysis: start the search and follow its progress."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow
from src.workflow import guide

page_setup()
wf = Workflow()

st.title("Run the analysis")
guide.page_intro(
    "Start the analysis and follow its progress. Every selected acquisition is searched in turn; "
    "with **label-free quantification** on, the samples are then compared. You can close this page "
    "and come back later with the same workspace.\n\n"
    "Each run is kept as an experiment on **Explore results**."
)
st.info(
    f"{wf.analysis_summary_text()} Change it on **Describe experiment**.",
    icon="🧭",
)

wf.show_execution_section()

guide.next_steps(
    wf,
    (guide.EXPLORE, "Look at the matches, spectra, coverage and abundances of each run."),
    (guide.EXPORT, "Take the result files and the processing report of every run with you."),
)
