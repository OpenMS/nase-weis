"""Explore results."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow
from src.workflow import guide

page_setup()
wf = Workflow()

st.title("Explore results")
guide.page_intro(
    "Choose a saved analysis to compare sample abundances and inspect each acquisition's "
    "matches, spectra and sequence coverage."
)
wf.show_results_section()

guide.next_steps(
    wf,
    (guide.EXPORT, "Take the result files and the processing report of a run with you."),
    (guide.DESCRIBE, "Search more acquisitions, or the same ones with other settings."),
)
