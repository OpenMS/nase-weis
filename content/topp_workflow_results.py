"""Explore your results."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Explore your results")
st.write(
    "Choose a saved analysis to compare sample abundances and inspect each acquisition's matches, spectra and sequence coverage."
)
Workflow().show_results_section()
