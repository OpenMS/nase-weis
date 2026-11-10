"""Explore your results."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Explore your results")
st.write(
    "Choose a saved run to inspect matches, sequence coverage, optional MS1 quantification and downloadable results."
)
Workflow().show_results_section()
