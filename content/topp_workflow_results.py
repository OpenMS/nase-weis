"""Explore your results."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Explore your results")
st.write(
    "Choose a completed run to inspect matches, sequence coverage and downloadable results."
)
Workflow().show_results_section()
