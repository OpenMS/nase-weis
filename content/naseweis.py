"""Search configuration, using the current OpenMS parameter widgets."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Configure your search")
st.write(
    "Select the input files, digestion and mass accuracy. Parameters are saved automatically."
)
Workflow().show_parameter_section()
