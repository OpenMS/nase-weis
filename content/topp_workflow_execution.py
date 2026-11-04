"""Run the RNA search."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Run the RNA search")
st.write("Name this run and start the search. Progress and tool logs appear below.")
Workflow().show_execution_section()
