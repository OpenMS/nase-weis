"""Standalone database preparation with the complete DecoyDatabase tool."""

import streamlit as st

from src.common.common import page_setup
from src.DecoyWorkflow import DecoyWorkflow

page_setup()
st.title("Build a target-decoy database")
st.write(
    "Prepare RNA or protein FASTA databases, including configurable shuffling, decoy ratios, protein entrapment and neighbor searches."
)
workflow = DecoyWorkflow()
# A step change must rerun the page: Streamlit tabs render every pane once,
# which otherwise leaves Configure's file list stale after a fragmented upload.
section = st.radio(
    "Database preparation",
    ["Upload", "Configure", "Run", "Results"],
    horizontal=True,
    key="decoy-database-step",
)
if section == "Upload":
    workflow.show_file_upload_section()
elif section == "Configure":
    workflow.show_parameter_section()
elif section == "Run":
    workflow.show_execution_section()
else:
    workflow.show_results_section()
