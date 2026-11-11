"""Upload your data."""

import os

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Upload your data")
st.write("Add MS/MS acquisitions and an RNA FASTA database to this workspace.")
st.caption(
    "Current files lists all available workspace inputs, including PRIDE downloads. "
    "Choose which acquisitions to analyze on Configure; adding a file does not select it for analysis."
)
if data_folder := os.environ.get("LOCAL_DATA_HOST_DIR"):
    st.info(
        f"Your data folder: `{data_folder}`. Move or copy your input files there, "
        "then click **Refresh folder**, select files and click **Add selected** below. "
        "The app reads these files in place; keep them in this folder to view saved runs. "
        "You can also upload files from your computer."
    )
Workflow().show_file_upload_section()
st.page_link("content/naseweis.py", label="Select samples for analysis", icon=":material/tune:")
