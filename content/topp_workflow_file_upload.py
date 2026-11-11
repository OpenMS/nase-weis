"""Upload your data."""

import os

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Upload your data")
st.write("Add an mzML run and an RNA FASTA database to this workspace.")
if data_folder := os.environ.get("LOCAL_DATA_HOST_DIR"):
    st.info(
        f"Your data folder: `{data_folder}`. Move or copy your input files there, "
        "then click **Refresh folder**, select files and click **Add selected** below. "
        "The app reads these files in place; keep them in this folder to view saved runs. "
        "You can also upload files from your computer."
    )
Workflow().show_file_upload_section()
