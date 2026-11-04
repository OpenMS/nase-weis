"""Upload your data."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup()
st.title("Upload your data")
st.write("Add an mzML run and an RNA FASTA database to this workspace.")
Workflow().show_file_upload_section()
