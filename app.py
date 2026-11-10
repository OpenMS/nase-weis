"""NASEWEIS: an OpenMS web app for oligonucleotide identification."""

import json
from pathlib import Path

import streamlit as st

if "settings" not in st.session_state:
    st.session_state.settings = json.loads(
        Path("settings.json").read_text(encoding="utf-8")
    )

pages = {
    "NASEWEIS": [
        st.Page("content/quickstart.py", title="Overview", icon=":material/home:"),
    ],
    "RNA search": [
        st.Page(
            "content/topp_workflow_file_upload.py",
            title="Upload",
            icon=":material/upload_file:",
        ),
        st.Page("content/naseweis.py", title="Configure", icon=":material/tune:"),
        st.Page(
            "content/topp_workflow_execution.py",
            title="Run",
            icon=":material/play_arrow:",
        ),
        st.Page(
            "content/topp_workflow_results.py",
            title="Results",
            icon=":material/analytics:",
        ),
        st.Page(
            "content/naseweis_raw_viewer.py",
            title="Raw data viewer",
            icon=":material/query_stats:",
        ),
    ],
    "Database tools": [
        st.Page(
            "content/decoy_database.py",
            title="DecoyDatabase",
            icon=":material/genetics:",
        )
    ],
}
st.navigation(pages).run()
