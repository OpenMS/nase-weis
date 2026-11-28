"""NASEWEIS: an OpenMS web app for oligonucleotide identification."""

import json
from pathlib import Path

import streamlit as st

from src.Workflow import Workflow
from src.workflow.guide import DESCRIBE, EXPLORE, EXPORT, RUN, step_page

if "settings" not in st.session_state:
    st.session_state.settings = json.loads(
        Path("settings.json").read_text(encoding="utf-8")
    )

pages = {
    "NASEWEIS": [
        st.Page("content/quickstart.py", title="Overview", icon=":material/home:"),
    ],
    "RNA search": [
        step_page(Workflow.configure_page, DESCRIBE),
        step_page(Workflow.run_page, RUN),
        step_page(Workflow.results_page, EXPLORE),
        step_page(Workflow.export_page, EXPORT),
    ],
    "Tools": [
        st.Page(
            "content/naseweis_raw_viewer.py",
            title="Raw data viewer",
            icon=":material/query_stats:",
        ),
        st.Page(
            "content/decoy_database.py",
            title="DecoyDatabase",
            icon=":material/genetics:",
        )
    ],
}
st.navigation(pages).run()
