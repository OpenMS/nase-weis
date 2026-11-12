"""Explore uploaded acquisitions with the template's pyOpenMS-viz viewers."""

import streamlit as st

from src.common.common import page_setup
from src.naseweis_viewers import render_raw_viewer
from src.Workflow import Workflow

page_setup()
st.title("Explore your acquisition")
st.write(
    "Browse MS1 signal, individual MS1/MS2 spectra and chromatograms. Use Results → Spectra for matched RNA fragment assignments, or Results → Quantification for LFQ feature evidence."
)
st.page_link(
    "content/topp_workflow_results.py",
    label="View identified spectra and assigned RNA ions in Results",
    icon=":material/analytics:",
)
workflow = Workflow()
render_raw_viewer(
    workflow.workflow_dir,
    workflow.params,
    st.session_state.get("settings", {}).get("online_deployment", False),
)
