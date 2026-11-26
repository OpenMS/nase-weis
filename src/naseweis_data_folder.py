"""Where the data folder of a container installation is, for Describe experiment."""

import os

import streamlit as st

from src.workflow.host_paths import folder_file_uri, folder_reference


def data_folder_notice() -> None:
    """Say where to put files that 📂 Browse reads in place (LOCAL_DATA_HOST_DIR)."""
    data_folder = os.environ.get("LOCAL_DATA_HOST_DIR")
    if not data_folder:
        return
    st.info(
        f"Your data folder: {folder_reference(data_folder)}. Move or copy your input files there, "
        "then open **📂 Browse**, click **Refresh folder**, select files and click **Add selected**. "
        "The app reads these files in place; keep them in this folder to view saved runs. "
        "You can also upload files from your computer in the same dialog."
    )
    if folder_file_uri(data_folder):
        st.caption(
            "If the link does not open your file browser, copy the path below and "
            "paste it into the file browser's address or Go to Folder field."
        )
    else:
        st.caption(
            "Copy this folder path. Relative paths such as ./data are relative "
            "to the folder containing docker-compose.yml."
        )
    st.code(data_folder, language=None)
