# Test-only: mounted over src/common/captcha_.py so a script can pass the captcha.
import streamlit as st
def captcha_control(*args, **kwargs):
    st.session_state["controllo"] = True
