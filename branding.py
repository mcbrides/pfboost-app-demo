"""Small shared visual identity for the Streamlit interface."""

import base64
import html
from pathlib import Path

import streamlit as st

from appearance import is_dark

NAVY = "#173e5b"
BLUE = "#327ba8"
LIGHT_BLUE = "#9bc9e4"
GOLD = "#c9952f"
INK = "#202a33"
COMPARE_COLORS = [NAVY, GOLD, "#589bc3", INK, "#408b85", "#b77536"]

MARK_PATH = Path(__file__).resolve().parent / "assets" / "pfboost_logo_YD_LBC_l.svg"
DARK_MARK_PATH = MARK_PATH.with_name("pfboost_logo_YDa_LBC_d.svg")
MARK_URI = "data:image/svg+xml;base64," + base64.b64encode(MARK_PATH.read_bytes()).decode("ascii")

DARK_MARK_URI = "data:image/svg+xml;base64," + base64.b64encode(DARK_MARK_PATH.read_bytes()).decode("ascii")


def mark_uri():
    path = DARK_MARK_PATH if is_dark() else MARK_PATH
    try:
        return "data:image/svg+xml;base64," + base64.b64encode(path.read_bytes()).decode("ascii")
    except OSError:
        return DARK_MARK_URI if is_dark() else MARK_URI


def sidebar_brand():
    st.markdown(
        f'<div class="pf-brand"><img src="{mark_uri()}" alt="" width="52" height="52">'
        '<div><strong>PF database</strong><span>Protein protection factors</span></div></div>',
        unsafe_allow_html=True,
    )


def page_header(title: str, subtitle: str):
    st.markdown(
        '<div class="pf-page-header"><div><div class="pf-eyebrow">PROTEIN PROTECTION FACTORS</div>'
        f'<h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p></div>'
        f'<img class="pf-header-mark" src="{mark_uri()}" alt="Folded ribbon illustration" width="98" height="98"></div>',
        unsafe_allow_html=True,
    )
