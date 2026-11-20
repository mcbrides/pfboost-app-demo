"""Per-viewer appearance for the pinned Streamlit 1.65 frontend.

The native frontend theme also colors canvas dataframes and widget popovers.
The same-window host theme message changes only this browser session; no global
Python theme option or data-collection file is changed by the switch.
"""
from __future__ import annotations

import json
import streamlit as st

LIGHT = {
    "background": "#ffffff", "sidebar": "#f1f6fa", "surface": "#f1f7fb",
    "secondary": "#f1f6fa", "text": "#202f3b", "muted": "#5d7080",
    "heading": "#173e5b", "accent": "#327ba8", "link": "#256d9c",
    "border": "#dce6ee", "grid": "#eaf0f5", "profile": "#202a33",
    "viewer": "#ffffff", "hover": "#edf5fa",
    "compare": ["#173e5b", "#c9952f", "#589bc3", "#202a33", "#408b85", "#b77536"],
}
DARK = {
    "background": "#242629", "sidebar": "#2d3034", "surface": "#2e3237",
    "secondary": "#34383e", "text": "#e5e9ef", "muted": "#b0bac6",
    "heading": "#b9def4", "accent": "#78b8de", "link": "#96cfee",
    "border": "#484e56", "grid": "#3b4047", "profile": "#e5e9ef",
    "viewer": "#3c4147", "hover": "#3c4249",
    "compare": ["#9bc9e4", "#e5b649", "#6ba9ce", "#e5e9ef", "#78c2b4", "#dea86d"],
}


def is_dark():
    return bool(st.session_state.get("dark_mode", False))


def colors(dark=None):
    return DARK if (is_dark() if dark is None else dark) else LIGHT


def native_theme(dark):
    palette = colors(dark)
    return {
        "base": 1 if dark else 0,
        "primaryColor": palette["accent"],
        "backgroundColor": palette["background"],
        "secondaryBackgroundColor": palette["secondary"],
        "textColor": palette["text"], "linkColor": palette["link"],
        "borderColor": palette["border"],
        "dataframeHeaderBackgroundColor": palette["surface"],
        "codeBackgroundColor": palette["surface"], "showWidgetBorder": True,
        "sidebar": {"backgroundColor": palette["sidebar"],
                    "secondaryBackgroundColor": palette["secondary"],
                    "textColor": palette["text"], "borderColor": palette["border"]},
    }


def apply_appearance():
    dark = is_dark()
    palette = colors(dark)
    variables = ";".join(f"--pf-{key}:{value}" for key, value in palette.items()
                         if isinstance(value, str))
    css = """
    .block-container {max-width:1450px;padding-top:2.2rem;}
    h1 {letter-spacing:-.035em;color:var(--pf-heading);}
    h2,h3 {letter-spacing:-.02em;color:var(--pf-heading);}
    div[data-testid="stMetric"] {background:var(--pf-surface);border:1px solid var(--pf-border);border-radius:14px;padding:16px 20px;}
    div[data-testid="stMetricValue"] {color:var(--pf-heading);}
    section[data-testid="stSidebar"] {background:var(--pf-sidebar);border-right:1px solid var(--pf-border);}
    .pf-brand {display:flex;align-items:center;gap:12px;margin:0 0 28px;}
    .pf-brand strong {display:block;font-size:1.65rem;letter-spacing:-.055em;color:var(--pf-heading);line-height:1.2;}
    .pf-brand span {display:block;font-size:.78rem;color:var(--pf-muted);margin-top:4px;}
    .pf-page-header {display:flex;align-items:center;justify-content:space-between;gap:24px;margin-bottom:24px;padding:4px 0 22px;border-bottom:1px solid var(--pf-border);}
    .pf-page-header h1 {font-size:2.5rem;line-height:1.16;margin:8px 0 10px;padding:0;font-weight:700;color:var(--pf-heading);}
    .pf-page-header p {color:var(--pf-muted);font-size:.98rem;line-height:1.6;margin:0;max-width:760px;}
    .pf-eyebrow {color:var(--pf-accent);font-size:.68rem;font-weight:650;letter-spacing:.14em;}
    .pf-header-mark {flex-shrink:0;}
    @media(max-width:700px) {.pf-page-header h1 {font-size:1.9rem;} .pf-header-mark {width:66px;height:66px;} .pf-page-header {gap:12px;}}
    .entry-label {color:var(--pf-accent);font-size:.78rem;font-weight:600;letter-spacing:.09em;text-transform:uppercase;}
    .metadata {font-size:.95rem;margin:.3rem 0 .9rem;}
    .metadata-label {color:var(--pf-muted);font-size:.8rem;margin-bottom:.15rem;}
    .pf-colorbar {border:1px solid var(--pf-border);}
    .pf-data-status {color:var(--pf-muted);font-size:.75rem;line-height:1.5;}
    .pf-data-status summary {width:fit-content;cursor:pointer;padding:0;font-weight:400;}
    .pf-data-status summary:hover {color:var(--pf-text);}
    .pf-data-status summary:focus-visible {outline:2px solid var(--pf-accent);outline-offset:3px;}
    .pf-data-status ul {margin:6px 0 0;padding-left:18px;}
    .pf-data-status li {margin:4px 0;overflow-wrap:anywhere;white-space:pre-wrap;}
    .pf-data-status p {font-size:inherit;margin:6px 0 0;}
    .pf-updates {display:grid;grid-template-columns:auto minmax(0,1fr);gap:4px 18px;border-top:1px solid var(--pf-border);padding-top:8px;margin-top:12px;}
    .pf-updates-title {color:var(--pf-muted);font-size:.75rem;font-weight:600;line-height:1.5;padding:2px 0;}
    .pf-update-row {display:flex;align-items:baseline;gap:8px;padding:2px 0;line-height:1.5;}
    .pf-update-dot {width:4px;height:4px;border-radius:50%;background:var(--pf-accent);flex-shrink:0;align-self:center;}
    .pf-update-copy {min-width:0;flex:1;font-size:.75rem;}
    .pf-update-copy strong {font-weight:500;color:var(--pf-text);}
    .pf-update-copy span {color:var(--pf-muted);margin-left:8px;overflow-wrap:anywhere;}
    .pf-update-row time {font-size:.7rem;color:var(--pf-muted);white-space:nowrap;}
    @media(max-width:700px) {.pf-updates {grid-template-columns:minmax(0,1fr);} .pf-update-row {flex-wrap:wrap;} .pf-update-copy {flex-basis:70%;} .pf-update-row time {margin-left:12px;}}
    """
    st.html(f"<style>:root {{{variables};color-scheme:{'dark' if dark else 'light'};}}{css}</style>")
    # Streamlit's host protocol is version-bound and browser-local. Allow-listed
    # app origins are configured in .streamlit/config.toml.
    payload = json.dumps({"stCommVersion": 1, "type": "SET_CUSTOM_THEME_CONFIG",
                          "themeName": "PF database dark" if dark else "PF database light",
                          "themeInfo": native_theme(dark)}, separators=(",", ":"))
    st.html("""<script>(function () {
      const message = __MESSAGE__;
      const signature = JSON.stringify(message.themeInfo);
      if (window.__pfdbThemeSignature === signature) return;
      window.__pfdbThemeSignature = signature;
      const apply = () => window.postMessage(message, window.location.origin);
      apply();
      // The host-message listener may still be initializing on a first load.
      setTimeout(() => {if (window.__pfdbThemeSignature === signature) apply();}, 200);
    })();</script>""".replace("__MESSAGE__", payload), unsafe_allow_javascript=True)
