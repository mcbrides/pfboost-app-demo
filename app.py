"""PFBoost’s read-only Streamlit browser for the bundled data snapshot.

Interface organization adapted from the lab's PFDB example.
"""

from __future__ import annotations

import hashlib
import html
from datetime import datetime
from functools import partial
from pathlib import Path
from urllib.parse import quote
from zoneinfo import ZoneInfo

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from branding import BLUE, GOLD, LIGHT_BLUE, mark_uri, page_header, sidebar_brand
from appearance import apply_appearance, colors, is_dark
from structure_viewer import structure_html
from state_labels import STATE_LABELS, display_binding
from sequence_io import sequence_records, fasta_bytes

from data_access import (
    PF_COLOR_MIN, PF_COLOR_MAX, PF_COLORSCALE, UNMAPPED_COLOR,
    as_dict, catalog_display, clear_caches, default_data_root, entry_folder, experimental_method_label, load_catalog,
    load_entry, mapped_colors, number, read_structure, recent_index_updates, select_pdb_model,
    structure_files, summarize, text,
)

st.set_page_config(page_title="PF database", page_icon=mark_uri(), layout="wide")
apply_appearance()

ROOT = default_data_root()
PLOT_CONFIG = {"displaylogo": False}


def label(entry_id: str, catalog: pd.DataFrame) -> str:
    row = catalog.loc[catalog.entry_id == entry_id].iloc[0]
    temperature = number(row.temperature_c)
    condition = f" · {temperature:.1f} °C" if temperature is not None else ""
    variant = f" · {row.variant_label}" if text(row.variant_label) else ""
    state = f" · {row.state_label}" if text(row.state_label) else ""
    return f"{entry_id} · {row.protein}{variant}{state}{condition}"


def metadata_field(label_text: str, value):
    rendered = html.escape(text(value) or "Not recorded")
    st.markdown(
        f'<div class="metadata"><div class="metadata-label">{html.escape(label_text)}</div>{rendered}</div>',
        unsafe_allow_html=True,
    )


def show_recent_updates(updates):
    labels = {"added": ("Dataset added", "Datasets added"),
              "structure": ("Structure updated", "Structures updated"),
              "pf": ("Measurement entry updated", "Measurement entries updated"),
              "entry": ("Entry updated", "Entries updated")}
    updates = [event for event in updates if event.get("kind") in labels and event.get("count", 0) > 0]
    if not updates:
        return
    rows = []
    for event in updates[:2]:
        count = event.get("count", 0)
        singular, plural = labels[event["kind"]]
        title = singular if count <= 1 else f"{count} {plural.lower()}"
        entries = event.get("entries", [])
        if count == 1 and entries:
            entry = entries[0]
            detail = " · ".join(filter(None, [entry.get("id"), entry.get("protein")]))
        else:
            detail = ", ".join(entry.get("id", "") for entry in entries)
            if count > len(entries):
                detail += f" +{count-len(entries)}"
        stamp = datetime.fromtimestamp(event["at"], ZoneInfo("America/New_York"))
        rows.append(
            '<div class="pf-update-row"><span class="pf-update-dot" aria-hidden="true"></span>'
            f'<div class="pf-update-copy"><strong>{html.escape(title)}</strong>'
            + (f'<span>{html.escape(detail)}</span>' if detail else '')
            + f'</div><time datetime="{stamp.isoformat()}">{stamp.strftime("%b %d · %H:%M %Z")}</time></div>'
        )
    st.markdown('<section class="pf-updates" aria-label="Recent updates">'
                '<div class="pf-updates-title">Recent updates</div><div class="pf-updates-list">'
                + ''.join(rows) + '</div></section>', unsafe_allow_html=True)


def show_data_status(issues):
    if not issues:
        return
    count = len(issues)
    label = f"Data status · {count} {'notice' if count == 1 else 'notices'}"
    messages = ''.join(f'<li>{html.escape(str(issue))}</li>' for issue in issues[:20])
    remaining = f'<p>Showing 20 of {count} notices.</p>' if count > 20 else ''
    st.markdown(
        '<details class="pf-data-status">'
        f'<summary>{label}</summary><ul>{messages}</ul>{remaining}</details>',
        unsafe_allow_html=True,
    )


def catalog_numeric_columns():
    return {
        "pH": st.column_config.NumberColumn(width=58, format="%.2f", alignment="left", help="Numeric pH_input used to calculate the standardized PF values."),
        "Temp": st.column_config.NumberColumn(width=58, format="%.1f", alignment="left", help="Exchange temperature in °C."),
        "Sites": st.column_config.NumberColumn(width=60, format="%d", alignment="left", help="Unique construct positions with retained PF measurements."),
    }


def format_plot(figure):
    palette = colors()
    figure.update_layout(
        template="plotly_dark" if is_dark() else "plotly_white",
        font={"family": "sans-serif", "color": palette["muted"]},
        margin={"l": 15, "r": 20, "t": 30, "b": 20}, legend_title_text="",
        hovermode="closest", colorway=palette["compare"],
        paper_bgcolor=palette["background"], plot_bgcolor=palette["background"],
        hoverlabel={"bgcolor": palette["surface"], "font_color": palette["text"],
                    "bordercolor": palette["border"]},
    )
    figure.update_xaxes(showgrid=False, color=palette["muted"])
    figure.update_yaxes(gridcolor=palette["grid"], zerolinecolor=palette["border"], color=palette["muted"])
    return figure


def plot_profile(values: pd.DataFrame, show_published=False):
    figure = go.Figure()
    for column, name, color in [
        ("log10_pf_standard", "Standardized PF", colors()["profile"]),
        ("log10_pf_reported", "Published PF", GOLD),
    ]:
        frame = values.copy()
        if column == "log10_pf_reported":
            if not show_published or "pf_reported" not in frame:
                continue
            import math
            frame[column] = frame.pf_reported.map(
                lambda value: math.log10(value) if number(value) is not None and value > 0 else float("nan")
            )
        if column not in frame:
            continue
        frame = frame.dropna(subset=["sequence_index", column])
        if frame.empty:
            continue
        # Markers preserve gaps where no measurement exists.
        figure.add_trace(go.Scatter(
            x=frame.sequence_index, y=frame[column], mode="markers", name=name,
            marker={"color": color, "size": 7, "opacity": .9,
                    "symbol": "diamond-open" if column == "log10_pf_reported" else "circle"},
            customdata=frame[["residue"]].to_numpy(),
            hovertemplate="%{customdata[0]}%{x}<br>log₁₀(PF): %{y:.3f}<extra>%{fullData.name}</extra>",
        ))
    if not figure.data:
        st.info("No numerical protection factors are recorded for this entry.")
        return
    figure.update_layout(height=365, xaxis_title="Construct residue position", yaxis_title="log₁₀(PF)")
    st.plotly_chart(format_plot(figure), width="stretch", config=PLOT_CONFIG)


def show_pf_legend():
    stops = ", ".join(f"{color} {position * 100:g}%" for position, color in PF_COLORSCALE)
    ticks = "".join(f"<span>{value}</span>" for value in [0, 2, 4, 6, 8, 10])
    st.markdown(
        '<div style="max-width:520px;margin:0.4rem 0 0.8rem">'
        '<div style="font-size:0.85rem;margin-bottom:5px">log₁₀(PF)</div>'
        f'<div class="pf-colorbar" role="img" aria-label="White to dark blue, log10 PF zero to ten" '
        f'style="height:15px;border-radius:3px;background:linear-gradient(to right, {stops})"></div>'
        f'<div style="display:flex;justify-content:space-between;font-size:0.8rem;margin-top:3px">{ticks}</div>'
        f'<div style="font-size:0.85rem;margin-top:8px"><span style="display:inline-block;width:12px;height:12px;'
        f'background:{UNMAPPED_COLOR};border-radius:2px;margin-right:6px"></span>Gray: no mapped PF</div></div>',
        unsafe_allow_html=True,
    )
    st.caption("Same 0–10 scale for every protein. Values outside the scale use the endpoint colors.")


def show_structure(folder: Path, metadata: dict, values: pd.DataFrame, mapping: pd.DataFrame):
    candidates = structure_files(folder, metadata)
    if not candidates:
        st.info("No local structure file is available for this entry yet.")
        return
    if len(candidates) > 1:
        version = hashlib.sha256("|".join(p.name for p in candidates).encode()).hexdigest()[:10]
        selected = st.selectbox("Structure file", candidates, format_func=lambda p: p.name,
                                key=f"structure_file_{folder.name}_{version}")
    else:
        selected = candidates[0]
    if selected.name != "structure.pdb":
        st.caption(f"Reference structure · {selected.name}")
    if as_dict(metadata.get("structure")).get("representation") == "separate_domain_references":
        st.caption("Separate domain reference; no interdomain geometry is defined.")
    pdb_text = read_structure(selected)
    pdb_text, model = select_pdb_model(pdb_text, as_dict(metadata.get("structure")).get("selected_model"))
    styles, measured_count = mapped_colors(values, mapping, pdb_text, model, coordinate_file=selected.name)
    st.iframe(structure_html(pdb_text, styles, dark=is_dark()), height=485)
    show_pf_legend()
    if styles:
        st.caption(f"{measured_count} measured construct positions colored using the residue mapping · model {model}")
    else:
        st.caption("Showing the structure without PF coloring; no matching mapped PF values are available.")


def show_sequence(metadata, entry_id):
    st.markdown("**Analyzed protein sequence**")
    try:
        records = sequence_records(metadata, entry_id)
    except ValueError as exc:
        st.info("The saved sequence needs review before FASTA export.")
        st.caption(str(exc))
        return
    if not records:
        st.info("The exact experimental construct sequence has not been resolved for this entry.")
        review = as_dict(metadata.get("unmapped_construct_pf_review"))
        if review.get("unresolved_sequence_evidence"):
            st.caption(review["unresolved_sequence_evidence"])
        reference = text(review.get("numbering_reference_sequence"))
        if reference:
            with st.expander("Source numbering reference"):
                st.caption("Reference sequence used to identify reported sites; this is not the unresolved experimental construct sequence.")
                st.code(reference, language=None)
        return
    if metadata.get("sequence_definition"):
        st.caption(text(metadata["sequence_definition"]))
    else:
        st.caption("Saved analyzed sequence and numbering frame. A reviewed core may omit unreported terminal tags; the full qualifications remain in the metadata.")
    if len(records) > 1:
        st.caption("Separate polypeptide chains are exported as separate FASTA records. Their documented offsets link back to the measurements.")
    for record in records:
        entity = f" · {record['entity_id']}" if len(records) > 1 else ""
        st.markdown(f"{len(record['sequence'])} amino acids{entity}")
        st.code(fasta_bytes([record]).decode(), language=None, wrap_lines=True)
    st.download_button("Download sequence FASTA", fasta_bytes(records),
                       file_name=f"{entry_id}.fasta", mime="text/plain", key=f"sequence_download_{entry_id}")
    if metadata.get("construct"):
        metadata_field("Full construct definition", metadata["construct"])
    if metadata.get("construct_qualification"):
        metadata_field("Sequence scope / qualifications", metadata["construct_qualification"])


def show_entry(root: Path, entry_id: str, visual_only=False):
    try:
        metadata, values, mapping = load_entry(root, entry_id)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        st.info("This entry is being updated or its files are not yet complete. Refresh shortly.")
        with st.expander("File details"):
            st.code(str(exc), language=None)
        return
    summary = summarize(metadata)
    folder = entry_folder(root, entry_id)
    st.divider()
    st.markdown(f'<div class="entry-label">{html.escape(entry_id)}</div>', unsafe_allow_html=True)
    st.subheader(summary["protein"] or entry_id)
    if summary["protein_name"] and summary["protein_name"] != summary["protein"]:
        st.caption(summary["protein_name"])
    links = []
    for name, value, base in [
        ("PubMed", summary["pmid"], "https://pubmed.ncbi.nlm.nih.gov/"),
        ("UniProt", summary["uniprot"], "https://www.uniprot.org/uniprotkb/"),
        ("PDB baseline" if as_dict(metadata.get("structure")).get("baseline_reference")
            or as_dict(metadata.get("structure")).get("pdb_id_is_template_reference") else "PDB",
         summary["pdb_id"], "https://www.rcsb.org/structure/"),
        ("DOI", text(metadata.get("DOI")), "https://doi.org/"),
    ]:
        if value:
            links.append(f"[{name}: {value}]({base}{quote(value, safe='')})")
    if links:
        st.markdown(" · ".join(links))
    columns = st.columns([2, 2, 1])
    with columns[0]:
        metadata_field("Construct / variant", summary["variant_label"])
        if summary["expression_tag"]:
            metadata_field("Expression tag", summary["expression_tag"])
        if summary["modifications"]:
            metadata_field("Modifications / tags", summary["modifications"])
    with columns[1]:
        displayed_binding, binding_assumed = display_binding(metadata)
        if binding_assumed:
            metadata_field("Binding state", "Apo (assumed; binding not recorded)")
        else:
            metadata_field("Binding state", STATE_LABELS.get(displayed_binding, displayed_binding))
        if metadata.get("binding_partner_variant"):
            metadata_field("Partner variant", metadata["binding_partner_variant"])
        for key, title in [("conformation", "Conformation"), ("redox_state", "Redox state"), ("assembly", "Assembly")]:
            value = text(metadata.get(key))
            if value:
                metadata_field(title, STATE_LABELS.get(value, value))
        ph = summary["ph_input"]
        metadata_field("pH", f"{ph:.2f}" if ph is not None else None)
    with columns[2]:
        temperature = summary["temperature_c"]
        metadata_field("Temperature", f"{temperature:.1f} °C" if temperature is not None else None)
        metadata_field("Measured positions", values.sequence_index.nunique())
    if visual_only:
        mode = st.radio("Display", ["PF profile", "Histogram", "Sequence map", "3D structure"],
                        horizontal=True, key=f"mode_{entry_id}")
        show_visual(mode, folder, metadata, values, mapping)
        return
    profile_tab, sequence_tab, measurements_tab, structure_tab, details_tab = st.tabs(
        ["PF profile", "Sequence", "Measurements", "Structure", "Experiment details"]
    )
    with profile_tab:
        published = st.checkbox("Overlay published PF values", value=False, key=f"published_{entry_id}")
        plot_profile(values, published)
    with sequence_tab:
        show_sequence(metadata, entry_id)
    with measurements_tab:
        numeric_columns = {
            column: st.column_config.NumberColumn(alignment="left")
            for column in values.select_dtypes(include="number").columns
        }
        numeric_columns["sequence_index"] = st.column_config.NumberColumn(
            "residue_index", width=115, format="%d", alignment="left",
            help="Residue position in the analyzed construct sequence.")
        numeric_columns["residue"] = st.column_config.TextColumn(width=75, alignment="left")
        st.dataframe(values, hide_index=True, width="stretch",
                     column_config=numeric_columns)
        st.download_button("Download measurements", values.to_csv(index=False),
                           file_name=f"{entry_id}_residues.csv", mime="text/csv", key=f"download_{entry_id}")
    with structure_tab:
        # Structure files are loaded only on request, keeping routine browsing quick.
        if st.toggle("Show 3D structure", key=f"show_structure_{entry_id}"):
            try:
                show_structure(folder, metadata, values, mapping)
            except (OSError, ValueError) as exc:
                st.info(f"The structure file is being updated or cannot be read: {exc}")
        else:
            st.caption("Open the local structure and color it by protection factor.")
    with details_tab:
        conditions = as_dict(metadata.get("conditions"))
        metadata_field("Experimental method", experimental_method_label(metadata))
        metadata_field("Publication source", text(metadata.get("publication_source")))
        source = summary["structure_source"]
        pending_reason = text(as_dict(metadata.get("structure")).get("pending_reason"))
        if source == "No structure selected" and pending_reason:
            source += ". " + pending_reason
        metadata_field("Structure source", source)
        if summary["structure_preparation"]:
            metadata_field("Structure preparation", summary["structure_preparation"])
        if binding_assumed:
            metadata_field("Binding-state display assumption", "APO is the database display default when binding is unrecorded. This is an assumption; Native describes conformation and does not establish ligand absence.")
        if summary["construct"] != summary["variant_label"]:
            metadata_field("Full construct definition", summary["construct"])
        if metadata.get("construct_notes"):
            metadata_field("Construct notes", metadata["construct_notes"])
        if metadata.get("state_notes"):
            metadata_field("State notes", metadata["state_notes"])
        if metadata.get("state_source_label"):
            metadata_field("Source state description", metadata["state_source_label"])
        state_assignment = as_dict(metadata.get("state_assignment"))
        if state_assignment.get("note"):
            metadata_field("State assignment", state_assignment["note"])
        if summary["source_construct_id"]:
            metadata_field("Source construct ID", summary["source_construct_id"])
        if summary["construct_source_label"]:
            metadata_field("Original construct description", summary["construct_source_label"])
        if summary["ph_label"]:
            metadata_field("Source pH / pD note", summary["ph_label"])
        metadata_field("Buffer", conditions.get("buffer") or conditions.get("buffer_missing_reason")
                       or "Not specified in the available experimental metadata")
        fraction = number(conditions.get("d2o_fraction"))
        metadata_field("D₂O", f"{fraction * 100:g}%" if fraction is not None else None)
        metadata_field("Denaturant", conditions.get("denaturant"))
        metadata_field("Measured construct sequence", metadata.get("sequence"))
        if conditions.get("assumption_note"):
            metadata_field("Condition notes", conditions["assumption_note"])
        qualifications = metadata.get("qualifications")
        if qualifications:
            metadata_field("Measurement notes", qualifications)
        st.caption("Standardized and published quantities are kept in separate columns. Unrecorded values remain blank.")


def show_visual(mode, folder, metadata, values, mapping):
    if mode == "PF profile":
        plot_profile(values)
    elif mode == "Histogram":
        if "log10_pf_standard" not in values:
            st.info("No standardized PF values are available.")
            return
        frame = values.dropna(subset=["log10_pf_standard"])
        fig = px.histogram(frame, x="log10_pf_standard", nbins=25,
                           color_discrete_sequence=[LIGHT_BLUE], labels={"log10_pf_standard": "log₁₀(PF)"})
        fig.update_traces(marker_line_color=BLUE, marker_line_width=.6)
        fig.update_layout(height=365, bargap=.08, yaxis_title="Residue measurements")
        st.plotly_chart(format_plot(fig), width="stretch", config=PLOT_CONFIG)
    elif mode == "Sequence map":
        sequence = text(metadata.get("sequence"))
        if not sequence:
            st.info("No construct sequence is recorded for this entry.")
            return
        sequence_frame = pd.DataFrame({"sequence_index": range(1, len(sequence) + 1), "residue": list(sequence)})
        if "log10_pf_standard" in values:
            if values.sequence_index.duplicated().any():
                st.info("This entry has multiple observations at a position; use the PF profile to inspect them.")
                return
            sequence_frame = sequence_frame.merge(values[["sequence_index", "log10_pf_standard"]],
                                                  on="sequence_index", how="left")
        else:
            sequence_frame["log10_pf_standard"] = float("nan")
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=sequence_frame.sequence_index, y=[0] * len(sequence_frame), mode="markers",
            marker={"symbol": "square", "size": 13, "color": UNMAPPED_COLOR},
            customdata=sequence_frame[["residue"]].to_numpy(), showlegend=False,
            hovertemplate="%{customdata[0]}%{x}<br>No numerical PF<extra></extra>",
        ))
        measured = sequence_frame.dropna(subset=["log10_pf_standard"])
        fig.add_trace(go.Scatter(
            x=measured.sequence_index, y=[0] * len(measured), mode="markers",
            marker={"symbol": "square", "size": 13, "color": measured.log10_pf_standard,
                    "colorscale": PF_COLORSCALE, "cmin": PF_COLOR_MIN, "cmax": PF_COLOR_MAX,
                    "cauto": False, "showscale": True,
                    "line": {"color": UNMAPPED_COLOR, "width": 0.25},
                    "colorbar": {"title": "log₁₀(PF)", "tickvals": [0, 2, 4, 6, 8, 10]}},
            customdata=measured[["residue", "log10_pf_standard"]].to_numpy(), showlegend=False,
            hovertemplate="%{customdata[0]}%{x}<br>log₁₀(PF): %{customdata[1]:.3f}<extra></extra>",
        ))
        fig.update_yaxes(visible=False)
        fig.update_layout(height=240, xaxis_title="Construct residue position")
        st.plotly_chart(format_plot(fig), width="stretch", config=PLOT_CONFIG)
        st.caption("Gray: no numerical PF. White → light blue → dark blue: log₁₀(PF) 0–10, fixed across all proteins. Uses the measured construct sequence.")
    else:
        try:
            show_structure(folder, metadata, values, mapping)
        except (OSError, ValueError) as exc:
            st.info(f"The structure file is being updated or cannot be read: {exc}")


def table_selected(table_key, identifiers):
    selection = st.session_state[table_key].selection.rows
    if selection and selection[0] < len(identifiers):
        st.session_state["browse_entry"] = identifiers[selection[0]]


def show_browse(catalog):
    page_header("Browse datasets", "Search proteins, constructs, publications, and identifiers. Select an experiment to explore its measurements.")
    search_col, protein_col, state_col, publication_col = st.columns([2, 1.2, 1.2, 1])
    query = search_col.text_input("Search", placeholder="Protein, mutation, PMID, UniProt, PDB…")
    proteins = sorted(set(catalog.protein) - {""})
    protein = protein_col.selectbox("Protein", ["All proteins"] + proteins)
    state = state_col.selectbox("State", ["All states"] + sorted(set(catalog.state_label) - {""}))
    publication = publication_col.selectbox(
        "Publication", ["All publications"] + sorted(set(catalog.publication) - {""})
    )
    filtered = catalog.copy()
    if query.strip():
        haystack = filtered[["entry_id", "protein", "protein_name", "construct", "construct_source_label", "expression_tag", "source_construct_id", "modifications", "pmid", "doi", "publication", "uniprot", "pdb_id", "structure_source", "state", "state_label"]].fillna("").agg(" ".join, axis=1)
        filtered = filtered[haystack.str.contains(query.strip(), case=False, regex=False)]
    if protein != "All proteins":
        filtered = filtered[filtered.protein == protein]
    if state != "All states":
        filtered = filtered[filtered.state_label == state]
    if publication != "All publications":
        filtered = filtered[filtered.publication == publication]
    st.caption(f"{len(filtered):,} of {len(catalog):,} datasets")
    if filtered.empty:
        st.info("No datasets match these filters.")
        return
    ids = filtered.entry_id.tolist()
    table_key = "browse_table_compact_" + hashlib.sha256("|".join(ids).encode()).hexdigest()[:12]
    st.dataframe(catalog_display(filtered), hide_index=True, width="stretch",
                 height=310, key=table_key, on_select=partial(table_selected, table_key, ids),
                 selection_mode="single-row", column_config={
                     **catalog_numeric_columns(),
                     "Entry": st.column_config.TextColumn(width=96, pinned=True),
                     "Protein": st.column_config.TextColumn(width=145),
                     "Variant": st.column_config.TextColumn(width=145, help="Mutations, domain, truncation or engineered construct."),
                     "State": st.column_config.TextColumn(width=150, help="Binding state and distinct conformational/redox states. Routine Native/Folded labels and assembly are kept in entry details."),
                     "PMID": st.column_config.TextColumn(width=84),
                     "UniProt": st.column_config.TextColumn(width=80),
                     "PDB": st.column_config.TextColumn(width=62),
                 })
    if st.session_state.get("browse_entry") not in ids:
        st.session_state["browse_entry"] = ids[0]
    selected = st.selectbox("Open dataset", ids, format_func=lambda eid: label(eid, filtered), key="browse_entry")
    show_entry(ROOT, selected)


def show_compare(catalog):
    page_header("Compare experiments", "Explore PF profiles across related constructs and experimental conditions.")
    groups = catalog.apply(lambda row: row.uniprot or row.protein, axis=1)
    group_labels = {}
    for group in sorted(set(groups) - {""}):
        subset = catalog[groups == group]
        if len(subset) >= 2:
            group_labels[group] = f"{subset.iloc[0].protein} · {group} · {len(subset)} entries"
    if not group_labels:
        st.info("No protein has multiple entries to compare yet.")
        return
    group = st.selectbox("Protein", list(group_labels), format_func=group_labels.get)
    subset = catalog[groups == group]
    ids = subset.entry_id.tolist()
    chosen = st.multiselect("Experiments", ids, default=ids[:2],
                            format_func=lambda eid: label(eid, subset), max_selections=6, key=f"compare_{group}")
    if not chosen:
        st.info("Choose one or more experiments.")
        return
    fig = go.Figure()
    sequences = set()
    for trace_index, eid in enumerate(chosen):
        try:
            metadata, values, _ = load_entry(ROOT, eid)
        except (OSError, ValueError, pd.errors.ParserError):
            st.caption(f"{eid} is being updated; refresh shortly.")
            continue
        sequences.add(text(metadata.get("sequence")))
        if "log10_pf_standard" not in values:
            continue
        frame = values.dropna(subset=["sequence_index", "log10_pf_standard"])
        fig.add_trace(go.Scatter(x=frame.sequence_index, y=frame.log10_pf_standard,
                                 mode="markers", name=label(eid, subset),
                                 marker={"size": 7, "color": colors()["compare"][trace_index],
                                         "symbol": ["circle", "diamond", "square", "triangle-up", "cross", "x"][trace_index]}))
    if len(sequences) > 1:
        st.info("These constructs have different sequences. Positions refer to each construct; they have not been aligned.")
    if fig.data:
        fig.update_layout(height=440, xaxis_title="Construct residue position", yaxis_title="log₁₀(PF)")
        st.plotly_chart(format_plot(fig), width="stretch", config=PLOT_CONFIG)
    st.dataframe(catalog_display(subset[subset.entry_id.isin(chosen)]), hide_index=True, width="stretch",
                 column_config=catalog_numeric_columns())


with st.sidebar:
    sidebar_brand()
    page = st.radio("Navigation", ["Home", "Browse", "Visualize", "Compare"], index=1)
    st.divider()
    auto_refresh = st.toggle("Automatic refresh", value=True, key="auto_refresh")
    st.toggle("Dark mode", key="dark_mode")
    if st.button("Refresh now", width="stretch"):
        clear_caches()
    st.caption("Checks for changes every 30 seconds while open." if auto_refresh else "Use Refresh now to load changes.")


@st.fragment(run_every="30s" if auto_refresh else None)
def render_page():
    try:
        with st.spinner("Loading the collection…"):
            catalog, issues = load_catalog(ROOT)
    except (OSError, ValueError):
        st.title("PF database")
        st.error("The bundled data folder is not accessible.")
        st.caption("Check that the data folder is present next to app.py.")
        return
    if catalog.empty:
        st.info("No complete entry folders have been found yet.")
        show_data_status(issues)
        return
    updates = recent_index_updates(ROOT)
    if page == "Home":
        page_header("The PF database", "Explore experimental protection factors, from individual residues to complete protein structures.")
        metrics = st.columns(3)
        metrics[0].metric("Datasets", f"{len(catalog):,}")
        metrics[1].metric("Parent proteins", f"{catalog.loc[catalog.protein != '', 'protein'].nunique():,}",
                          help="Each parent protein is counted once across variants, states and conditions.")
        total = pd.to_numeric(catalog.residue_count, errors="coerce").sum()
        metrics[2].metric("Residue measurements", f"{total:,.0f}")
        st.subheader("About the collection")
        st.markdown(
            r"Hydrogen-exchange protection factors and HX-derived opening free energies "
            r"($\Delta G_{\mathrm{op}}$), with experimental conditions, source references, sequences "
            "and structures where available."
        )
        st.write("The collection is almost entirely from conditions without added chemical denaturant. "
                 "A small number of denaturant-containing conditions and values extrapolated to zero "
                 "denaturant are also included.")
        st.subheader("PF standardization")
        st.markdown(
            "Standardized PFs use the updated intrinsic exchange-rate parameters of "
            "[Nguyen et al. (2018)](https://doi.org/10.1007/s13361-018-2021-z). "
            "Measured exchange rates were used when available; for PF-only studies, rates were "
            "back-calculated from published PFs using original intrinsic-rate values or documented "
            "historical models. PFs were then recalculated with the updated parameters. "
            "Published values and reconstruction assumptions are retained."
        )
        st.subheader("Explore the database")
        st.markdown(
            "**Browse** — Search or filter the collection, then select an entry for its PF profile, "
            "sequence, measurements, structure and experiment details.\n\n"
            "**Visualize** — Explore a single dataset as a PF profile, histogram, sequence map or 3D structure.\n\n"
            "**Compare** — Overlay PF profiles from related experiments and conditions."
        )
        show_recent_updates(updates)
    elif page == "Browse":
        show_browse(catalog)
    elif page == "Visualize":
        page_header("Visualize protection factors", "Explore residue profiles, PF distributions, and mapped protein structures.")
        eid = st.selectbox("Dataset", catalog.entry_id.tolist(), format_func=lambda value: label(value, catalog))
        show_entry(ROOT, eid, visual_only=True)
    else:
        show_compare(catalog)
    checked = datetime.now(ZoneInfo("America/New_York")).strftime("%H:%M:%S %Z")
    st.caption(f"Collection checked at {checked} · {'Automatic refresh on' if auto_refresh else 'Manual refresh'}")
    show_data_status(issues)


render_page()
