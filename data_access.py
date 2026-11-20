"""Read-only scientific data access, plus a separate app-local index activity cache."""

from __future__ import annotations

import json
import hashlib
import tempfile
import threading
import io
import math
import os
import re
import time
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import pandas as pd

from construct_labels import normalize_construct
from state_labels import overview_state
from display_labels import display_variant

ENTRY_PATTERN = re.compile(r"pfb\d+")
NUMERIC_COLUMNS = [
    "sequence_index", "kex_s_inv", "kint_reported_s_inv", "pf_reported",
    "kint_standard_s_inv", "pf_standard", "log10_pf_standard", "dG_open_kJ_mol",
]

# One absolute log10(PF) scale for every protein and every value-colored view.
PF_COLOR_MIN = 0.0
PF_COLOR_MAX = 10.0
UNMAPPED_COLOR = "#d3d7db"
PF_COLORSCALE = [
    (0.0, "#ffffff"),
    (0.2, "#deebf7"),
    (0.4, "#9ecae1"),
    (0.55, "#6baed6"),
    (0.7, "#3182bd"),
    (0.85, "#08519c"),
    (1.0, "#08306b"),
]


def default_data_root() -> Path:
    return Path(__file__).resolve().parent / "data"


def fingerprint(path: Path) -> tuple:
    stat = path.stat()
    return (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size, stat.st_ino)


def _current_bytes(filename: str, signature: tuple) -> bytes:
    """Cache only stable reads of the local files."""
    path = Path(filename)
    if fingerprint(path) != signature:
        raise OSError("File changed while being read")
    with path.open("rb", buffering=0) as handle:
        payload = handle.read()
    if fingerprint(path) != signature:
        raise OSError("File changed while being read")
    return payload


def _read_current(loader, path: Path):
    # A normal save can truncate or replace a file between stat() and read().
    # Retry against the new version; only successful stable reads enter caches.
    for attempt in range(3):
        try:
            return loader(str(path), fingerprint(path))
        except (OSError, ValueError, UnicodeError):
            if attempt == 2:
                raise
            time.sleep(0.02)


@lru_cache(maxsize=24)
def _json_cached(filename: str, signature: tuple) -> dict:
    value = json.loads(_current_bytes(filename, signature))
    if not isinstance(value, dict):
        raise ValueError("Entry metadata must be a JSON object")
    return value


def read_metadata(path: Path) -> dict:
    # Signatures, including ctime, refresh an edited file even when its size is unchanged.
    return _read_current(_json_cached, path)


@lru_cache(maxsize=48)
def _csv_cached(filename: str, signature: tuple) -> pd.DataFrame:
    # Preserve literal chain IDs such as 'NA', as well as blank chains/insertion codes.
    return pd.read_csv(io.BytesIO(_current_bytes(filename, signature)), dtype=str, keep_default_na=False)


def read_csv(path: Path) -> pd.DataFrame:
    return _read_current(_csv_cached, path).copy()


def read_structure(path: Path) -> str:
    def load(filename, signature):
        pdb_text = _current_bytes(filename, signature).decode("utf-8", errors="replace")
        if not pdb_atom_lookup(pdb_text):
            raise ValueError("Structure file has no readable atoms yet")
        return pdb_text
    return _read_current(load, path)


def text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    if isinstance(value, (list, tuple)):
        return ", ".join(text(item) for item in value if text(item))
    return str(value).strip()


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def pf_color(value) -> str:
    """Interpolate the shared white-to-blue scale, clipping only the display color."""
    numeric = number(value)
    if numeric is None:
        return UNMAPPED_COLOR
    fraction = max(0.0, min(1.0, (numeric - PF_COLOR_MIN) / (PF_COLOR_MAX - PF_COLOR_MIN)))
    for (low, start), (high, end) in zip(PF_COLORSCALE, PF_COLORSCALE[1:]):
        if fraction <= high:
            weight = (fraction - low) / (high - low)
            rgb = [round(int(start[i:i + 2], 16) * (1 - weight) + int(end[i:i + 2], 16) * weight)
                   for i in (1, 3, 5)]
            return "#" + "".join(f"{channel:02x}" for channel in rgb)
    return PF_COLORSCALE[-1][1]


def as_dict(value) -> dict:
    return value if isinstance(value, dict) else {}


def experimental_method_label(metadata: dict) -> str:
    """Separate the HX readout technique from kint models and PDB methods."""
    candidates = [metadata.get(key) for key in
                  ("experimental_method", "measurement_method", "technique")]
    candidates.append(as_dict(metadata.get("experiment")).get("method"))
    for candidate in candidates:
        if not text(candidate):
            continue
        token = re.sub(r"[^A-Z0-9]", "", text(candidate).upper())
        if token in {"NMR", "HXNMR", "HDXNMR"}:
            return "HX/NMR"
        if token in {"MS", "HXMS", "HDXMS", "MASSSPECTROMETRY"}:
            return "HX/MS"
        # An unfamiliar explicit method must not be overridden by a legacy default.
        return text(candidate)
    provenance = as_dict(metadata.get("dictionary_source"))
    if (Path(text(provenance.get("path"))).name == "nmr_residue_measurements.csv"
            and as_dict(metadata.get("generator")).get("name") == "standardize_nmr.py"):
        return "HX/NMR"
    return "Not recorded"


def modification_text(metadata: dict) -> str:
    modifications = metadata.get("sequence_modifications", metadata.get("modifications"))
    if not modifications:
        return ""
    if not isinstance(modifications, list):
        return text(modifications) if not isinstance(modifications, dict) else json.dumps(modifications)
    labels = []
    for item in modifications:
        if not isinstance(item, dict):
            labels.append(text(item))
            continue
        label = text(item.get("label") or item.get("description") or item.get("type")).replace("_", " ")
        if item.get("sequence"):
            label += ": " + text(item["sequence"])
        if item.get("construct_positions"):
            label += " (construct positions " + text(item["construct_positions"]) + ")"
        labels.append(label)
    return "; ".join(labels)


def structure_source_label(metadata: dict) -> str:
    """Keep predictor names separate from real experimental PDB identifiers."""
    structure = as_dict(metadata.get("structure"))
    baseline = as_dict(structure.get("baseline_reference"))
    details = as_dict(metadata.get("experimental_details"))
    source = text(metadata.get("structure_source") or structure.get("structure_source")
                  or details.get("structure_source") or baseline.get("structure_source"))
    token = re.sub(r"[^a-z0-9]", "", source.lower())
    if token in {"af3", "alphafold3"}:
        return "AF3"
    if token in {"af2", "alphafold2"}:
        return "AF2"
    origin = text(baseline.get("origin") or structure.get("origin")).lower()
    if "solution_nmr" in origin or token in {"nmr", "experimentalsolutionnmr"}:
        return "NMR"
    if token in {"experimentalpdb", "crystal", "xray", "experimentalxray"} or "xray" in origin.replace("_", ""):
        return "Crystal"
    if source:
        return source
    if origin == "predicted":
        return "Predicted"
    if origin in {"experimental", "locally_prepared"}:
        return "Experimental PDB"
    if structure.get("status") == "pending" and not baseline and not structure.get("coordinate_path"):
        return "No structure selected"
    return ""


def structure_preparation_label(metadata: dict) -> str:
    """Show concise preparation status without changing the recorded provenance."""
    if structure_source_label(metadata) in {"AF2", "AF3"}:
        return ""
    structure = as_dict(metadata.get("structure"))
    baseline = as_dict(structure.get("baseline_reference"))
    mutations = baseline.get("required_mutations") or structure.get("required_mutations")
    if baseline and mutations:
        labels = [text(item.get("label")) for item in mutations if isinstance(item, dict)]
        suffix = ": " + ", ".join(label for label in labels if label) if any(labels) else ""
        return "Mutation pending" + suffix
    method = text(structure.get("preparation_method")).lower()
    evidence = text(as_dict(structure.get("preparation_evidence")).get("note")).lower()
    description = method + " " + evidence
    if "rosetta" in description:
        if re.search(r"point[- ]mut|mutation|repacking", description):
            return "Crystal + Rosetta mutation" if structure_source_label(metadata) == "Crystal" else "Rosetta mutation"
        return "Crystal + Rosetta preparation" if structure_source_label(metadata) == "Crystal" else "Rosetta preparation"
    if structure.get("status") == "pending_preparation":
        return "Atom completion pending"
    if "component selection only" in method:
        return "Component selection"
    if method and structure.get("origin") == "locally_prepared":
        return "Crystal-derived preparation" if structure_source_label(metadata) == "Crystal" else "Local preparation"
    return ""


def summarize(metadata: dict) -> dict:
    conditions = as_dict(metadata.get("conditions"))
    structure = as_dict(metadata.get("structure"))
    labels = normalize_construct(metadata.get("construct", ""))
    kelvin = number(conditions.get("temperature_K"))
    pmid = text(metadata.get("PMID"))
    doi = text(metadata.get("DOI"))
    publication = text(metadata.get("publication_source")) or pmid or doi
    return {
        "protein": text(metadata.get("protein_id") or metadata.get("protein_name")),
        "protein_name": text(metadata.get("protein_name")),
        "pmid": pmid,
        "doi": doi,
        "publication": publication,
        "uniprot": text(metadata.get("uniprot_accession")),
        "construct": labels["construct"],
        "variant_label": display_variant(labels["construct"]),
        "construct_source_label": text(metadata.get("construct_source_label") or labels.get("construct_source_label")),
        "expression_tag": text(metadata.get("expression_tag") or labels.get("expression_tag")),
        "source_construct_id": text(metadata.get("source_construct_id") or labels.get("source_construct_id")),
        "modifications": modification_text(metadata),
        "state": text(metadata.get("state")),
        "state_label": overview_state(metadata),
        "temperature_c": kelvin - 273.15 if kelvin is not None else None,
        "ph_input": number(conditions.get("pH_input")),
        "ph_label": text(conditions.get("reported_pH_label")),
        "pdb_id": text(structure.get("pdb_id") or (";".join(ref["pdb_id"] for ref in structure.get("domain_references", [])) if structure.get("representation") == "separate_domain_references" else as_dict(structure.get("baseline_reference")).get("pdb_id"))),
        "structure_source": structure_source_label(metadata),
        "structure_preparation": structure_preparation_label(metadata),
        "residue_count": number(as_dict(metadata.get("checks")).get("quantitative_residues")),
    }


def catalog_display(catalog: pd.DataFrame) -> pd.DataFrame:
    """Display calculation inputs as numbers; source prose belongs in entry notes."""
    displayed = catalog[["entry_id", "protein", "construct", "state_label", "ph_input", "temperature_c",
                         "pmid", "uniprot", "pdb_id", "residue_count"]].copy()
    displayed["construct"] = catalog["construct"].map(display_variant)
    displayed["ph_input"] = pd.to_numeric(displayed["ph_input"], errors="coerce").round(2)
    displayed["temperature_c"] = pd.to_numeric(displayed["temperature_c"], errors="coerce").round(1)
    return displayed.rename(columns={
        "entry_id": "Entry", "protein": "Protein", "construct": "Variant", "state_label": "State",
        "pmid": "PMID", "uniprot": "UniProt", "temperature_c": "Temp",
        "ph_input": "pH", "pdb_id": "PDB", "residue_count": "Sites",
    })


@lru_cache(maxsize=2048)
def _summary_cached(filename: str, signature: tuple) -> dict:
    # Cache small summaries, not the large evidence sections in every metadata file.
    metadata = json.loads(_current_bytes(filename, signature))
    if not isinstance(metadata, dict):
        raise ValueError("Entry metadata must be a JSON object")
    return summarize(metadata)


def load_catalog(root: Path) -> tuple[pd.DataFrame, list[str]]:
    if not root.is_dir():
        raise FileNotFoundError(f"Training data folder not found: {root}")
    indexed = {}
    issues = []
    if (root / "index.csv").is_file():
        try:
            for row in read_csv(root / "index.csv").to_dict("records"):
                identifier = text(row.get("dataset_id"))
                if ENTRY_PATTERN.fullmatch(identifier):
                    indexed[identifier] = row
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            issues.append(f"Index is being updated or cannot be read: {exc}")
    # Discover complete new folders even before index.csv has been regenerated.
    with os.scandir(root) as entries:
        folders = sorted(
            (Path(entry.path) for entry in entries
             if ENTRY_PATTERN.fullmatch(entry.name) and entry.is_dir()), key=lambda p: p.name
        )
    rows = []
    for folder in folders:
        metadata_path = folder / "metadata.json"
        try:
            summary = _read_current(_summary_cached, metadata_path).copy()
        except (OSError, ValueError) as exc:
            issues.append(f"{folder.name}: metadata temporarily unavailable ({exc})")
            if folder.name not in indexed:
                continue
            old = indexed[folder.name]
            summary = summarize({
                "protein_id": old.get("protein_id"), "PMID": old.get("PMID"),
                "publication_source": old.get("publication_source"), "DOI": old.get("DOI"),
                "uniprot_accession": old.get("uniprot_accession"),
                "construct": old.get("construct"), "state": old.get("state"),
                "conditions": {"temperature_K": old.get("temperature_K"), "pH_input": old.get("pH_input")},
                "structure": {"pdb_id": old.get("pdb_id")},
            })
        if summary.get("residue_count") is None:
            summary["residue_count"] = number(indexed.get(folder.name, {}).get("quantitative_residues"))
        rows.append({"entry_id": folder.name, **summary})
    columns = ["entry_id", "protein", "protein_name", "pmid", "doi", "publication", "uniprot", "construct", "variant_label",
               "construct_source_label", "expression_tag", "source_construct_id",
               "modifications", "state", "state_label", "temperature_c", "ph_input", "ph_label", "pdb_id",
               "structure_source", "structure_preparation", "residue_count"]
    return pd.DataFrame(rows, columns=columns), issues


def entry_folder(root: Path, entry_id: str) -> Path:
    if not ENTRY_PATTERN.fullmatch(entry_id):
        raise ValueError("Invalid PFBoost entry ID")
    path = (root / entry_id).resolve()
    if path.parent != root.resolve():
        raise ValueError("Entry is outside the training data folder")
    return path


def load_entry(root: Path, entry_id: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    folder = entry_folder(root, entry_id)
    metadata = read_metadata(folder / "metadata.json")
    values = read_csv(folder / "residues.csv")
    if not {"sequence_index", "residue"}.issubset(values.columns):
        raise ValueError("residues.csv needs sequence_index and residue columns")
    for column in NUMERIC_COLUMNS:
        if column in values:
            values[column] = pd.to_numeric(values[column], errors="coerce")
            values[column] = values[column].replace([float("inf"), -float("inf")], float("nan"))
    values = values.sort_values("sequence_index")
    mapping_path = folder / "residue_mapping.csv"
    mapping = read_csv(mapping_path) if mapping_path.is_file() else pd.DataFrame()
    return metadata, values, mapping


def structure_files(folder: Path, metadata: dict) -> list[Path]:
    structure = as_dict(metadata.get("structure"))
    if structure.get("representation") == "separate_domain_references":
        refs = structure.get("domain_references", [])
        names = [ref.get("coordinate_file", "") for ref in refs]
        if not names or len(names) != len(set(names)) or any(Path(name).name != name or not name.endswith(".pdb") for name in names):
            raise ValueError("Separate domain references require unique entry-relative PDB files")
        paths = [folder / name for name in names]
        if not all(path.is_file() for path in paths):
            raise ValueError("A declared separate domain coordinate file is missing")
        return paths
    candidates = {path.name: path for path in folder.glob("structure*.pdb") if path.is_file()}
    structure = as_dict(metadata.get("structure"))
    named = text(structure.get("coordinate_file") or structure.get("coordinate_path"))
    basename = Path(named).name if named else ""
    if basename.endswith(".pdb") and (folder / basename).is_file():
        candidates[basename] = folder / basename
    ordered = []
    for name in ["structure.pdb", basename]:
        if name in candidates:
            ordered.append(candidates.pop(name))
    return ordered + [candidates[name] for name in sorted(candidates)]


def select_pdb_model(pdb_text: str, preferred_model=None) -> tuple[str, int]:
    """Render one model, preserving atom serials and original chain/residue labels."""
    models = {}
    current = None
    for line in pdb_text.splitlines():
        if line.startswith("MODEL "):
            try:
                current = int(line[10:14])
            except ValueError:
                current = len(models) + 1
            models.setdefault(current, [])
        elif line.startswith("ENDMDL"):
            current = None
        elif current is not None:
            models[current].append(line)
    if not models:
        return pdb_text, 1
    preferred = number(preferred_model)
    selected = int(preferred) if preferred in models else next(iter(models))
    return "\n".join(models[selected]) + "\nEND\n", selected


def pdb_atom_lookup(pdb_text: str) -> dict:
    """Resolve exact PDB keys to atom serials, including blank chains and insertion codes."""
    atoms = defaultdict(set)
    serial_keys = defaultdict(set)
    for line in pdb_text.splitlines():
        if line[:6] not in {"ATOM  ", "HETATM"}:
            continue
        try:
            key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip())
            serial = int(line[6:11])
        except (ValueError, IndexError):
            continue
        atoms[key].add(serial)
        serial_keys[serial].add(key)
    # Do not color an unrelated residue in malformed files with repeated atom serials.
    ambiguous = {serial for serial, keys in serial_keys.items() if len(keys) > 1}
    return {key: serials - ambiguous for key, serials in atoms.items()}


def mapped_colors(values: pd.DataFrame, mapping: pd.DataFrame, pdb_text: str,
                  selected_model: int, *, coordinate_file: str | None = None) -> tuple[list[dict], int]:
    if "coordinate_file" in mapping.columns:
        if not coordinate_file:
            raise ValueError("Select a coordinate file before coloring separate domain references")
        mapping = mapping[mapping["coordinate_file"].fillna("") == coordinate_file]
    required = {"sequence_index", "pdb_chain", "pdb_residue_number", "pdb_insertion_code"}
    if "log10_pf_standard" not in values or not required.issubset(mapping.columns):
        return [], 0
    lookup = pdb_atom_lookup(pdb_text)
    observations = {}
    for row in values.to_dict("records"):
        index, pf = number(row.get("sequence_index")), number(row.get("log10_pf_standard"))
        if index is not None and pf is not None:
            observations.setdefault(int(index), []).append((pf, text(row.get("residue"))))
    targets = defaultdict(list)
    for row in mapping.to_dict("records"):
        index, residue_number = number(row.get("sequence_index")), number(row.get("pdb_residue_number"))
        if index is None or residue_number is None or int(index) != index or int(residue_number) != residue_number:
            continue
        model = number(row.get("pdb_model"))
        if model is not None and model != selected_model:
            continue
        key = (text(row.get("pdb_chain")), int(residue_number), text(row.get("pdb_insertion_code")))
        if not lookup.get(key):
            continue
        for pf, residue in observations.get(int(index), []):
            mapped_residue = text(row.get("residue"))
            if mapped_residue and residue and mapped_residue != residue:
                continue
            targets[key].append((int(index), pf))
    styles = []
    mapped_indices = set()
    for key, candidates in targets.items():
        factors = {pf for _, pf in candidates}
        if len(factors) != 1:
            continue
        pf = next(iter(factors))
        styles.append({"selection": {"model": 0, "serial": sorted(lookup[key])},
                       "color": pf_color(pf), "coordinate_key": key})
        mapped_indices.update(index for index, _ in candidates)
    return styles, len(mapped_indices)


def clear_caches():
    _json_cached.cache_clear()
    _csv_cached.cache_clear()
    _summary_cached.cache_clear()


# A compact Home activity strip, driven only by the main index. Its one cache
# file lives in the app, never in the scientific collection or either index.
_RECENT_UPDATE_LOCK = threading.RLock()
_RECENT_UPDATE_MEMORY = {}
_STRUCTURE_INDEX_FIELDS = {"pdb_id", "structure_origin", "structure_file", "structure_path",
                           "structure_name", "pdb_path", "coordinate_file", "coordinate_path"}
_PF_INDEX_FIELDS = {"values_csv", "quantitative_residues", "source_rate_dataset", "source_pf_dataset",
                    "pf_sha256", "measurements_sha256", "residues_sha256"}


def recent_index_updates(root: Path, *, state_path: Path | None = None) -> list[dict]:
    source = Path(root).resolve() / "index.csv"
    if state_path is None:
        identifier = hashlib.sha256(str(source).encode()).hexdigest()[:12]
        state_path = Path(__file__).resolve().parent / ".runtime" / f"recent-updates-{identifier}.json"
    state_path = Path(state_path)
    with _RECENT_UPDATE_LOCK:
        try:
            state = read_metadata(state_path)
        except (OSError, ValueError, UnicodeError):
            state = _RECENT_UPDATE_MEMORY.get(str(state_path), {})
        if state.get("source") != str(source):
            state = {}
        # Ignore the old file-refresh placeholder, including on the unchanged
        # signature fast path. Only changes to actual entry rows are activity.
        updates = [event for event in state.get("updates", [])
                   if event.get("kind") in {"added", "structure", "pf", "entry"}
                   and event.get("count", 0) > 0][:2]
        try:
            signature = fingerprint(source)
            if list(signature) == state.get("signature"):
                return updates
            frame = read_csv(source)
            if "dataset_id" not in frame:
                return updates
            rows = {text(row["dataset_id"]): {key: text(value) for key, value in row.items()}
                    for row in frame.to_dict("records") if ENTRY_PATTERN.fullmatch(text(row["dataset_id"]))}
            # Ignore empty/partial index saves; existing entries remain usable.
            if not rows or fingerprint(source) != signature:
                return updates
        except (OSError, ValueError, UnicodeError, pd.errors.ParserError):
            return updates
        changed_at = signature[0] / 1_000_000_000
        previous = state.get("rows", {})
        if not previous:
            # Establish the first snapshot silently; it is not an update.
            updates = []
        else:
            groups = defaultdict(list)
            for entry_id, row in rows.items():
                old = previous.get(entry_id)
                if old is None:
                    kind = "added"
                else:
                    changed = {key for key in row.keys() | old.keys() if row.get(key, "") != old.get(key, "")}
                    if not changed:
                        continue
                    kind = ("structure" if changed & _STRUCTURE_INDEX_FIELDS else
                            "pf" if changed & _PF_INDEX_FIELDS else "entry")
                groups[kind].append({"id": entry_id, "protein": row.get("protein_id", "")})
            new = [{"kind": kind, "count": len(groups[kind]), "entries": groups[kind][:2], "at": changed_at}
                   for kind in ("added", "structure", "pf", "entry") if groups[kind]]
            updates = (new + updates)[:2]
        # Keep known rows absent from a partial save, so restoring that file
        # cannot generate a burst of false 'new dataset' announcements.
        state = {"source": str(source), "signature": list(signature),
                 "rows": {**previous, **rows}, "updates": updates}
        _RECENT_UPDATE_MEMORY[str(state_path)] = state
        temporary = None
        try:
            state_path.parent.mkdir(parents=True, exist_ok=True)
            descriptor, temporary = tempfile.mkstemp(prefix=".recent-updates-", dir=state_path.parent)
            with os.fdopen(descriptor, "w") as handle:
                json.dump(state, handle, separators=(",", ":"))
                handle.write("\n")
            os.chmod(temporary, 0o644)
            os.replace(temporary, state_path)
        except OSError:
            pass  # A read-only app folder must not prevent viewing the data.
        finally:
            if temporary and os.path.exists(temporary):
                os.unlink(temporary)
        return updates
