"""Readable, conservative experimental-state labels for the data browser."""

STATE_LABELS = {
    "APO": "Apo",
    "LACTOSE": "Lactose bound",
    "APO_DOMINANT_MIXTURE": "Apo-dominant mixture",
    "NATIVE": "Native",
    "NATIVE_LIKE": "Native-like",
    "MOLTEN_GLOBULE": "Molten globule",
    "UNFOLDED": "Unfolded",
    "ACID_ACTIVATED": "Acid-activated",
    "STATE_1": "State 1",
    "STATE_2": "State 2",
    "STATE_1_2_ENSEMBLE": "State 1 / state 2 ensemble",
    "STATE_1_PREDOMINANT": "Predominantly state 1",
    "PREDOMINANTLY_FOLDED": "Predominantly folded",
    "FOLDED": "Folded",
    "ACID_A_STATE": "Acid A state",
    "ACID_INTERMEDIATE_I": "Acid intermediate I",
    "ACID_ACTIVATED_PARTIALLY_UNFOLDED": "Acid-activated, partially unfolded",
    "ACIDIC_IIB_IIC_UNFOLDED_ENSEMBLE": "Acidic IIb / IIc / unfolded ensemble",
    "ALKALINE_HEME_LIGATION_ENSEMBLE": "Alkaline heme-ligation ensemble",
    "KINETIC_BURST_INTERMEDIATE": "Kinetic burst-phase intermediate",
    "THERMALLY_UNFOLDING": "Thermally unfolding during HX",
    "CLOSED": "Closed",
    "TRANSITIONAL": "Transitional",
    "REDUCED": "Reduced",
    "OXIDIZED": "Oxidized",
    "MONOMER": "Monomer",
    "DIMER": "Dimer",
    "TETRAMER": "Tetramer",
    "TRIMER": "Trimer",
    "HEPTAMER": "Heptamer",
    "HETERODIMER": "Heterodimer",
    "PREDOMINANTLY_DIMER": "Predominantly dimeric",
    "MONOMER_DIMER_EQUILIBRIUM": "Monomer / dimer equilibrium",
    "INITIALLY_TETRAMERIC": "Initially tetrameric",
    "IRON": "Iron",
    "HEME": "Heme",
    "HEME_CO": "Heme / carbon monoxide",
    "CALCIUM": "Calcium",
    "CALCIUM_PDTP": "Calcium / pdTp",
    "MAGNESIUM": "Magnesium",
    "FOLATE": "Folate",
    "NADP_FOLATE": "NADP⁺ / folate",
    "BETA_CYCLODEXTRIN": "β-Cyclodextrin",
    "PALMITOYL_COA": "Palmitoyl-CoA",
    "PHOSPHORYLCHOLINE": "Phosphorylcholine",
    "AP5A": "Ap5A",
    "FMN": "FMN",
    "GPPNHP": "GppNHp",
    "OLEATE": "Oleate",
    "GLYCOCHOLATE": "Glycocholate",
    "CLIP": "CLIP peptide",
    "CP_THIOLATE": "Catalytic-cysteine thiolate",
}

OVERVIEW_LABELS = {
    **STATE_LABELS, "APO": "APO", "APO_DOMINANT_MIXTURE": "Apo-dominant mix",
    "LACTOSE": "Lactose", "HEME_CO": "Heme / CO",
    "STATE_1_2_ENSEMBLE": "States 1/2", "STATE_1_PREDOMINANT": "Mostly state 1",
    "PREDOMINANTLY_FOLDED": "Mostly folded", "PREDOMINANTLY_DIMER": "Mostly dimeric",
    "ACID_ACTIVATED_PARTIALLY_UNFOLDED": "Acid-activated / partially unfolded",
    "ACIDIC_IIB_IIC_UNFOLDED_ENSEMBLE": "IIb/IIc/unfolded",
    "ALKALINE_HEME_LIGATION_ENSEMBLE": "Alkaline ligation ensemble",
    "KINETIC_BURST_INTERMEDIATE": "Burst-phase intermediate",
    "MONOMER_DIMER_EQUILIBRIUM": "Monomer/dimer equilibrium",
    "THERMALLY_UNFOLDING": "Thermally unfolding",
}

ASSEMBLY_STATES = {
    "MONOMER", "DIMER", "TETRAMER", "TRIMER", "HEPTAMER", "HETERODIMER",
    "PREDOMINANTLY_DIMER", "MONOMER_DIMER_EQUILIBRIUM", "INITIALLY_TETRAMERIC",
}
ROUTINE_CONFORMATIONS = {"NATIVE", "FOLDED"}
BINDING_STATES = {
    "APO", "LACTOSE", "APO_DOMINANT_MIXTURE", "IRON", "HEME", "HEME_CO",
    "CALCIUM", "CALCIUM_PDTP", "MAGNESIUM", "FOLATE", "NADP_FOLATE",
    "BETA_CYCLODEXTRIN", "PALMITOYL_COA", "PHOSPHORYLCHOLINE", "AP5A",
    "FMN", "GPPNHP", "OLEATE", "GLYCOCHOLATE", "CLIP",
}
REDOX_STATES = {"REDUCED", "OXIDIZED"}
CONFORMATION_STATES = set(STATE_LABELS) - BINDING_STATES - ASSEMBLY_STATES - REDOX_STATES - {"CP_THIOLATE"}


def display_binding(metadata):
    """Return (label, assumed); the owner's default never rewrites source metadata."""
    binding = metadata.get("binding_state")
    if isinstance(binding, str) and binding.strip():
        return binding.strip(), False
    legacy = metadata.get("state")
    legacy = legacy.strip() if isinstance(legacy, str) else ""
    token = legacy.upper()
    if token in BINDING_STATES:
        return token, False
    # Preserve short explicit partner tokens in older records, such as DNA/HVEM.
    if legacy.isupper() and len(legacy) <= 32 and legacy.replace("_", "").isalnum() and token not in STATE_LABELS:
        return legacy, False
    return "APO", True


def overview_state(metadata):
    """Keep binding, non-native conformations and redox differences visible."""
    def value(key):
        item = metadata.get(key)
        return item.strip() if isinstance(item, str) else ""

    binding, _ = display_binding(metadata)
    conformation, redox, assembly = [value(key) for key in (
        "conformation", "redox_state", "assembly")]
    legacy = value("state")
    if not conformation and legacy.upper() in CONFORMATION_STATES:
        conformation = legacy.upper()
    if not redox and legacy.upper() in REDOX_STATES:
        redox = legacy.upper()
    parts = [OVERVIEW_LABELS.get(binding, binding)] if binding else []
    partner_variant = value("binding_partner_variant")
    if binding and partner_variant:
        parts[0] += f" ({partner_variant})"
    if conformation and conformation not in ROUTINE_CONFORMATIONS:
        parts.append(OVERVIEW_LABELS.get(conformation, conformation))
    if redox:
        parts.append(OVERVIEW_LABELS.get(redox, redox))
    # Preserve an explicit chemical state without treating it as a binding claim.
    if legacy == "CP_THIOLATE":
        parts.append(OVERVIEW_LABELS[legacy])
    return " · ".join(parts)

