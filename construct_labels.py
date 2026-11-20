"""Reviewed, conservative construct-label normalization for PFBoost.

Only label fields are changed. Original descriptions remain available as
construct_source_label; source identifiers and expression tags have their own fields.
Unknown constructs are never assumed to be wild type.
"""

import re

POLICY = "pfboost_construct_labels_v1"

LABEL_REPLACEMENTS = {
    "Small archaeal modifier protein 1 from Methanosarcina acetivorans": "SAMP1 + N-terminal GH",
    "galectin-8_NTD": "Galectin-8 NTD",
    "galectin-8_CTD": "Galectin-8 CTD",
    "galectin-1": "Galectin-1",
    "galectin-7": "Galectin-7",
    "PSE-4_mature_enzyme": "Mature PSE-4",
    "recombinant uniformly15N-enriched": "Reference construct",
    "15N labeled": "Reference construct",
    "human DHFR L80F": "L80F",
    "human DHFR WT": "WT",
    "human DHFR D153V": "D153V",
    "C5V/C55A; [30-51;14-38]V5A55": "C5V/C55A",
    "C14A/C30A/C38A/C51A; [5-55]Ala": "C14A/C30A/C38A/C51A",
    "M7W/K98I/N99R/H102N/R106G (Rd-apocyt b562)": "M7W/K98I/N99R/H102N/R106G",
    "wild type; two native disulfides": "WT",
    "RNase H* C13A/C63A/C133A": "C13A/C63A/C133A (RNase H*)",
    "variant A; recombinant labeled protein": "Variant A",
    "wild type; 15N-labeled": "WT",
    "WT (not F41W)": "WT",
    "KI-FHA domain; native KAPP numbering": "KI-FHA domain",
    "C69A pseudo-wild type": "C69A (pseudo-WT)",
    "wild-type adenylate kinase": "WT",
    "C54T/C97A (WT*)": "C54T/C97A (WT*)",
    "WT* C54T/C97A": "C54T/C97A (WT*)",
    "N-terminalMet0;native99aa": "Mature 1–99 + N-terminal Met",
    "VH-(Gly4Ser)3-VL; domain H": "scFv VH domain",
    "VH-(Gly4Ser)3-VL; domain L": "scFv VL domain",
    "predominantly apo (~80%); heterogeneous N-terminal methionine/prosthetic group": "Apo-dominant ACP",
    "hFGF-1 residues 15-154; numbering of full 154-residue protein": "FGF1 (15–154)",
    "recombinant human alpha-lactalbumin; mature numbering1-123; no mutation reported": "Mature 1–123",
    "SN-OB residues1-103; V66L/G88V": "V66L/G88V (1–103)",
    "C31S; N-terminal His6 tag MHHHHHHLVPRGM": "C31S",
    "C102T; [U-15N] labeled": "C102T",
    "C102T; unlabeled": "C102T",
    "recombinant wild-type GMH4; 147 residues; [U-15N]; reconstituted with b-type heme": "WT",
    "F2 (P88K/S92K); 15N-labeled": "F2 (P88K/S92K)",
    "F4 (P88K/A90L/S92K/A94L); 15N-labeled": "F4 (P88K/A90L/S92K/A94L)",
    "AhpF residues 1-202; wild-type isolated NTD": "WT NTD (1–202)",
    "wild type; Sigma type VI": "WT",
    "T70N; mature 130-residue human lysozyme": "T70N (mature 130 aa)",
    "GCN4-p2N; N-acetylated, C-amidated, Cys(-1)-linked disulfide homodimer": "GCN4-p2N",
    "Wild-type mature HdeA residues1-89": "WT (mature 1–89)",
    "untagged; no internal sequence mutation reported": "Untagged reference construct",
    "wild type; Ala54; mature 131-residue chain": "WT (A54; mature 1–131)",
    "IKBa(67-206),free,2H/15N-labeled": "IκBα (67–206)",
    "21aaNterminalHislinker plus43aa design": "HEEH_TK_rd5_0341",
    "full-length residues79-146 with N-terminal Gly-Ala remnant; numbering as Feng1994": "SH3 (79–146) + N-terminal GA",
    "histidine-tagged wild-type Im7 (Im7*)": "WT (Im7*)",
    "wild type;GST-cleaved;223aa verified human core;exact scar unreported": "WT (223 aa)",
    "recombinant MpDHFR;162aa MIV...K experimental sequence;Met1 verified in TableS24 and BMRB17946": "MpDHFR (162 aa)",
    "recombinant EcDHFR;159aa reference-verified core;terminal extensions not independently resolved": "EcDHFR (159 aa core)",
    "C69A; uniformly 15N-labelled; source/PDB verified 179-residue core": "C69A (179 aa core)",
    "V173M/T174M; mature core81-224; N- and C-terminal His6 tags reported": "V173M/T174M (81–224)",
    "commercial15N/13C ubiquitin; no mutation/tag specified": "Ubiquitin reference construct",
    "T2Q (source WT GB1)": "T2Q",
    "T2Q (source GB1)": "T2Q",
    "truncated CI2; full-length20-83; source numbering1-64": "CI2 (20–83)",
    "truncated CI2 I20V (full-length I39V)": "I20V (CI2 20–83; full-length I39V)",
    "WT core1-258; extra N-terminal GSH after affinity-tag cleavage": "WT (1–258) + N-terminal GSH",
    "H-Ras1-166 WT": "WT (1–166)",
    "H-Ras1-166 Y157Q": "Y157Q (1–166)",
    "H-Ras1-166 L120A": "L120A (1–166)",
    "H-Ras1-166 Q99A": "Q99A (1–166)",
    "H-Ras1-166 H27G": "H27G (1–166)",
    "H-Ras1-166 T35S": "T35S (1–166)",
    "recombinant rhIL-1ra;verified152aa maturecore;paper numbering=core+1;extra N-terminal identity not specified": "Mature IL-1ra (152 aa core)",
    "human FABP1 Ala94 sequence variant; MRGS N-terminal and LVPR C-terminal extensions; 135 residues including 127-residue core": "A94 variant (MRGS–FABP1–LVPR)",
    "monomeric bovine seminal RNase N67D;124-residue core;N-terminal Met0 and Cys31/Cys32 glutathione attachments from cited construct provenance": "N67D (monomeric BS-RNase)",
    "wild-type E.coli CspA;paper-cited core author2-70;N-terminal processing not independently verified": "WT CspA (2–70 core)",
    "wild-type97-residueGroES;native free homoheptamer": "WT (97 aa)",
    "WT full-length hAQP1(1-269)-linker-Myc-His6;292aa;U13C/U15N": "WT AQP1 (1–269)",
    "WT recombinant mature human gammaS;His-SUMO cleaved;177aa reference core": "WT (mature 177 aa)",
    "N14D recombinant mature human gammaS;His-SUMO cleaved;177aa reference core": "N14D (mature 177 aa)",
    "N76D recombinant mature human gammaS;His-SUMO cleaved;177aa reference core": "N76D (mature 177 aa)",
    "WT;human T141-E335 core195aa;HX tag retention unreported": "WT (141–335)",
    "E317A;human T141-E335 core195aa;HX tag retention unreported": "E317A (141–335)",
    "WT;GSH+humanDJ1core1-189": "WT (1–189) + N-terminal GSH",
    "M26I;GSH+humanDJ1core1-189": "M26I (1–189) + N-terminal GSH",
    "M51V;initiator Met1 removed;author2-88": "M51V (2–88)",
    "WT/CLIP102-120M107W;alpha1-192/beta1-198": "WT HLA-DR1 / CLIP M107W (102–120)",
    "betaN82A/CLIP102-120M107W;alpha1-192/beta1-198": "βN82A HLA-DR1 / CLIP M107W (102–120)",
    "WT mature 153-aa core": "WT (mature 153 aa)",
    "PM23 circular permutant; TAQT linker; WT residue numbering; exact full sequence unassigned": "PM23 circular permutant",
    "PM65 circular permutant; TAQT linker; WT residue numbering; exact full sequence unassigned": "PM65 circular permutant",
    "PM76 circular permutant; TAQT linker; WT residue numbering; exact full sequence unassigned": "PM76 circular permutant",
    "PM142 circular permutant; TAQT linker; WT residue numbering; exact full sequence unassigned": "PM142 circular permutant",
    "SCI-b;B28D;B29P;A8H;A14E;six-residue EEGPRR linker": "SCI-b (B28D/B29P/A8H/A14E)",
    "lambda-YA (nominal residues 6-85); Y22W,Q33Y,G46A,G48A in paper numbering": "λ-YA (Y22W/Q33Y/G46A/G48A; 6–85)",
    "WT human PRC1 dimerization domain;native1-65;C-terminal His-tag described,exact tail unresolved": "WT dimerization domain (1–65)",
    "WT Thermus thermophilus ribosomal protein S6;101aa": "WT (101 aa)",
    "Thermus thermophilus S6 circular permutant P54-55;96aa;newM1 and ASTTPG linker": "P54–55 circular permutant (96 aa)",
}

EXPLICIT_TAGS = {
    "C31S; N-terminal His6 tag MHHHHHHLVPRGM": "His6",
    "untagged; no internal sequence mutation reported": "None",
    "21aaNterminalHislinker plus43aa design": "His-tag/linker (21 aa)",
    "histidine-tagged wild-type Im7 (Im7*)": "His-tag (length unspecified)",
    "V173M/T174M; mature core81-224; N- and C-terminal His6 tags reported": "N- and C-terminal His6",
    "WT full-length hAQP1(1-269)-linker-Myc-His6;292aa;U13C/U15N": "Myc–His6",
    "WT recombinant mature human gammaS;His-SUMO cleaved;177aa reference core": "His–SUMO removed",
    "N14D recombinant mature human gammaS;His-SUMO cleaved;177aa reference core": "His–SUMO removed",
    "N76D recombinant mature human gammaS;His-SUMO cleaved;177aa reference core": "His–SUMO removed",
    "WT human PRC1 dimerization domain;native1-65;C-terminal His-tag described,exact tail unresolved": "C-terminal His-tag (tail unspecified)",
}


def normalize_construct(label: str) -> dict:
    original = str(label or "").strip()
    result = {"construct": original}
    match = re.fullmatch(r"(No-His|6xHis);\s*HX_ID=(\S+)", original, re.IGNORECASE)
    if match:
        source_id = match[2]
        mutation_match = re.fullmatch(r"(.+?)((?:_[A-Z]\d+[A-Z])+)", source_id)
        result["construct"] = (
            f"{mutation_match[1]} ({'/'.join(mutation_match[2].strip('_').split('_'))})"
            if mutation_match else source_id
        )
        result["expression_tag"] = "None" if match[1].lower() == "no-his" else "His6"
        result["source_construct_id"] = source_id
    else:
        result["construct"] = LABEL_REPLACEMENTS.get(original, re.sub(r"\bwild[ -]?type\b", "WT", original, flags=re.IGNORECASE))
        if original in EXPLICIT_TAGS:
            result["expression_tag"] = EXPLICIT_TAGS[original]
    if result["construct"] != original:
        result["construct_source_label"] = original
        result["construct_label_policy"] = POLICY
    return result
