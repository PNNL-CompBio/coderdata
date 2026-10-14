#!/usr/bin/env python3
"""
cnf_utils.py: shared specimen classification utilities for the cNF build scripts.

Used by 01-samples-cnf.py, 02-omics-cnf.py, and 04-experiments-cnf.py to
canonicalize specimen strings so omics data, drug measurements, and the sample
table all use the same improve_sample_id keys.
"""

import re


# ----------------------------------------------------------------------------
# Sample type → coderdata model_type enum
# ----------------------------------------------------------------------------
MODEL_TYPE_MAP = {
    "organoid":         "patient derived organoid",
    "treated_organoid": "patient derived organoid",
    "tumor_tissue":     "tumor",
    "normal_skin":      "normal_tissue",
}




# ----------------------------------------------------------------------------
# Regex constants
# ----------------------------------------------------------------------------
# Bare tumor organoid: NFxxxx_Tx (4 patient digits, T + tumor digits)
TUMOR_RE          = re.compile(r"^NF\d{4}_T\d+$")
# Patient prefix at the start of any specimen string
PATIENT_PREFIX_RE = re.compile(r"^(NF\d{4})")
# Normal skin marker: NFxxxx[_-]skin..., allowing trailing date/replicate junk
SKIN_RE           = re.compile(r"^(NF\d{4})[_-]skin", re.IGNORECASE)
# Primary tumor tissue (not organoid culture): NFxxxx_Tx_tissue
TISSUE_RE         = re.compile(r"^(NF\d{4}_T\d+)_tissue$", re.IGNORECASE)
# Treated organoid: NFxxxx_Tx_<treatment_label>, e.g. NF0021_T1_Onalespid_1uM
TREATED_RE        = re.compile(r"^NF\d{4}_T\d+_.+$")


# ----------------------------------------------------------------------------
# Classification
# ----------------------------------------------------------------------------
def classify_specimen(spec) -> tuple[str, str] | None:
    """Classify a raw specimen string into (canonical_id, sample_type) or None.

    sample_type ∈ {"organoid", "treated_organoid", "tumor_tissue", "normal_skin"}.

    Canonical IDs:
        organoid         → NFxxxx_Tx          (suffix stripped)
        treated_organoid → NFxxxx_Tx_<treatment_label> (label preserved so
                                              treated samples don't collide
                                              with the untreated organoid)
        tumor_tissue     → NFxxxx_Tx_tissue   (suffix preserved)
        normal_skin      → NFxxxx_skin        (replicates collapsed per patient)

    Unrecognized prefixes and malformed inputs return None.
    """
    if spec is None:
        return None
    s = str(spec).strip()
    if not s or s.lower() == "nan":
        return None

    # Convert dot-separated form (e.g. RNA file's "NF0017.T1.organoid") to
    # underscore-separated.
    s_norm = s.replace(".", "_")

    # Skin: any string starting with NFxxxx that contains a "skin" marker
    skin_m = SKIN_RE.match(s_norm)
    if skin_m:
        return (f"{skin_m.group(1)}_skin", "normal_skin")

    # Tumor tissue: NFxxxx_Tx_tissue
    tissue_m = TISSUE_RE.match(s_norm)
    if tissue_m:
        return (f"{tissue_m.group(1)}_tissue", "tumor_tissue")

    # Strip organoid suffix (handles "_organoid" and "_organoids")
    s_stripped = re.sub(r"_organoids?$", "", s_norm, flags=re.IGNORECASE)

    # Bare organoid (untreated): NFxxxx_Tx exactly
    if TUMOR_RE.match(s_stripped):
        return (s_stripped, "organoid")

    # Treated organoid: NFxxxx_Tx_<anything else>
    if TREATED_RE.match(s_stripped):
        return (s_stripped, "treated_organoid")

    return None


def patient_from_specimen(specimen: str) -> str:
    """Extract NFxxxx patient ID from a canonical specimen ID."""
    m = PATIENT_PREFIX_RE.match(specimen)
    return m.group(1) if m else ""


def canonicalize_specimen_column(df, specimen_col: str):
    """Apply classify_specimen to a dataframe's specimen column in place.

    Returns the same df with:
      - rows whose specimen classifies dropped if None
      - specimen_col values replaced with their canonical IDs
      - a new column 'sample_type' added

    Used by build_omics.py and build_exp.py to map raw Synapse specimen
    strings (e.g. "NF0018.T1.organoid", "NF0021.T1.Onalespid.1uM") to
    canonical IDs that match the cnf_samples.csv `other_id` column.
    """
    classifications = df[specimen_col].apply(classify_specimen)
    keep = classifications.notna()
    df = df[keep].copy()
    classifications = classifications[keep]
    df[specimen_col] = [c[0] for c in classifications]
    df["sample_type"] = [c[1] for c in classifications]
    return df
