#!/usr/bin/env python3
"""
04-experiments-cnf.py: generate cnf_experiments.tsv.

Discovers every drug-screen viability file under three Synapse parent folders,
attaches a specimen ID by reading the file's `specimenID` annotation (with a
filename-pattern fallback), then per (specimen, drug):

  * Multi-concentration measurements are run through the standard
    coderbuild/utils/fit_curve.py utility, producing fit_auc / fit_ic50 /
    fit_einf / fit_hs / fit_r2.
  * Single-concentration measurements at 1 μM are recorded with
    metric = 'uM_viability' and value = Viability_percentage / 100.
    Schema dependency: 'uM_viability' must be in the ResponseMetric enum.
    See schema_patch.md.

Specimen attribution policy:
  1. Trust the file's `specimenID` annotation if present.
  2. Fall back to parsing the filename (e.g. "NF0017_T2_Viabilities.csv"
     → "NF0017_T2").
  3. Run the result through cnf_utils.classify_specimen as a safety net
     (handles dot-separators, suffixes, case variants). Files whose
     specimen can't be classified are skipped with a warning.

Dual-mapping: one viability file's results should be attributed to two
specimens. Handled via SPECIMEN_DUAL_MAPPINGS at the top of the script;
add new entries as more cases come up.

Usage:
    python 04-experiments-cnf.py <cnf_samples.csv> <cnf_drugs.tsv>
                        [--output /tmp/cnf_experiments.tsv]
                        [--parents syn51301414,syn51301420,syn51301426]
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd
import synapseclient

from cnf_utils import classify_specimen


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
# Synapse folders to walk for drug-screen viability files
DEFAULT_DRUG_SCREEN_PARENTS = [
    "syn51301414",
    "syn51301420",
    "syn51301426",
]

TIME              = 120
TIME_UNIT         = "hours"
STUDY             = "cnf"
SOURCE            = "Synapse"
SINGLE_DOSE_UM    = 1.0
SINGLE_DOSE_METRIC = "uM_viability"

FIT_CURVE_SCRIPT  = os.environ.get(
    "CODERDATA_FIT_CURVE_SCRIPT", "/coderbuild/utils/fit_curve.py")

DEFAULT_OUTPUT = "/tmp/cnf_experiments.tsv"

# Regex to recover specimen from filenames like "NF0017_T2_Viabilities.csv"
# or "NF0021_T1_Onalespid_1uM_Viabilities.csv"
FILENAME_SPECIMEN_RE = re.compile(
    r"^(NF\d{4}(?:_T\d+)?(?:_[A-Za-z0-9._-]+?)?)"   # specimen prefix
    r"_(?:Viabilities|viabilities|Viability|drug_screen)"  # known label
    r"\.(?:csv|tsv|txt)$",
    re.IGNORECASE,
)

# ----------------------------------------------------------------------------
# DUAL-MAPPING TABLE
# ----------------------------------------------------------------------------
# Maps a canonical specimen to any additional specimens that should receive
# the same drug-response rows. Add new entries as more cases come up.
#
# Known case: NF0021_T1_Onalespid_1uM was screened against the same drug
# panel as NF0021_T1 and should be available under both labels.
SPECIMEN_DUAL_MAPPINGS: dict[str, list[str]] = {
    "NF0021_T1_Onalespid_1uM": ["NF0021_T1"],
}


def configure_logging():
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


# ----------------------------------------------------------------------------
# Specimen attribution
# ----------------------------------------------------------------------------
def specimen_from_annotation(syn, file_id: str) -> str | None:
    """Return the file's specimenID annotation, or None if missing."""
    try:
        anno = syn.getAnnotations(file_id)
    except Exception as exc:
        logging.warning("getAnnotations failed for %s: %s", file_id, exc)
        return None
    val = anno.get("specimenID")
    if isinstance(val, list):
        val = val[0] if val else None
    if val is None:
        return None
    s = str(val).strip()
    return s or None


def specimen_from_filename(name: str) -> str | None:
    """Parse 'NF0017_T2_Viabilities.csv' → 'NF0017_T2'."""
    m = FILENAME_SPECIMEN_RE.match(name)
    if m:
        return m.group(1)
    return None


def resolve_specimen(syn, file_id: str, file_name: str) -> str | None:
    """Resolve a file to a canonical specimen ID, trying annotation first.

    Returns canonical_id (e.g. "NF0017_T2") or None if both annotation and
    filename parsing fail or the result doesn't classify as a known sample
    type.
    """
    raw = specimen_from_annotation(syn, file_id)
    if not raw:
        raw = specimen_from_filename(file_name)
        if not raw:
            return None

    # Safety-net classification — handles dot-separators, suffixes, case
    result = classify_specimen(raw)
    if result is None:
        # Annotation or filename produced something that doesn't fit any
        # known sample pattern.
        logging.warning(
            "Could not classify specimen '%s' (file %s, %s); skipping",
            raw, file_id, file_name,
        )
        return None
    return result[0]


def expand_dual_mappings(specimen: str) -> list[str]:
    """Return [specimen] + any additional specimens it's dual-mapped to."""
    return [specimen] + SPECIMEN_DUAL_MAPPINGS.get(specimen, [])


# ----------------------------------------------------------------------------
# Synapse walking
# ----------------------------------------------------------------------------
def walk_for_viability_files(syn, parent_id: str) -> list[dict]:
    """Recursively yield File entities under a Synapse folder.

    Returns a list of {id, name} dicts for files whose name suggests they
    contain viability data (drug-screen CSV/TSV).
    """
    try:
        import synapseutils
    except ImportError:
        logging.error("synapseutils not available; cannot walk %s", parent_id)
        return []

    files = []
    for dirpath, _subfolders, file_list in synapseutils.walk(syn, parent_id):
        for fname, fid in file_list:
            # Filter to plausible viability files. Accept anything ending in
            # csv/tsv that has 'viab' or 'drug' in the name; fall back to
            # 'csv' if neither matches but the file is small enough.
            lower = fname.lower()
            if not lower.endswith((".csv", ".tsv", ".txt")):
                continue
            if any(tok in lower for tok in ("viab", "drug", "screen")):
                files.append({"id": fid, "name": fname,
                              "parent_path": dirpath})
            # Don't auto-grab arbitrary CSVs — we only want drug screens.
    return files


def discover_drug_files(syn, parents: list[str]) -> pd.DataFrame:
    """Walk all parent folders and dedupe by Synapse ID."""
    all_files = []
    for parent in parents:
        logging.info("Walking %s ...", parent)
        files = walk_for_viability_files(syn, parent)
        logging.info("  found %d candidate viability files", len(files))
        all_files.extend(files)

    if not all_files:
        return pd.DataFrame(columns=["id", "name"])

    df = pd.DataFrame(all_files).drop_duplicates(subset="id").reset_index(drop=True)
    logging.info("Discovered %d unique viability files across %d parents",
                 len(df), len(parents))
    return df


# ----------------------------------------------------------------------------
# Read and combine drug-screen files
# ----------------------------------------------------------------------------
def pull_screen_files(syn, parents: list[str]) -> pd.DataFrame:
    """Walk parents, read each viability CSV, attach canonical specimen.

    Files whose specimen can't be resolved are skipped with a warning.
    Files attributed to a dual-mapped specimen produce duplicate rows for
    every additional specimen in the mapping.
    """
    files_df = discover_drug_files(syn, parents)
    if files_df.empty:
        raise RuntimeError("No drug-screen files found under any parent")

    frames = []
    for _, row in files_df.iterrows():
        fid, fname = row["id"], row["name"]

        canonical = resolve_specimen(syn, fid, fname)
        if not canonical:
            continue

        try:
            df = pd.read_csv(syn.get(fid).path)
        except Exception as exc:
            logging.warning("Skipping file %s: %s", fid, exc)
            continue

        # Verify required columns
        required = {"Drug", "Concentration_uM", "Viability_percentage"}
        missing = required - set(df.columns)
        if missing:
            logging.warning("File %s missing %s; skipping",
                            fname, sorted(missing))
            continue

        # Attribute to canonical + any dual-mapped specimens
        for specimen in expand_dual_mappings(canonical):
            attributed = df.copy()
            attributed["specimen"]    = specimen
            attributed["source_file"] = fid
            frames.append(attributed)

    if not frames:
        raise RuntimeError("No drug-screen files could be read successfully")

    combined = pd.concat(frames, ignore_index=True)
    n_specimens = combined["specimen"].nunique()
    n_drugs     = combined["Drug"].nunique()
    logging.info(
        "Combined %d drug-screen rows across %d specimens × %d drugs",
        len(combined), n_specimens, n_drugs,
    )
    return combined


# ----------------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------------
def load_sample_map(samples_path: str) -> dict[str, int]:
    df = pd.read_csv(samples_path)
    df["improve_sample_id"] = df["improve_sample_id"].astype(int)
    return dict(zip(df["other_id"].astype(str), df["improve_sample_id"]))


def load_drug_map(drugs_path: str) -> dict[str, str]:
    df = pd.read_csv(drugs_path, sep="\t")
    if "chem_name" not in df.columns or "improve_drug_id" not in df.columns:
        raise KeyError(
            "drugs file missing required columns "
            f"(have: {list(df.columns)})"
        )
    mapping = {
        str(name).strip().lower(): str(did)
        for name, did in zip(df["chem_name"], df["improve_drug_id"])
    }
    logging.info("Loaded %d drug name → improve_drug_id mappings", len(mapping))
    return mapping


def lookup_drug_id(name, drug_map: dict[str, str]) -> str | None:
    if pd.isna(name):
        return None
    return drug_map.get(str(name).strip().lower())


# ----------------------------------------------------------------------------
# Single vs multi dose split
# ----------------------------------------------------------------------------
def split_single_vs_multi(combined: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per (specimen, drug), separate single-dose vs multi-dose."""
    counts = (
        combined.groupby(["specimen", "Drug"])["Concentration_uM"]
        .nunique()
        .reset_index(name="n_conc")
    )
    multi_keys  = counts[counts["n_conc"] > 1]
    single_keys = counts[counts["n_conc"] == 1]

    multi  = combined.merge(multi_keys[["specimen", "Drug"]],
                              on=["specimen", "Drug"], how="inner")
    single = combined.merge(single_keys[["specimen", "Drug"]],
                              on=["specimen", "Drug"], how="inner")
    # Single-dose rows must be at the canonical 1 μM concentration.
    single = single[single["Concentration_uM"] == SINGLE_DOSE_UM].copy()

    logging.info("Multi-dose rows: %d (across %d sample-drug pairs)",
                 len(multi), len(multi_keys))
    logging.info("Single-dose rows: %d (across %d sample-drug pairs at %g μM)",
                 len(single), len(single_keys), SINGLE_DOSE_UM)
    return multi, single


# ----------------------------------------------------------------------------
# Curve fitting (multi-dose) and single-dose formatting
# ----------------------------------------------------------------------------
def run_curve_fit(multi: pd.DataFrame, sample_map, drug_map) -> pd.DataFrame:
    if multi.empty:
        return pd.DataFrame()

    multi = multi.copy()
    multi["improve_sample_id"] = multi["specimen"].map(sample_map)
    multi["improve_drug_id"]   = multi["Drug"].apply(
        lambda n: lookup_drug_id(n, drug_map))

    n_drop = (multi["improve_sample_id"].isna() |
               multi["improve_drug_id"].isna()).sum()
    if n_drop:
        logging.warning("Dropping %d multi-dose rows with unmapped sample/drug",
                        n_drop)
    multi = multi.dropna(subset=["improve_sample_id", "improve_drug_id"])
    if multi.empty:
        return pd.DataFrame()

    # fit_curve.py expects DOSE in μM and GROWTH as viability percentage
    multi["DOSE"]   = multi["Concentration_uM"].astype(float) + 1e-4
    multi["GROWTH"] = multi["Viability_percentage"].astype(float)
    multi["time"]   = TIME
    multi["time_unit"] = TIME_UNIT
    multi["study"]  = STUDY
    multi["source"] = SOURCE

    cols = ["DOSE", "GROWTH", "study", "source", "improve_sample_id",
            "Drug", "time", "time_unit"]
    multi_for_fit = multi[cols].copy()
    multi_for_fit["improve_sample_id"] = (
        multi_for_fit["improve_sample_id"].astype(int))

    with tempfile.TemporaryDirectory() as tmpdir:
        in_path    = os.path.join(tmpdir, "drug_response.tsv")
        out_prefix = os.path.join(tmpdir, "cnf_curve")
        multi_for_fit.to_csv(in_path, sep="\t", index=False)

        cmd = ["python", FIT_CURVE_SCRIPT,
               "--input", in_path, "--output", out_prefix]
        logging.info("Running fit_curve: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

        # fit_curve.py writes <prefix>.0 in the reference snippet's behavior
        out_path = out_prefix + ".0"
        if not os.path.exists(out_path):
            alt = out_prefix
            if os.path.exists(alt):
                out_path = alt
            else:
                raise FileNotFoundError(
                    f"fit_curve.py did not produce output at {out_path}"
                )
        fitted = pd.read_csv(out_path, sep="\t")

    # fit_curve.py renames Drug → improve_drug_id, so the column is present
    # but holds the raw drug name, not the SMI_ ID. Always remap through drug_map.
    if "improve_drug_id" in fitted.columns:
        name_col = "improve_drug_id"
    elif "Drug" in fitted.columns:
        name_col = "Drug"
    else:
        raise KeyError(
            f"fit_curve output has neither improve_drug_id nor Drug "
            f"column (have: {list(fitted.columns)})"
        )
    fitted["improve_drug_id"] = fitted[name_col].apply(
        lambda n: lookup_drug_id(n, drug_map))

    fitted = fitted.dropna(subset=["improve_drug_id"])
    n_metrics = (fitted["dose_response_metric"].nunique()
                  if "dose_response_metric" in fitted.columns else 0)
    logging.info("Curve-fit produced %d rows across %d metrics",
                 len(fitted), n_metrics)
    return fitted


def format_single_dose(single: pd.DataFrame, sample_map, drug_map) -> pd.DataFrame:
    if single.empty:
        return pd.DataFrame()

    df = single.copy()
    df["improve_sample_id"] = df["specimen"].map(sample_map)
    df["improve_drug_id"]   = df["Drug"].apply(
        lambda n: lookup_drug_id(n, drug_map))

    n_drop = (df["improve_sample_id"].isna() |
               df["improve_drug_id"].isna()).sum()
    if n_drop:
        logging.warning("Dropping %d single-dose rows with unmapped sample/drug",
                        n_drop)
    df = df.dropna(subset=["improve_sample_id", "improve_drug_id"]).copy()

    df["dose_response_metric"] = SINGLE_DOSE_METRIC
    df["dose_response_value"]  = df["Viability_percentage"].astype(float) / 100.0
    df["time"]      = TIME
    df["time_unit"] = TIME_UNIT
    df["study"]     = STUDY
    df["source"]    = SOURCE
    df["improve_sample_id"] = df["improve_sample_id"].astype(int)

    cols = ["source", "improve_sample_id", "improve_drug_id", "study",
            "time", "time_unit", "dose_response_metric", "dose_response_value"]
    df = df[cols]
    logging.info("Formatted %d single-dose rows", len(df))
    return df


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", help="cnf_samples.csv from build_samples")
    parser.add_argument("drugs",   help="cnf_drugs.tsv from build_drugs")
    parser.add_argument("--output",  default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--parents", default=",".join(DEFAULT_DRUG_SCREEN_PARENTS),
        help="Comma-delimited Synapse folder IDs to walk for viability files",
    )
    args = parser.parse_args()

    configure_logging()

    parents = [p.strip() for p in args.parents.split(",") if p.strip()]
    if not parents:
        logging.error("No parent folders provided")
        sys.exit(2)

    sample_map = load_sample_map(args.samples)
    drug_map   = load_drug_map(args.drugs)

    syn = synapseclient.Synapse()
    syn.login()

    combined = pull_screen_files(syn, parents)
    multi, single = split_single_vs_multi(combined)

    fitted    = run_curve_fit(multi, sample_map, drug_map)
    single_df = format_single_dose(single, sample_map, drug_map)

    final = pd.concat([fitted, single_df], ignore_index=True)
    cols = ["source", "improve_sample_id", "improve_drug_id", "study",
            "time", "time_unit", "dose_response_metric", "dose_response_value"]
    for c in cols:
        if c not in final.columns:
            final[c] = pd.NA
    final = final[cols]

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    final.to_csv(args.output, sep="\t", index=False)
    logging.info("Wrote %s (%d rows)", args.output, len(final))


if __name__ == "__main__":
    sys.exit(main() or 0)
