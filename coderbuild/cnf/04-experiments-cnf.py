#!/usr/bin/env python3
"""
build_exp.py — generate cnf_experiments.tsv.

For every (sample, drug) pair in the cNF drug screen, this script produces
one or more experiment rows:

  * Drugs measured at multiple concentrations are run through
    coderbuild/utils/fit_curve.py to produce fit_auc / fit_ic50 / fit_einf /
    fit_hs / fit_r2 metrics.
  * Drugs measured at a single concentration (1 μM) are kept as-is with
    metric = 'uM_viability' and value = Viability_percentage / 100.
    NOTE: 'uM_viability' must be added to the ResponseMetric enum in
    schema/coderdata.yaml. See README and the schema patch note.

Constants:
    time = 120, time_unit = 'hours', study = 'cnf', source = 'Synapse:syn69947322'

Usage:
    python build_exp.py <cnf_samples.csv> <cnf_drugs.tsv>
                        [--output cnf_experiments.tsv]
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile

import pandas as pd
import synapseclient


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
DRUG_SCREEN_TABLE = "syn51301431"
DRUG_SCREEN_PARENT = "syn69947322"     # parent for source attribution
TIME              = 120
TIME_UNIT         = "hours"
STUDY             = "cnf"
SOURCE            = f"Synapse:{DRUG_SCREEN_PARENT}"
SINGLE_DOSE_UM    = 1.0                # treat 1 μM as the canonical single-dose
SINGLE_DOSE_METRIC = "uM_viability"
FIT_CURVE_SCRIPT  = os.environ.get(
    "CODERDATA_FIT_CURVE_SCRIPT", "/coderbuild/utils/fit_curve.py")


def configure_logging():
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


# ----------------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------------
def load_sample_map(samples_path: str) -> dict[str, int]:
    """Build other_id (specimen string) → improve_sample_id (int) lookup."""
    df = pd.read_csv(samples_path)
    df["improve_sample_id"] = df["improve_sample_id"].astype(int)
    return dict(zip(df["other_id"].astype(str), df["improve_sample_id"]))


def load_drug_map(drugs_path: str) -> dict[str, str]:
    """Build chem_name → improve_drug_id lookup. Names are case-insensitive."""
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


def lookup_drug_id(name: str, drug_map: dict[str, str]) -> str | None:
    if pd.isna(name):
        return None
    return drug_map.get(str(name).strip().lower())


# ----------------------------------------------------------------------------
# Pull and split drug-screen data
# ----------------------------------------------------------------------------
def pull_screen_files(syn) -> pd.DataFrame:
    """Concatenate every drug-screen file with specimenID attached as 'specimen'."""
    query = (
        "select id, individualID, specimenID "
        f"from {DRUG_SCREEN_TABLE} "
        "where dataType='drug screen'"
    )
    file_index = syn.tableQuery(query).asDataFrame()
    logging.info("Found %d drug-screen files", len(file_index))

    frames = []
    for _, row in file_index.iterrows():
        try:
            df = pd.read_csv(syn.get(row["id"]).path)
        except Exception as exc:
            logging.warning("Skipping file %s: %s", row["id"], exc)
            continue
        df["specimen"] = row["specimenID"]
        frames.append(df)

    if not frames:
        raise RuntimeError("No drug-screen files could be read")

    combined = pd.concat(frames, ignore_index=True)
    logging.info("Combined %d drug-screen rows", len(combined))
    return combined


def split_single_vs_multi(combined: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per (specimen, drug), separate single-dose (n_conc==1) vs multi-dose (n_conc>1)."""
    # Count concentrations per (specimen, drug)
    counts = (
        combined.groupby(["specimen", "Drug"])["Concentration_uM"]
        .nunique()
        .reset_index(name="n_conc")
    )
    multi_keys = counts[counts["n_conc"] > 1]
    single_keys = counts[counts["n_conc"] == 1]

    multi = combined.merge(multi_keys[["specimen", "Drug"]],
                            on=["specimen", "Drug"], how="inner")
    single = combined.merge(single_keys[["specimen", "Drug"]],
                             on=["specimen", "Drug"], how="inner")
    # Single-dose rows must be at the canonical single-dose concentration.
    single = single[single["Concentration_uM"] == SINGLE_DOSE_UM].copy()
    logging.info("Multi-dose rows: %d (across %d sample-drug pairs)",
                 len(multi), len(multi_keys))
    logging.info("Single-dose rows: %d (across %d sample-drug pairs)",
                 len(single), len(single_keys))
    return multi, single


# ----------------------------------------------------------------------------
# Curve fitting (multi-dose) and single-dose formatting
# ----------------------------------------------------------------------------
def run_curve_fit(multi: pd.DataFrame, sample_map, drug_map) -> pd.DataFrame:
    """Write a temp DOSE/GROWTH TSV, call fit_curve.py, read its output."""
    if multi.empty:
        return pd.DataFrame()

    # Map specimen → improve_sample_id and Drug → improve_drug_id
    multi = multi.copy()
    multi["improve_sample_id"] = multi["specimen"].map(sample_map)
    multi["improve_drug_id"]   = multi["Drug"].apply(
        lambda n: lookup_drug_id(n, drug_map))

    n_drop = (multi["improve_sample_id"].isna() |
              multi["improve_drug_id"].isna()).sum()
    if n_drop:
        logging.warning("Dropping %d multi-dose rows with unmapped sample or drug",
                        n_drop)
    multi = multi.dropna(subset=["improve_sample_id", "improve_drug_id"])

    # fit_curve.py expects DOSE in μM (or whatever the unit it documents) and
    # GROWTH in fractional viability. The reference snippet uses
    # GROWTH = Viability_percentage; matching that here for compatibility,
    # plus a small offset on DOSE so log-transform is well-defined at zero.
    multi["DOSE"]  = multi["Concentration_uM"].astype(float) + 1e-4
    multi["GROWTH"] = multi["Viability_percentage"].astype(float)
    multi["time"] = TIME
    multi["time_unit"] = TIME_UNIT
    multi["study"] = STUDY
    multi["source"] = SOURCE
    multi = multi.rename(columns={"Drug": "Drug"})  # explicit, kept

    cols = ["DOSE", "GROWTH", "study", "source", "improve_sample_id",
            "Drug", "time", "time_unit"]
    multi_for_fit = multi[cols].copy()
    multi_for_fit["improve_sample_id"] = multi_for_fit["improve_sample_id"].astype(int)

    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = os.path.join(tmpdir, "drug_response.tsv")
        out_prefix = os.path.join(tmpdir, "cnf_curve")
        multi_for_fit.to_csv(in_path, sep="\t", index=False)

        cmd = ["python", FIT_CURVE_SCRIPT,
               "--input", in_path, "--output", out_prefix]
        logging.info("Running fit_curve: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

        # fit_curve.py writes <prefix>.0 (long format) per the user reference snippet.
        out_path = out_prefix + ".0"
        if not os.path.exists(out_path):
            # try alternative output naming
            alt = out_prefix
            if os.path.exists(alt):
                out_path = alt
            else:
                raise FileNotFoundError(
                    f"fit_curve.py did not produce expected output at {out_path}"
                )
        fitted = pd.read_csv(out_path, sep="\t")

    # The fitter outputs in long format; remap Drug → improve_drug_id
    if "improve_drug_id" not in fitted.columns:
        if "Drug" in fitted.columns:
            fitted["improve_drug_id"] = fitted["Drug"].apply(
                lambda n: lookup_drug_id(n, drug_map))
        else:
            raise KeyError(
                f"fit_curve output has neither improve_drug_id nor Drug column "
                f"(have: {list(fitted.columns)})"
            )

    fitted = fitted.dropna(subset=["improve_drug_id"])
    logging.info("Curve-fit produced %d rows across %d metrics",
                 len(fitted),
                 fitted["dose_response_metric"].nunique() if
                 "dose_response_metric" in fitted.columns else 0)
    return fitted


def format_single_dose(single: pd.DataFrame, sample_map, drug_map) -> pd.DataFrame:
    """Format single-dose rows with metric=uM_viability."""
    if single.empty:
        return pd.DataFrame()

    df = single.copy()
    df["improve_sample_id"] = df["specimen"].map(sample_map)
    df["improve_drug_id"]   = df["Drug"].apply(
        lambda n: lookup_drug_id(n, drug_map))

    n_drop = (df["improve_sample_id"].isna() |
              df["improve_drug_id"].isna()).sum()
    if n_drop:
        logging.warning("Dropping %d single-dose rows with unmapped sample or drug",
                        n_drop)
    df = df.dropna(subset=["improve_sample_id", "improve_drug_id"]).copy()

    df["dose_response_metric"] = SINGLE_DOSE_METRIC
    df["dose_response_value"]  = df["Viability_percentage"].astype(float) / 100.0
    df["time"] = TIME
    df["time_unit"] = TIME_UNIT
    df["study"] = STUDY
    df["source"] = SOURCE
    df["improve_sample_id"] = df["improve_sample_id"].astype(int)

    cols = ["source", "improve_sample_id", "improve_drug_id", "study",
            "time", "time_unit", "dose_response_metric", "dose_response_value"]
    df = df[cols]
    logging.info("Formatted %d single-dose rows", len(df))
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", help="cnf_samples.csv from build_samples")
    parser.add_argument("drugs", help="cnf_drugs.tsv from build_drugs")
    parser.add_argument("--output", default="/tmp/cnf_experiments.tsv")
    args = parser.parse_args()

    configure_logging()

    sample_map = load_sample_map(args.samples)
    drug_map   = load_drug_map(args.drugs)

    syn = synapseclient.Synapse()
    syn.login()

    combined = pull_screen_files(syn)
    multi, single = split_single_vs_multi(combined)

    fitted     = run_curve_fit(multi, sample_map, drug_map)
    single_df  = format_single_dose(single, sample_map, drug_map)

    final = pd.concat([fitted, single_df], ignore_index=True)
    # Make sure required columns are present and ordered
    cols = ["source", "improve_sample_id", "improve_drug_id", "study",
            "time", "time_unit", "dose_response_metric", "dose_response_value"]
    for c in cols:
        if c not in final.columns:
            final[c] = pd.NA
    final = final[cols]

    final.to_csv(args.output, sep="\t", index=False)
    logging.info("Wrote %s (%d rows)", args.output, len(final))


if __name__ == "__main__":
    sys.exit(main() or 0)