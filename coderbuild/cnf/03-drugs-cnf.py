#!/usr/bin/env python3
"""
03-drugs-cnf.py: generate cnf_drugs.tsv and cnf_drug_descriptors.tsv.

Pulls unique drug names from the cNF drug screen on Synapse, then calls the
standard coderdata utilities:

    coderbuild/utils/pubchem_retrieval.py
        -> cnf_drugs.tsv

    coderbuild/utils/build_drug_descriptor_table.py
        -> cnf_drug_descriptors.tsv

Important implementation note:
    pubchem_retrieval.py does not expose a command-line interface in the
    provided version. It defines update_dataframe_and_write_tsv(), but running
    `python pubchem_retrieval.py --input ... --output ...` exits 0 without
    doing any work. Therefore this script imports pubchem_retrieval.py directly
    and calls update_dataframe_and_write_tsv().

The descriptor builder is still run as a subprocess. It may call R/biomaRt
internally, which can fail when Ensembl is unavailable. To keep the drugs file
produced even when descriptors fail, this script:

  1. Retries external/remote work with exponential backoff.
  2. Sets a generous timeout.
  3. Supports --skip-descriptors.
  4. Supports --only-descriptors for retrying descriptor generation later.

Usage:
    python 03-drugs-cnf.py --prev_drugs <file1.tsv,file2.tsv,...>
                          [--out_drugs /tmp/cnf_drugs.tsv]
                          [--out_desc  /tmp/cnf_drug_descriptors.tsv]
                          [--skip-descriptors]
                          [--only-descriptors]
                          [--retries 3]
                          [--timeout 1800]
"""

import argparse
import importlib.util
import logging
import os
import signal
import subprocess
import sys
import tempfile
import time
import pandas as pd
import synapseclient


# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------

DRUG_SCREEN_TABLE = "syn51301431"

PUBCHEM_SCRIPT = os.environ.get(
    "CODERDATA_PUBCHEM_SCRIPT",
    "/coderbuild/utils/pubchem_retrieval.py",
)

DESCRIPTOR_SCRIPT = os.environ.get(
    "CODERDATA_DESCRIPTOR_SCRIPT",
    "/coderbuild/utils/build_drug_descriptor_table.py",
)

DEFAULT_OUT_DRUGS = "/tmp/cnf_drugs.tsv"
DEFAULT_OUT_DESC = "/tmp/cnf_drug_descriptors.tsv"

DEFAULT_RETRIES = 3
DEFAULT_TIMEOUT = 1800  # seconds per attempt
RETRY_BACKOFF_S = 30    # 30s, 60s, 120s, ...


# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------

def configure_logging() -> None:
    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


# ----------------------------------------------------------------------------
# Generic subprocess wrapper with retry
# ----------------------------------------------------------------------------

def run_with_retries(
    cmd: list[str],
    retries: int,
    timeout: int,
    label: str,
) -> bool:
    """
    Run a subprocess, retrying on failure with exponential backoff.

    Returns True on success, False if all retries are exhausted.
    """
    for attempt in range(1, retries + 1):
        logging.info("[%s] attempt %d/%d", label, attempt, retries)

        try:
            subprocess.run(cmd, check=True, timeout=timeout)
            logging.info("[%s] succeeded on attempt %d", label, attempt)
            return True

        except subprocess.TimeoutExpired:
            logging.warning(
                "[%s] timed out after %ds on attempt %d",
                label,
                timeout,
                attempt,
            )

        except subprocess.CalledProcessError as exc:
            logging.warning(
                "[%s] failed with exit %d on attempt %d",
                label,
                exc.returncode,
                attempt,
            )

        except Exception as exc:
            logging.warning(
                "[%s] unexpected error on attempt %d: %s",
                label,
                attempt,
                exc,
            )

        if attempt < retries:
            backoff = RETRY_BACKOFF_S * (2 ** (attempt - 1))
            logging.info("[%s] sleeping %ds before retry", label, backoff)
            time.sleep(backoff)

    logging.error("[%s] exhausted %d retries", label, retries)
    return False


# ----------------------------------------------------------------------------
# PubChem utility loading/calling
# ----------------------------------------------------------------------------

def load_pubchem_module(path: str):
    """
    Load pubchem_retrieval.py as a Python module.

    The provided pubchem_retrieval.py has no command-line entry point, so it
    must be imported and called directly.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"PubChem retrieval script not found: {path}")

    spec = importlib.util.spec_from_file_location(
        "coderdata_pubchem_retrieval",
        path,
    )

    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create import spec for {path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "update_dataframe_and_write_tsv"):
        raise AttributeError(
            f"{path} does not define update_dataframe_and_write_tsv()"
        )

    return module


def call_pubchem_retrieval(
    names: list[str],
    prev_drug_files: list[str],
    out_drugs: str,
    retries: int,
    timeout: int,
) -> bool:
    """
    Generate the drugs TSV by calling pubchem_retrieval.update_dataframe_and_write_tsv().

    This intentionally does not use subprocess because the supplied
    pubchem_retrieval.py has no argparse/main block and exits 0 without writing
    anything when run as a script.
    """
    prev_arg = ",".join(prev_drug_files) if prev_drug_files else None

    for attempt in range(1, retries + 1):
        logging.info("[pubchem_retrieval] attempt %d/%d", attempt, retries)

        ignore_chems = None

        try:
            pubchem = load_pubchem_module(PUBCHEM_SCRIPT)

            with tempfile.NamedTemporaryFile(
                "w",
                suffix="_cnf_ignore_chems.txt",
                delete=False,
            ) as tmp_ignore:
                ignore_chems = tmp_ignore.name

            df = pubchem.update_dataframe_and_write_tsv(
                unique_names=names,
                output_filename=out_drugs,
                ignore_chems=ignore_chems,
                batch_size=1,
                isname=True,
                time_limit=timeout,
                prev_drug_filepaths=prev_arg,
                restrict_to_raw_names=names,
            )

            # pubchem_retrieval.py sets signal.alarm(time_limit), but does not
            # clear it. Clear it here so the descriptor step is not interrupted.
            try:
                signal.alarm(0)
            except Exception:
                pass

            if getattr(pubchem, "should_continue", True) is False:
                raise TimeoutError(
                    f"pubchem_retrieval reached its time limit of {timeout}s"
                )

            if not os.path.exists(out_drugs):
                raise FileNotFoundError(
                    f"pubchem_retrieval returned but did not create {out_drugs}"
                )

            if os.path.getsize(out_drugs) == 0:
                raise RuntimeError(
                    f"pubchem_retrieval created {out_drugs}, but it is empty"
                )

            n_rows = len(df) if df is not None else "unknown"

            logging.info(
                "[pubchem_retrieval] succeeded on attempt %d; wrote %s rows to %s",
                attempt,
                n_rows,
                out_drugs,
            )

            return True

        except Exception as exc:
            try:
                signal.alarm(0)
            except Exception:
                pass

            logging.warning(
                "[pubchem_retrieval] failed on attempt %d/%d: %s",
                attempt,
                retries,
                exc,
            )

            if attempt < retries:
                backoff = RETRY_BACKOFF_S * (2 ** (attempt - 1))
                logging.info(
                    "[pubchem_retrieval] sleeping %ds before retry",
                    backoff,
                )
                time.sleep(backoff)

        finally:
            if ignore_chems and os.path.exists(ignore_chems):
                try:
                    os.unlink(ignore_chems)
                except OSError:
                    pass

    logging.error("[pubchem_retrieval] exhausted %d retries", retries)
    return False


# ----------------------------------------------------------------------------
# Drug name collection
# ----------------------------------------------------------------------------

def collect_drug_names(syn) -> list[str]:
    """
    Pull every cNF drug-screen file from Synapse and collect unique drug names.
    """
    query = (
        "select id, individualID, specimenID "
        f"from {DRUG_SCREEN_TABLE} "
        "where dataType='drug screen'"
    )

    files = syn.tableQuery(query).asDataFrame()

    logging.info("Reading %d drug-screen files for drug names", len(files))

    names: set[str] = set()

    for _, row in files.iterrows():
        try:
            entity = syn.get(row["id"])
            df = pd.read_csv(entity.path)

        except Exception as exc:
            logging.warning("Skipping file %s: %s", row["id"], exc)
            continue

        if "Drug" not in df.columns:
            logging.warning("File %s has no 'Drug' column; skipping", row["id"])
            continue

        names.update(
            df["Drug"]
            .dropna()
            .astype(str)
            .str.strip()
            .tolist()
        )

    names = sorted(n for n in names if n and n.lower() != "nan")

    logging.info("Collected %d unique drug names", len(names))

    return names


# ----------------------------------------------------------------------------
# Descriptor utility
# ----------------------------------------------------------------------------

def call_descriptor_table(
    drug_file: str,
    out_desc: str,
    retries: int,
    timeout: int,
) -> bool:
    """
    Run the descriptor table builder as a subprocess.

    This function assumes build_drug_descriptor_table.py has a CLI accepting
    --input and --output.
    """
    cmd = [
        "python",
        DESCRIPTOR_SCRIPT,
        "--input",
        drug_file,
        "--output",
        out_desc,
    ]

    ok = run_with_retries(
        cmd,
        retries,
        timeout,
        "descriptor_table",
    )

    if ok and not os.path.exists(out_desc):
        logging.error(
            "descriptor_table exited successfully but did not create %s",
            out_desc,
        )
        return False

    return ok


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--prev_drugs",
        default="",
        help="Comma-delimited list of existing drug TSV files for ID reuse.",
    )

    parser.add_argument(
        "--out_drugs",
        default=DEFAULT_OUT_DRUGS,
        help=f"Output drugs TSV. Default: {DEFAULT_OUT_DRUGS}",
    )

    parser.add_argument(
        "--out_desc",
        default=DEFAULT_OUT_DESC,
        help=f"Output drug descriptors TSV. Default: {DEFAULT_OUT_DESC}",
    )

    parser.add_argument(
        "--skip-descriptors",
        action="store_true",
        help=(
            "Skip the descriptor step entirely. Use when the descriptor "
            "utility or its biomaRt/Ensembl dependency is failing. The drugs "
            "table is still produced."
        ),
    )

    parser.add_argument(
        "--only-descriptors",
        action="store_true",
        help=(
            "Only run the descriptor step, expecting --out_drugs to already "
            "exist. Useful for retrying after a previous descriptor failure."
        ),
    )

    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help=f"Number of retry attempts. Default: {DEFAULT_RETRIES}",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"Timeout in seconds per attempt. Default: {DEFAULT_TIMEOUT}",
    )

    args = parser.parse_args()

    configure_logging()

    if args.skip_descriptors and args.only_descriptors:
        logging.error("Cannot combine --skip-descriptors with --only-descriptors")
        return 2

    for path in (args.out_drugs, args.out_desc):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    # ------------------------------------------------------------------------
    # Drugs step
    # ------------------------------------------------------------------------

    if not args.only_descriptors:
        syn = synapseclient.Synapse()
        syn.login()

        names = collect_drug_names(syn)

        if not names:
            logging.error("No drug names found; aborting")
            return 1

        prev = [
            p.strip()
            for p in args.prev_drugs.split(",")
            if p.strip()
        ]

        ok = call_pubchem_retrieval(
            names=names,
            prev_drug_files=prev,
            out_drugs=args.out_drugs,
            retries=args.retries,
            timeout=args.timeout,
        )

        if not ok:
            logging.error("pubchem_retrieval failed; aborting")
            return 1

        if not os.path.exists(args.out_drugs):
            logging.error(
                "pubchem_retrieval did not produce %s",
                args.out_drugs,
            )
            return 1

        logging.info("Wrote %s", args.out_drugs)

    # ------------------------------------------------------------------------
    # Descriptor step
    # ------------------------------------------------------------------------

    if args.skip_descriptors:
        logging.warning(
            "Skipping descriptor step. Re-run with --only-descriptors to "
            "produce %s when descriptor dependencies are reachable.",
            args.out_desc,
        )
        return 0

    if not os.path.exists(args.out_drugs):
        logging.error(
            "Drugs file %s missing; cannot run descriptors",
            args.out_drugs,
        )
        return 1

    ok = call_descriptor_table(
        drug_file=args.out_drugs,
        out_desc=args.out_desc,
        retries=args.retries,
        timeout=args.timeout,
    )

    if not ok:
        logging.error(
            "Descriptor step failed after %d attempts. The drugs file %s is "
            "complete; re-run with --only-descriptors when descriptor "
            "dependencies are reachable to produce %s.",
            args.retries,
            args.out_drugs,
            args.out_desc,
        )

        # Exit 0 because the primary drugs table was produced. Descriptors can
        # be retried later. Use --skip-descriptors to suppress this explicitly.
        return 0

    logging.info("Wrote %s", args.out_desc)

    return 0


if __name__ == "__main__":
    sys.exit(main())
