import pandas as pd
import os
import argparse
import time
import synapseclient as sc
import pubchem_retrieval as pr

###figshare link:

filelink='https://aacr.figshare.com/ndownloader/files/39996295'
#synid = 'syn64333325'
##get third tab and drugsa re listeda cross top

## Supplementary drug list (drug names are the column headers).
##
## Source: Tiriac H, Belleau P, Engle DD, Plenker D, Deschenes A, Somerville TDD,
## et al. "Organoid Profiling Identifies Common Responders to Chemotherapy in
## Pancreatic Cancer." Cancer Discovery (2018) 8(9):1112-1129.
## doi:10.1158/2159-8290.CD-18-0349
##
## Fetched from Figshare rather than the AACR CDN. The aacr.silverchair-cdn.com
## link that used to be hardcoded here is a CloudFront SIGNED url carrying an
## Expires timestamp, so it stops working on a fixed date and returns 403
## forever after. The one committed here expired 2025-01-27 and took the
## pancreatic drugs step down with it; a freshly generated replacement expires
## about five weeks after it is issued, so pasting in a new one only resets the
## clock. Figshare download links do not expire.
##
## Verified: this file has 4 sheets (Key, Chemo, Targeted, Targeted for
## Chem-refractory) and yields the same 26 drug names the AACR copy did.
TABLINK = 'https://ndownloader.figstatic.com/files/39996295'

## Number of drug names the supplementary table is expected to yield. If the
## file upstream is replaced or restructured the count moves, and silently
## building pancreatic from a different drug set would be worse than failing.
EXPECTED_DRUG_COUNT = 26

RETRY_SLEEPS = (60, 180, 600, 900)


def getDrugNames(token=""):
    last_error = None
    for attempt in range(len(RETRY_SLEEPS) + 1):
        try:
            # One fetch, both sheets. Reading the url twice downloaded this
            # 5.3MB workbook twice.
            with pd.ExcelFile(TABLINK) as book:
                ctab = pd.read_excel(book, sheet_name=1, skiprows=1)   # chemo
                ttab = pd.read_excel(book, sheet_name=2, skiprows=1)   # targeted
            break
        except Exception as e:                    # noqa: BLE001 - network shapes vary
            last_error = e
            if attempt < len(RETRY_SLEEPS):
                wait = RETRY_SLEEPS[attempt]
                print(f"  fetching {TABLINK} failed ({e}); attempt "
                      f"{attempt + 1}/{len(RETRY_SLEEPS) + 1}, retrying in "
                      f"{wait // 60} min", flush=True)
                time.sleep(wait)
    else:
        raise RuntimeError(
            f"Could not fetch the pancreatic supplementary drug table from "
            f"{TABLINK} after {len(RETRY_SLEEPS) + 1} attempts: {last_error}")

    drugs = [a.lower() for a in ctab.columns]+[a.lower() for a in ttab.columns]
    drugs = set(drugs)-set(['sample id','insensitive'])

    if len(drugs) != EXPECTED_DRUG_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_DRUG_COUNT} pancreatic drug names from the "
            f"supplementary table but found {len(drugs)}: {sorted(drugs)}. The "
            f"upstream file has changed; confirm it is still the Tiriac et al. "
            f"table and update EXPECTED_DRUG_COUNT deliberately.")
    return drugs


def main():
    parser = argparse.ArgumentParser(description='Download and match pancreatic drugs')
#    parser.add_argument('-p', '--pat',help='Synapse authentication token with permission to syn64333325')
    parser.add_argument('-d', '--prevDrugFile', default=None, help='Comma-delimited list of previous drug files')
    parser.add_argument('-o', '--output', default = '/tmp/pancreatic_drugs.tsv')

    args = parser.parse_args()
    newdrugnames = getDrugNames()
    print(f"Raw pancreatic drug names ({len(newdrugnames)}): {sorted(newdrugnames)}")

    final_df = pr.update_dataframe_and_write_tsv(
        unique_names=newdrugnames,
        output_filename=args.output,
        batch_size=50,
        isname=True,
        prev_drug_filepaths=args.prevDrugFile if args.prevDrugFile and args.prevDrugFile.strip() else None,
        restrict_to_raw_names=newdrugnames
    )

    if final_df.empty:
        print("Warning: no pancreatic drugs were found.")
    else:
        kept_ids = set(final_df.get('improve_drug_id', []))
        print(f"Retained {len(final_df)} rows across {len(kept_ids)} improve_drug_id(s).")
    
    
if __name__=='__main__':
    main()
