import synapseclient
import pandas as pd
import numpy as np
import argparse
import os
import pubchem_retrieval as pr


def create_novartis_pdx_drugs_file(synObject, prevDrugFilepath, outputPath):
    file = synObject.get('syn66276102')
    # read raw drug data from synapse
    rawDrugData = pd.read_csv(file.path)
    # split on + operator - there are 2- and one 3- way drug combos in this dataset
    sepDrugNames = pd.Series(rawDrugData['Treatment'].unique()).str.split("+", expand=True)
  
    # taking the drug names from the first and second column from the split - there is only one 
    # drug name in the 3rd column (onen 3-way combo) that is replicated in other treatments as well
    alldrugnames = pd.Series(pd.concat([sepDrugNames[0], sepDrugNames[1]]).dropna()).str.split('"', expand=True)[0].str.split("-", expand=True)[0]
    #nodoseinfo = pd.Series(alldrugnames.str.split("-", expand =True)[0])
    #combineddrugames = pd.concat([alldrugnames, nodoseinfo])
    finalDrugNames = pd.Series(alldrugnames.unique()).str.strip().unique()
    # get unique drugs
    newdrugnames = finalDrugNames[finalDrugNames != 'untreated']
    
    raw_names = {str(n).strip() for n in finalDrugNames if str(n).strip().lower() != 'untreated'}


    print(f"Novartis PDX raw drug names ({len(raw_names)}): {sorted(raw_names)}")

    # delegate to unified PubChem retrieval, restricting to only these drugs
    final_df = pr.update_dataframe_and_write_tsv(
        unique_names=raw_names,
        output_filename=outputPath,
        batch_size=50,
        isname=True,
        prev_drug_filepaths=prevDrugFilepath if prevDrugFilepath and str(prevDrugFilepath).strip() else None,
        restrict_to_raw_names=raw_names
    )
    
    if final_df is None or final_df.empty:
        print("Warning: no Novartis PDX drugs were found.")
    else:
        kept_ids = set(final_df.get('improve_drug_id', []))
        print(f"Retained {len(final_df)} rows across {len(kept_ids)} improve_drug_id(s).")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug data files for the Novartis PDX data")
    parser.add_argument('-d', '--prevDrugFilePath', help='Path to a previous drug file for bladderpdo', nargs="?", default = None)
    parser.add_argument('-o', '--outputPath', help='Output path for updated novartis drug file', default = "/tmp/novartis_drugs.tsv") 
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    print("after PAT assignment")
    synObject = synapseclient.login(authToken=PAT)
    print('after creating synObject')
    if args.prevDrugFilePath:
        previousDrugs = args.prevDrugFilePath
    else:
        previousDrugs = None
    create_novartis_pdx_drugs_file(synObject, previousDrugs, args.outputPath)