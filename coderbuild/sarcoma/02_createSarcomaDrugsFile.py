import synapseclient
import pandas as pd
import numpy as np
import argparse
import os
import pubchem_retrieval as pr

def create_sarcoma_drugs_file(synObject, prevDrugFilepath, outputPath):
    drug_query = synObject.tableQuery("select * from syn61892224")
    drug_data = drug_query.asDataFrame()

    # get unique drugs
    raw_names = set([str(n) for n in drug_data['Drug_Name'].dropna().unique()])

    print(f"Sarcoma raw drug names: {raw_names}")

    final_df = pr.update_dataframe_and_write_tsv(
        unique_names=raw_names,
        output_filename=outputPath,
        batch_size=50,
        isname=True,
        prev_drug_filepaths=prevDrugFilepath if prevDrugFilepath else None,
        restrict_to_raw_names=raw_names
    )

    if final_df.empty:
        print("Warning: no Sarcoma drugs were found.")
    else:
        kept_ids = set(final_df.get('improve_drug_id', []))
        print(f"Retained {len(final_df)} rows across {len(kept_ids)} improve_drug_id(s).")

    return final_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug data files for the Sarcoma PDO project")
    parser.add_argument('-d', '--prevDrugFilePath', help='Path to a previous drug file for sarcoma', default = None)
    parser.add_argument('-o', '--outputPath', help='Output path for updated sarcoma drug file', default = "/tmp/sarcoma_drugs.tsv") 
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)

    create_sarcoma_drugs_file(synObject, args.prevDrugFilePath, args.outputPath)

