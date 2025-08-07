import synapseclient
import pandas as pd
import numpy as np
import argparse
import os
# for testing locally
#from utils.pubchem_retrieval import update_dataframe_and_write_tsv
# for building in docker
from pubchem_retrieval import update_dataframe_and_write_tsv


def create_bladder_pdo_drugs_file(synObject, prevDrugFilepath, outputPath):
    bladder_dir = synObject.get('syn64765430')
    filenames = list(synObject.getChildren(parent='syn64765430', includeTypes=['file']))

    # '-4' - there are 4 nondrug files in this directory.
    bladder_drugs = pd.DataFrame({'drugNames': [str]})
    for i in range(len(filenames) - 4):
        bladder_drugs.loc[i, 'drugNames'] = filenames[i]['name'].split(")")[1].split("(")[0].split(".")[0].strip()

    # get unique drug names
    raw_names = [str(n) for n in bladder_drugs['drugNames'].dropna().unique() if str(n).strip()]

    if not raw_names:
        print("No bladderPDO drug names extracted; exiting.")
        return

    print(f"BladderPDO raw drug names: {raw_names}")

    #New pubchem call
    update_dataframe_and_write_tsv(
        unique_names=raw_names,
        output_filename=outputPath,
        prev_drug_filepaths=prevDrugFilepath if prevDrugFilepath else None,
        isname=True,
        batch_size=50,
        restrict_to_raw_names=raw_names
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug data files for the Lee Bladder PDO project")
    parser.add_argument('-d', '--prevDrugFilePath', help='Path to a previous drug file for sarcpdo', nargs="?", default = None)
    parser.add_argument('-o', '--outputPath', help='Output path for updated sarcpdo drug file', default = "/tmp/sarcpdo_drugs.tsv") 
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)
    if args.prevDrugFilePath:
        previousDrugs = args.prevDrugFilePath
    else:
        previousDrugs = None
    create_bladder_pdo_drugs_file(synObject, previousDrugs, args.outputPath)