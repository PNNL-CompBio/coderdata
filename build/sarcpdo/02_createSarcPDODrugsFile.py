import synapseclient
import pandas as pd
import numpy as np
import argparse
import os

from pubchem_retrieval import update_dataframe_and_write_tsv


def create_sarcpdo_drugs_file(synObject, prevDrugFilepath, outputPath):
    drug_query = synObject.tableQuery("select * from syn61892224")
    drug_data = drug_query.asDataFrame()
    # check status of previous drug file
    if not prevDrugFilepath:
        # if sarcpdo_drugs.tsv is null, create the empty dataframe. 
        empty_drugs = pd.DataFrame(columns = ['improve_drug_id', 'chem_name', 'pubchem_id', 'canSMILES', 'InChIKey', 'formula', 'weight'])
        empty_drugs.to_csv('outputPath', sep='\t', index=False)
    
    # get unique drugs 
    unique_drugs = drug_data['Drug_Name'].unique()
    # use helper functions in pubchem_retrieval.py 
    update_dataframe_and_write_tsv(unique_drugs, output_filename=outputPath, # specify ignore_chems as null?
                                   batch_size=1, isname=True, time_limit=48 * 60 * 60)


if __name__ == "__main__":
    print('in main')
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug data files for the Sarcoma PDO project")
    parser.add_argument('-d', '--prevDrugFilePath', help='Path to a previous drug file for sarcpdo',default=None)
    parser.add_argument('-o', '--outputPath', help='Output path for updated sarcpdo drug file', default = None) 
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)

    create_sarcpdo_drugs_file(synObject, args.prevDrugFilePath, args.outputPath)

    # command line testing: python3 02_createSarcPDODrugsFile.py -t $SYNAPSE_AUTH_TOKEN -d ../../../sarcpdo_drugs.csv -o sarcpdo_drugs.csv