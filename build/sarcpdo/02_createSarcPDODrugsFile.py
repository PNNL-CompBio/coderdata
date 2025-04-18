import synapseclient
import pandas as pd
import numpy as np
import argparse
import os

#from utils.pubchem_retrieval import update_dataframe_and_write_tsv
from pubchem_retrieval import update_dataframe_and_write_tsv


def create_sarcpdo_drugs_file(synObject, prevDrugFilepath, outputPath):
    drug_query = synObject.tableQuery("select * from syn61892224")
    drug_data = drug_query.asDataFrame()
    
    # get unique drugs 
    newdrugnames = drug_data['Drug_Name'].unique()
    # use helper functions in pubchem_retrieval.py 
    alldrugs = []
    if prevDrugFilepath is not None and prevDrugFilepath is not "":
        prevdrugs = [pd.read_csv(t,sep='\t') for t in prevDrugFilepath.split(',')]
        alldrugs = pd.concat(prevdrugs).drop_duplicates()

        imps = alldrugs[alldrugs.chem_name.isin(newdrugnames)]
        newdrugs = alldrugs[alldrugs.improve_drug_id.isin(imps.improve_drug_id)]
        
        ##write drugs
        newdrugs.to_csv(outputPath, sep='\t', index=False)

    if len(alldrugs)==0 or len(newdrugnames)>len(set(newdrugs.improve_drug_id)): #we have more names we didn't match
        print('Missing drugs in existing file, querying pubchem')
        update_dataframe_and_write_tsv(newdrugnames,outputPath)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug data files for the Sarcoma PDO project")
    parser.add_argument('-d', '--prevDrugFilePath', help='Path to a previous drug file for sarcpdo', default = None)
    parser.add_argument('-o', '--outputPath', help='Output path for updated sarcpdo drug file', default = "/tmp/sarcpdo_drugs.tsv") 
    parser.add_argument('-t', '--token', help='Synapse token')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)

    create_sarcpdo_drugs_file(synObject, args.prevDrugFilePath, args.outputPath)

