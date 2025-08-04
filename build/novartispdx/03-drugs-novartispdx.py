import synapseclient
import pandas as pd
import numpy as np
import argparse
import os
# for testing locally
from pubchem_retrieval import update_dataframe_and_write_tsv
# for building in docker
#from pubchem_retrieval import update_dataframe_and_write_tsv


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

    #print(finalDrugNames.tolist) 
    #newdrugnames = finalDrugNames.remove('untreated')
    print(2)
    print(newdrugnames)


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

    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug data files for the Novartis PDX data")
    parser.add_argument('-d', '--prevDrugFilePath', help='Path to a previous drug file for bladderpdo', nargs="?", default = None)
    parser.add_argument('-o', '--outputPath', help='Output path for updated novartispdx drug file', default = "/tmp/novartispdx_drugs.tsv") 
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