import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 
from pubchem_retrieval import update_dataframe_and_write_tsv

### get drug data
def download_synapse_data(synID:str, save_path:str = None, synToken:str = None):
    """ 
    Download drug data from Synapse. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    
    Parameters
    ----------
    synID : string
        SynapseID of dataset to download.
        
    save_path : string
        Local path where the downloaded file will be saved.

    synToken : string
        Synapse Personal Access Token of user.  Requires a Synapse account. More information at: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
        
    Returns
    -------
    str
        Filepath to downloaded excel file
    """   
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    downloaded_data = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    data_filepath = downloaded_data.path
    return(data_filepath)


### create drug csv
def create_crc_drug_data(fitted_drug_data_path:str, prevDrugFilepath:str, output_drug_data_path:str):
    # import fitted drug data and get drug names from DRUG_NAME column
    fitted_drug_df = pd.read_csv(fitted_drug_data_path)
    crc_drugs_df = pd.DataFrame(fitted_drug_df['DRUG_NAME'].unique())
    # if there is a prev drug file, check for new drugs
    if prevDrugFilepath is not None and prevDrugFilepath is not "":
        prev_drug_df = pd.read_csv(prevDrugFilepath)
        # get drugs that are only in the crc_drugs_df (aka new drugs only)
        new_drugs_df = crc_drugs_df[~crc_drugs_df.chem_name.isin(prev_drug_df.chem_name)]
    else:
        # if there's no prev drugs, then all drugs are new
        new_drugs_df = crc_drugs_df
    # get new drug names
    new_drug_names = new_drugs_df['chem_name'].unique()
    # call function that gets info for these drugs
    update_dataframe_and_write_tsv(new_drug_names,output_drug_data_path)


############################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    # arguments for what data to process
    parser.add_argument('-d', '--Download', action = 'store_true', default=False, help='Download drug data.')
    parser.add_argument('-t', '--Token', type=str, default=None, help='Synapse Token')
    parser.add_argument('-D', '--Drug', action = 'store_true', default=False, help='Generate drug data.')
    parser.add_argument('-p', '--PrevDrugs', type=str, default=None, help='Synapse Token')

    args = parser.parse_args()


    ###########################

    if args.Download:
        if args.Token is None:
            print("No synpase download tocken was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # download fitted and raw drug data from synapse
            fitted_drug_data_path = download_synapse_data(synID = "syn65452841", save_path = "/tmp/", synToken = args.Token)
    if args.Drug:
        if args.PrevDrugs is None or args.PrevDrugs=='':
            print("No previous drugs file provided.  Starting improve_drug_id from SMI_1. Running drug file generation")
            create_crc_drug_data(fitted_drug_data_path = "/tmp/fitted_data_GDSC_Org_restricted_11Mar25.csv", output_drug_data_path = "/tmp/crc_drugs.csv")
        else:
            print("Previous drugs file {} detected. Running drugs file generation and checking for duplicate IDs.".format(args.PrevDrugs))
            create_crc_drug_data(fitted_drug_data_path = "/tmp/fitted_data_GDSC_Org_restricted_11Mar25.csv", prevDrugFilepath = args.PrevDrugs, output_drug_data_path = "/tmp/crc_drugs.csv")

