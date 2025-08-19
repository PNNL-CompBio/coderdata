import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 
import pubchem_retrieval as pr
import warnings
warnings.filterwarnings("ignore")

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
def create_colorectal_drug_data(fitted_drug_data_path:str, prevDrugFilepath:str, output_drug_data_path:str):
    # import fitted drug data and get drug names from DRUG_NAME column
    fitted_drug_df = pd.read_csv(fitted_drug_data_path)
    raw_names = fitted_drug_df['DRUG_NAME'].unique()

    # prepare prev_drug_filepaths argument (None if empty)
    prev_arg = prevDrugFilepath if prevDrugFilepath and str(prevDrugFilepath).strip() != "" else None

    # call updated helper to fetch/merge and restrict to current Colorectal drugs
    final_df = pr.update_dataframe_and_write_tsv(
        unique_names=raw_names,
        output_filename=output_drug_data_path,
        batch_size=50,
        isname=True,
        prev_drug_filepaths=prev_arg,
        restrict_to_raw_names=raw_names
    )
    return final_df


############################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    # arguments for what data to process
    parser.add_argument('-d', '--Download', action = 'store_true', default=False, help='Download drug data.')
    parser.add_argument('-t', '--Token', type=str, default=None, help='Synapse Token')
    parser.add_argument('-D', '--Drug', action = 'store_true', default=False, help='Generate drug data.')
    parser.add_argument('-p', '--PrevDrugs', nargs='?', type=str, default='', const='', help='Previous drug file')

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
        prev_arg = args.PrevDrugs if args.PrevDrugs else None
        if not prev_arg:
            print("No previous drugs file provided. Starting improve_drug_id from SMI_1. Running drug file generation")
        else:
            print(f"Previous drugs file {args.PrevDrugs} detected. Running drugs file generation and checking for duplicate IDs.")
        create_colorectal_drug_data(
            fitted_drug_data_path="/tmp/fitted_data_GDSC_Org_restricted_11Mar25.csv",
            prevDrugFilepath=prev_arg if prev_arg is not None else "",
            output_drug_data_path="/tmp/colorectal_drugs.tsv"
        )