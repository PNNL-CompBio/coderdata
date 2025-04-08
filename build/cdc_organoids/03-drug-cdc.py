import pandas as pd
import numpy as np
import os
import math
import argparse

### get drug data
def download_synapse_data(synID:str , save_path:str = None, synToken:str = None):
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
def make_drug_data():


    return()


###


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    # arguments for what data to process
    parser.add_argument('-D', '--download', action = 'store_true', default=False, help='Download drug data.')
    parser.add_argument('-t', '--token', type=str, default=None, help='Synapse Token')

    args = parser.parse_args()


    ###########################

    if args.download:
        if args.token is None:
            print("No synpase download tocken was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # download fitted and raw drug data from synapse
            fitted_drug_data_path = download_synapse_data(synID = fitted_drug_data_synID, save_path = "/tmp/fitted_drug_data.csv", synToken = "syn65452841")
            raw_drug_data_path = download_synapse_data(synID = raw_drug_data_synID, save_path = "/tmp/raw_drug_data.csv", synToken = "syn65452842")

 