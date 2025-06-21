import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 


def download_experiments_data(synID:str , save_path:str = None, synToken:str = None):
    """ 
    Download omics data from Synapse at synapseID syn66401303. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    
    Parameters
    ----------
    synID : string
        SynapseID of dataset to download. Default is synapseID of the omics dataset.
        
    save_path : string
        atal path where the downloaded file will be saved.

    synToken : string
        Synapse Personal Access Token of user.  Requires a Synapse account. More information at: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
        
    Returns
    -------
    experiments_filepath : string
        Path to downloaded file

    """
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    downloaded_data = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    experiments_filepath = downloaded_data.path
    
    return(experiments_filepath)



def parse_experiments_excel_sheets(first_file_path, second_file_path):
    # read in the excel files
    first_exp_excel = pd.ExcelFile(open(first_file_path, 'rb'))
    first_experiments_dict = pd.read_excel(first_exp_excel, sheet_name=None, header=None)
    rest_exp_excel = pd.ExcelFile(open(second_file_path, 'rb'))
    rest_experiments_dict = pd.read_excel(rest_exp_excel, sheet_name=None, header=None)
    # use for loops to interate through the dictionaries, melt the df's into longer df's instead of matrices, and then concat
    list_of_exp_excels = [first_experiments_dict,rest_experiments_dict]
    full_df_list = []
    for dictionary in list_of_exp_excels:
        list_of_finished_dfs = []
        for experiment_key in dictionary.keys():
            one_sample_df = dictionary[experiment_key] # get 1 df from the df dictionary
            one_sample_df = one_sample_df.fillna(value={0:"concentration"}) # for many of the pages, they didn't write "concentration" but just left it blank. fill these na's with "concentration"
            list_of_dfs = [] # initiate empty list of df's for each conc type
            conc_indexes = one_sample_df[one_sample_df[0] == "concentration"].index.to_list() # get indexes of rows with concentrations in them (these will be column names)
            conc_indexes = conc_indexes + [one_sample_df.index[-1]+1]
            for index in range(0,(len(conc_indexes)-1)):
                one_conc_df = one_sample_df.loc[conc_indexes[index]:(conc_indexes[(index+1)]-1)]
                one_conc_df.columns = one_conc_df.iloc[0]
                one_conc_df = one_conc_df[1:]
                one_conc_df = pd.melt(one_conc_df, id_vars=['concentration'], value_vars=one_conc_df.columns[one_conc_df.columns != 'concentration'])
                one_conc_df = one_conc_df.rename(columns={"concentration":"drug_id",one_conc_df.columns[1]:"concentration","value":"count"})
                list_of_dfs.append(one_conc_df)
            elongated_df = pd.concat(list_of_dfs)
            elongated_df['sample_name'] = experiment_key
            list_of_finished_dfs.append(elongated_df)
        full_experiments_df = pd.concat(list_of_finished_dfs)
        full_df_list.append(full_experiments_df)
    experiments_df = pd.concat(full_df_list)
    return(experiments_df)