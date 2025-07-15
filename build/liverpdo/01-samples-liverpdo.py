# import required packages
import pandas as pd
import numpy as np
import os
import gzip
import requests
import argparse
import synapseclient 



## Download the samples data from synapse
def download_samples_data(synID:str , save_path:str = None, synToken:str = None):
    """ 
    Download samples data from Synapse at synapseID syn64961953. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    
    Parameters
    ----------
    synID : string
        SynapseID of dataset to download. Default is synapseID of the samples dataset.
        
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
    samples_filepath = downloaded_data.path
    return(samples_filepath)

def map_substring(s, dict_map):
    for key in dict_map.keys():
        if key in s: 
            return dict_map[key]
    return np.nan

 ### create sample sheet function 
def generate_sample_file(samples_data_path:str = None, prev_samples_path:str = "") -> pd.DataFrame:
    """
    Creates sample file from samples data excel file. Checks the input sample file against previous sample files to make sure 
    there are no clashing sample names and assigns improved ID's starting from where previous sample sheet left off.
    
    Parameters
    ----------
    samples_data_path : string
        Path to samples data from https://pmc.ncbi.nlm.nih.gov/articles/PMC10949980/#_ad93_. Supplementary Table S1-S13
        
    prev_samples_path : string
        Path to previous sample sheet.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the combined samples data.
    
    """
    # reading in samples excel file
    samples_excel = pd.ExcelFile(open(samples_data_path, 'rb'))
    clinical_info_df = pd.read_excel(samples_excel, 'S1. Clinical information') # table with samples information 
    
    # reading in previous sample file 
    if prev_samples_path != "":
        prev_samples = pd.read_csv(prev_samples_path)

    # formatting table
    clinical_info_df = clinical_info_df.iloc[3:68] # keep only these rows, since the rest are formatted oddly and there's no info there
    samples_df = pd.DataFrame({'other_id':clinical_info_df.iloc[3:68,0]}).reset_index(drop=True)   
    samples_df['common_name'] = samples_df['other_id'] 
    ctype_dict = {'HCC':'hepatocellular carcinoma',
               'ICC':'intrahepatic cholangiocarcinoma',
               'CHC': 'combined hepatocellular-cholangiocarcinoma',
               'HB': 'hepatoblastoma'}
    samples_df['cancer_type'] = samples_df['other_id'].apply(lambda x: map_substring(x, ctype_dict))
    samples_df['other_id_source'] = "Synapse"
    samples_df['species'] = "Homo sapiens (Human)"
    samples_df['model_type'] = "patient derived organoid"


    # check other_id doesn't clash with previous sample names
    if prev_samples_path != "":
        if prev_samples.other_id.values in samples_df.other_id.values:
            print("Duplicate id names detected. Cannot proceed with generating sample sheet until resolved.")
            exit()
    if prev_samples_path == "":
        maxval = 0
    else:
        maxval = max(prev_samples.improve_sample_id)
    samples_df['improve_sample_id'] = samples_df.index + maxval + 1 # take index plus 1 to create counter, start from max value
    return(samples_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    parser.add_argument('-D', '--download',action='store_true', default=False, help='Download RNA seq and sequencing data from GEO and supplemental materials from https://www.cell.com/cell/fulltext/S0092-8674(15)00373-6#mmc2')
    parser.add_argument('-t', '--token', type=str, default=None, help='Synapse Token')
    parser.add_argument('-i', '--synapseID', type=str, default="syn66593307", help='SynapseID of data to download')

    parser.add_argument('-s', '--samples', action = 'store_true', help='Only generate samples, requires previous samples',default=False)
    parser.add_argument('-p', '--prevSamples', nargs='?',type=str, default='', const='', help='Use this to provide previous sample file')



    args = parser.parse_args()

    if args.download:
        if args.token is None:
            print("No synpase download tocken was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # Download samples data
            samples_download_path = download_samples_data(synID = args.synapseID, synToken = args.token, save_path = "/tmp")

    if args.samples:
        if args.prevSamples is None or args.prevSamples=='':
            print("No previous samples file provided.  Starting improve_sample_id from 1. Running sample file generation")
            sample_sheet = generate_sample_file(samples_data_path = samples_download_path)
        else:
            print("Previous sample sheet {} detected. Running sample file generation and checking for duplicate IDs.".format(args.prevSamples))
            sample_sheet = generate_sample_file(samples_data_path = samples_download_path, prev_samples_path= args.prevSamples)
        sample_sheet.to_csv("/tmp/liverpdo_samples.csv", index=False)
    

