import pandas as pd
import numpy as np
import os
import math
import argparse

def download_parse_omics_novPDX(synID:str , save_path:str = None, synToken:str = None):
    """ 
    Download omics data from Synapse at synapseID syn66364488. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    Omics data is an excel file.  The excel file is then parsed for the RNAseq, copy number, and mutations data.
    
    Parameters
    ----------
    synID : string
        SynapseID of dataset to download. Default is synapseID of the sequencing dataset.
        
    save_path : string
        Local path where the downloaded file will be saved.

    synToken : string
        Synapse Personal Access Token of user.  Requires a Synapse account. More information at: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
        
    Returns
    -------
    mutations_data : pd.DataFrame
        A DataFrame containing mutations data.

    copy_number_data : pd.DataFrame
        A DataFrame containing copy number data.

    rnaseq_data : pd.DataFrame
        A DataFrame containing RNAseq data.
    """

    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    syn66364488 = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    sequencing_filepath = syn66364488.path
    all_omics_excel = pd.ExcelFile(open(sequencing_filepath, 'rb'))
    mutations_data = pd.read_excel(all_omics_excel, 'pdxe_mut_and_cn2') # table with somatic mutation information 
    copy_number_data = pd.read_excel(all_omics_excel, 'copy number') # table with copy number information
    rnaseq_data = pd.read_excel(all_omics_excel, 'RNAseq_fpkm')


    return(rnaseq_data, copy_number_data, mutations_data)