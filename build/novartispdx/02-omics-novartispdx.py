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


def map_copy_number_novPDX(copy_number_data, improve_id_data, entrez_data):
    """
    Maps copy number data to improved sample id's and entrez gene data. Also does some data formatting.
    
    Parameters
    ----------
    copy_number_data : pd.Dataframe OR string
        Pandas dataframe object with copy number data OR path to csv with copy number data

    improve_id_data : pd.Dataframe OR string
        Pandas dataframe object with improve id data OR path to csv with improve id data.  This is one of the outputs of parse_mmc2()

    entrez_data : pd.Dataframe OR string
        Pandas dataframe object with entrez gene data OR path to csv with entrez gene data.  Use this code to get this file: https://github.com/PNNL-CompBio/coderdata/tree/e65634b99d060136190ec5fba0b7798f8d140dfb/build/genes 

    Returns
    -------
    sample_entrez_cn_df : pd.DataFrame
        A DataFrame containing the mapped copy number data with columns: entrez_id, copy_number, copy_call, study, source ,improve_sample_id

    """
    # read in data
    if isinstance(copy_number_data, pd.DataFrame) == False:
        copy_number_data = pd.read_csv(copy_number_data)

    if isinstance(improve_id_data, pd.DataFrame) == False:
        improve_id_data = pd.read_csv(improve_id_data)
    
    if isinstance(entrez_data, pd.DataFrame) == False:
        entrez_data = pd.read_csv(entrez_data)
    
    # melt dataframe so that there is gene name and improve_sample_id per row
    long_cn_df = pd.melt(copy_number_data, id_vars=['Sample'], value_vars=copy_number_data.columns[copy_number_data.columns != 'Sample'])

    # get entrez id's from Sample
    entrez_cn_df = pd.merge(long_cn_df, entrez_data[['other_id','entrez_id']].drop_duplicates(), how = 'left', left_on= "Sample", right_on= "other_id")

    # get copy call from value column (aka copy number)
    entrez_cn_df['copy_call'] = [get_copy_call(a) for a in entrez_cn_df['value']]
    
    # get improve sample id
    improve_id_data['to_merge'] = improve_id_data['common_name'].str.replace("NIBR","")
    sample_entrez_cn_df = pd.merge(entrez_cn_df.drop_duplicates(), improve_id_data[['to_merge','improve_sample_id']].drop_duplicates(), how = 'left', left_on= "variable", right_on= "to_merge")

    # clean up columns and data types
    sample_entrez_cn_df = sample_entrez_cn_df.drop(columns=['Sample','variable','other_id','to_merge'])
    sample_entrez_cn_df['source'] = "CPDM"
    sample_entrez_cn_df['study'] = "novartispdx"
    sample_entrez_cn_df = sample_entrez_cn_df.rename(columns={'value':'copy_number'})
    sample_entrez_cn_df = sample_entrez_cn_df.astype({'entrez_id':'int','improve_sample_id':'int'})
    sample_entrez_cn_df = sample_entrez_cn_df[['entrez_id','copy_number','copy_call','study','source','improve_sample_id']]
    
    return(sample_entrez_cn_df)