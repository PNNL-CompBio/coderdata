import pandas as pd
import numpy as np
import wget
import os
import gzip
import requests
import argparse
import synapseclient 

###### NOTES ######
#   * need to change all paths to paths relevant to docker image
#   * add description to parser
#   * run functions in ipynb to test they are working

def download_rnaseq(geo_url:str = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65253&format=file&file=GSE65253%5Fcol%5Ftum%5Forg%5Fmerge%2Ecsv%2Egz", save_path:str = None): 
    """
    Retrieve data from a given GEO URL and identify the downloaded file by its name.

    This function uses the wget tool to download a file from the provided GEO URL.
    By comparing the directory contents before and after the download, 
    it identifies the newly downloaded file's name.

    Parameters
    ----------
    geo_url : str
        The GEO URL pointing to the data to be downloaded. Default is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65253
    
    save_path : string
        Local path where the downloaded file will be saved.

    Returns
    -------
    None
    """
    
    response = requests.get(geo_url)
    with open(save_path, 'wb') as f:
        f.write(response.content)
    return

def download_sequencing_data(synID:str , save_path:str = None, synToken:str = None):
    """ 
    Download sequencing data from Synapse at synapseID syn64961953. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    
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
    str
        Filepath to downloaded excel file
    """
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    syn64961953 = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    sequencing_filepath = syn64961953.path
    return sequencing_filepath

def generate_sample_file(sequencing_data_path:str = None, prev_samples_path:str = "") -> pd.DataFrame:
    """
    Creates sample file from sequencing data excel file. Checks the input sample file against previous sample files to make sure 
    there are no clashing sample names and assigns improved ID's starting from where previous sample sheet left off.
    
    Parameters
    ----------
    sequencing_data_path : string
        Path to sequencing data from https://www.cell.com/cell/fulltext/S0092-8674(15)00373-6#sec-4 . Supplementary Table S1
        
    prev_samples_path : string
        Path to previous sample sheet.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the combined samples data.
    
    """
    # reading in sequencing excel file
    sequencing_excel = pd.ExcelFile(open(sequencing_data_path, 'rb'))
    recurrent_mutations = pd.read_excel(sequencing_excel, 'TableS1I_Recurrent mutations') # table with recurrent mutation information 
    somatic_mutations = pd.read_excel(sequencing_excel, 'TableS1J-Somatic mutations') # table with somatic mutation information 
    
    # reading in previous sample file 
    if prev_samples_path != "":
        prev_samples = pd.read_csv(prev_samples_path)
    recurrent_tumor = pd.DataFrame({'other_id':recurrent_mutations['Tumor_Sample_Barcode'].str.split('-',n = 1,expand=True).iloc[:,1].unique()})
    recurrent_normal = pd.DataFrame({'other_id':recurrent_mutations['Matched_Norm_Sample_Barcode'].str.split('-',n = 1,expand=True).iloc[:,1].unique()})
    
    # merging somatic organoids too just in case recurrent excludes some
    somatic_tumor = pd.DataFrame({'other_id':somatic_mutations['Tumor_Sample_Barcode'].str.split('-',n = 1,expand=True).iloc[:,1].unique()})
    somatic_normal = pd.DataFrame({'other_id':somatic_mutations['Matched_Norm_Sample_Barcode'].str.split('-',n = 1,expand=True).iloc[:,1].unique()})
    samples_df = pd.concat([recurrent_tumor,recurrent_normal, somatic_tumor, somatic_normal])
    
    # formatting the table
    samples_df = samples_df.drop_duplicates('other_id')
    samples_df = samples_df.reset_index()
    samples_df['common_name'] = samples_df['other_id'].str.split('-', n = 1,expand=True).iloc[:,0] + "-"
    samples_df['model_type'] = ""
    for index, row in samples_df.iterrows():
        if "Tumor-Organoid" in samples_df.loc[index, 'other_id']:
            samples_df.loc[index, 'common_name'] = samples_df.loc[index, 'common_name'] + "T-O"
            samples_df.loc[index, 'model_type'] = "organoid"
        if "Tumor-Biopsy" in samples_df.loc[index, 'other_id']:
            samples_df.loc[index, 'common_name'] = samples_df.loc[index, 'common_name'] + "T-B"
            samples_df.loc[index, 'model_type'] = "biopsy"
        if "Normal-Organoid" in samples_df.loc[index, 'other_id']:
            samples_df.loc[index, 'common_name'] = samples_df.loc[index, 'common_name'] + "N-O"
            samples_df.loc[index, 'model_type'] = "organoid"
    samples_df['other_id_source'] = "vandeWetering_2015"
    samples_df['cancer_type'] = "Colorectal Carcinoma"
    samples_df['species'] = "Homo sapiens (Human)"

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
    samples_df = samples_df.drop(columns = 'index')
    return(samples_df)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    parser.add_argument('-D', '--download',action='store_true', default=False, help='Download RNA seq and sequencing data from GEO and supplemental materials from https://www.cell.com/cell/fulltext/S0092-8674(15)00373-6#mmc2')
    parser.add_argument('-t', '--token', type=str, default=None, help='Synapse Token')
    parser.add_argument('-i', '--synapseID', type=str, default="syn64961953", help='SynapseID of data to download')

    parser.add_argument('-s', '--samples', action = 'store_true', help='Only generate samples, requires previous samples',default=False)
    parser.add_argument('-p', '--prevSamples', nargs='?',type=str, default='', const='', help='Use this to provide previous sample file')



    args = parser.parse_args()


    ###########################

    if args.download:
        if args.token is None:
            print("No synpase download tocken was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # Download RNA seq data
            download_rnaseq(save_path = "/tmp/GSE65253_col_tum_org_merge.csv.gz")
            # Download sequencing data
            sequencing_download_path = download_sequencing_data(synID = args.synapseID, synToken = args.token, save_path = "/tmp/GSE65253_col_tum_org_merge.csv.gz")

    if args.samples:
        if args.prevSamples is None or args.prevSamples=='':
            print("No previous samples file provided.  Starting improve_sample_id from 1. Running sample file generation")
            sample_sheet = generate_sample_file(sequencing_data_path = sequencing_download_path)
        else:
            print("Previous sample sheet {} detected. Running sample file generation and checking for duplicate IDs.".format(args.prevSamples))
            sample_sheet = generate_sample_file(sequencing_data_path = sequencing_download_path, prev_samples_path= args.prevSamples)
        sample_sheet.to_csv("/tmp/cdc_samples.csv")
    

