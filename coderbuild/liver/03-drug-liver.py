import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 
from pubchem_retrieval import update_dataframe_and_write_tsv
import warnings
warnings.filterwarnings("ignore")


# function for loading data
def download_parse_drug_data(synID:str , save_path:str = None, synToken:str = None):
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
    mutations_data : pd.Dataframe
        A pandas dataframe containing mutations data

    copynum_data : pd.Dataframe
        A pandas dataframe containing copy number data

    proteomics_data : pd.Dataframe
        A pandas dataframe containing proteomics data
    """
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    downloaded_data = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    drugs_filepath = downloaded_data.path
    
    return(drugs_filepath)


# def create_liverpdo_drug_data(drug_info_path:str, prevDrugFilepath:str, output_drug_data_path:str):
#     # import fitted drug data and get drug names from DRUG_NAME column
#     drug_info_df = pd.read_csv(drug_info_path)
#     liverpdo_drugs_df = pd.DataFrame({"chem_name":drug_info_df['Drug'].unique()})
#     # if there is a prev drug file, check for new drugs
#     if prevDrugFilepath != "":
#         if prevDrugFilepath.__contains__(".tsv"):
#             prev_drug_df = pd.read_csv(prevDrugFilepath, sep='\t')
#         else:
#             prev_drug_df = pd.read_csv(prevDrugFilepath)
#         # get drugs that are only in the crcpdo_drugs_df (aka new drugs only)
#         new_drugs_df = liverpdo_drugs_df[~liverpdo_drugs_df.chem_name.isin(prev_drug_df.chem_name)]
#     else:
#         # if there's no prev drugs, then all drugs are new
#         new_drugs_df = liverpdo_drugs_df
#     # get new drug names
#     new_drug_names = new_drugs_df['chem_name'].unique()
#     # call function that gets info for these drugs
#     update_dataframe_and_write_tsv(unique_names = new_drug_names,output_filename = output_drug_data_path)


def create_liver_drug_data(drug_info_path: str, prevDrugFilepath: str, output_drug_data_path: str):
    # read current liver drug names
    drug_info_df = pd.read_csv(drug_info_path)
    raw_drug_names = [str(x) for x in pd.Series(drug_info_df["Drug"].unique()) if pd.notna(x)]

    # delegate to centralized PubChem retrieval logic, restricting output to only liver drugs
    update_dataframe_and_write_tsv(
        unique_names=raw_drug_names,
        output_filename=output_drug_data_path,
        prev_drug_filepaths=prevDrugFilepath if prevDrugFilepath else None,
        isname=True,
        batch_size=50,
        restrict_to_raw_names=raw_drug_names
    )
    
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
            print("No synpase download token was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # download fitted and raw drug data from synapse
            fitted_drug_data_path = download_parse_drug_data(synID = "syn66401300", save_path = "/tmp/", synToken = args.Token)
            drug_excel = pd.ExcelFile(open(fitted_drug_data_path, 'rb'))
            druginfo_df = pd.read_excel(drug_excel)
            druginfo_df.to_csv("/tmp/raw_druginfo.csv")
    if args.Drug:
        if args.PrevDrugs is None or args.PrevDrugs=='':
            print("No previous drugs file provided.  Starting improve_drug_id from SMI_1. Running drug file generation")
            create_liver_drug_data(drug_info_path = "/tmp/raw_druginfo.csv", output_drug_data_path = "/tmp/liver_drugs.tsv", prevDrugFilepath = "")
        else:
            print("Previous drugs file {} detected. Running drugs file generation and checking for duplicate IDs.".format(args.PrevDrugs))
            create_liver_drug_data(drug_info_path = "/tmp/raw_druginfo.csv", prevDrugFilepath = args.PrevDrugs, output_drug_data_path = "/tmp/liver_drugs.tsv")

