import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 
import re


def remove_zero_between_letter_and_digit(text):
    """
    Removes a '0' character that is immediately preceded by a letter
    and immediately followed by a digit. We use this when merging drugsids and sampleids into the experiments df.
    """
    # The regex pattern looks for:
    # (r'([a-zA-Z])0(\d)')
    # ([a-zA-Z]): Captures any single uppercase or lowercase letter (Group 1)
    # 0: Matches the literal '0'
    # (\d): Captures any single digit (Group 2)
    # The replacement string r'\1\2' puts Group 1 and Group 2 back together,
    # effectively removing the '0' in between.
    return re.sub(r'([a-zA-Z])0(\d)', r'\1\2', text)

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



### Parse Data Function
def parse_experiments_excel_sheets(first_file_path, second_file_path):
    # read in the excel files
    first_exp_excel = pd.ExcelFile(open(first_file_path, 'rb'))
    first_experiments_dict = pd.read_excel(first_exp_excel, sheet_name=None, header=None)
    rest_exp_excel = pd.ExcelFile(open(second_file_path, 'rb'))
    rest_experiments_dict = pd.read_excel(rest_exp_excel, sheet_name=None, header=None)
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
                # print("length before melt is:", len(one_conc_df))
                # print("end index is ",(conc_indexes[(index+1)]-1))
                one_conc_df.columns = one_conc_df.iloc[0]
                one_conc_df = one_conc_df[1:]
                one_conc_df = pd.melt(one_conc_df, id_vars=['concentration'], value_vars=one_conc_df.columns[one_conc_df.columns != 'concentration'])
                one_conc_df = one_conc_df.rename(columns={"concentration":"drug_id",one_conc_df.columns[1]:"DOSE","value":"count"})
                one_conc_df = one_conc_df.astype({"DOSE": 'float'})
                one_conc_df = one_conc_df.reset_index(drop=True)
                # now convert all counts to growth rates
                for drug in one_conc_df['drug_id'].unique():
                    # print("the drug name is",drug)
                    # print("the mean of the 0's is ",one_conc_df[(one_conc_df['drug_id'] == drug) & (one_conc_df['DOSE'] == 0)]['count'].mean())
                    mean_of_zeros = one_conc_df[(one_conc_df['drug_id'] == drug) & (one_conc_df['DOSE'] == 0)]['count'].mean()
                    one_conc_df.loc[one_conc_df['drug_id'] == drug, 'GROWTH'] = (one_conc_df[(one_conc_df['drug_id'] == drug)]['count']/mean_of_zeros)*100
                # print("sample is ", experiment_key)
                # print("index is ",index)
                # print("length of df is", len(one_conc_df))
                # print("number of unique drugs is", one_conc_df['drug_id'].nunique())
                list_of_dfs.append(one_conc_df)
            # print(list_of_dfs)
            elongated_df = pd.concat(list_of_dfs)
            elongated_df['sample_name'] = experiment_key
            # print(experiment_key)
            # print(elongated_df['drug_id'].nunique())
            list_of_finished_dfs.append(elongated_df)
        full_experiments_df = pd.concat(list_of_finished_dfs)
        full_df_list.append(full_experiments_df)
    experiments_df = pd.concat(full_df_list)
    return(experiments_df)
    

def merge_improve_samples_drugs(experiment_data:pd.DataFrame, samples_data_path:str, drugs_info_path:str, improve_drugs_path:str):
    # read in data
    experiments_df = experiment_data
    improve_sample_df = pd.read_csv(samples_data_path)
    improve_drug_df = pd.read_csv(improve_drugs_path, sep='\t')
    druginfo_df = pd.read_excel(drugs_info_path)
    # merging improve drug id's 
    drugnames_merged = pd.merge(experiments_df, druginfo_df[['Catalogue','Drug']], how = 'inner', left_on= "drug_id", right_on= "Catalogue")
    drugnames_merged['Drug'] = drugnames_merged['Drug'].str.lower()
    drugids_merged = pd.merge(drugnames_merged, improve_drug_df[['improve_drug_id','chem_name']], how = 'inner', left_on= "Drug", right_on= "chem_name")
    # merging improve sample id's 
    drugids_merged['sample_name'] = drugids_merged['sample_name'].apply(remove_zero_between_letter_and_digit) # need to apply this function bc some of the naming conventions for the sample names are inconsistent (ex: HCCO01 and HCCO1)
    all_merged = pd.merge(drugids_merged, improve_sample_df[['other_id','improve_sample_id']], how = 'left', left_on= "sample_name", right_on= "other_id")
    # now do some formatting
    all_merged = all_merged[all_merged['DOSE'] != 0] # get rid of any dose that is 0 bc that will cause issues during curve fitting
    all_merged['time'] = 72
    all_merged['time_unit'] = "hours"
    all_merged['study'] = "liver"
    all_merged['source'] = "synapse"
    all_merged = all_merged.drop(columns={'drug_id','count', 'sample_name','Catalogue','chem_name','other_id','Drug'})
    all_merged = all_merged.rename(columns={'improve_drug_id':'Drug'})
    
    # identify rows where improve_sample_id is NaN or non-finite
    all_merged['improve_sample_id'] = pd.to_numeric(all_merged['improve_sample_id'], errors='coerce')
    bad_mask = all_merged['improve_sample_id'].isna() | np.isinf(all_merged['improve_sample_id'])

    print(f"Rows before dropping bad improve_sample_id: {len(all_merged)}")
    if bad_mask.any():
        print(f"{bad_mask.sum()} rows with missing/non-finite improve_sample_id will be dropped")
        # drop and report after
        all_merged = all_merged.loc[~bad_mask].copy()
        print(f"Rows after dropping: {len(all_merged)}")

    # now safe to cast
    all_merged['improve_sample_id'] = all_merged['improve_sample_id'].astype(int)
    all_merged = all_merged[['study','time','DOSE','GROWTH','Drug','improve_sample_id','time_unit','source']]
    all_merged = all_merged.dropna() # drop na's bc that will also cause issues in curve fitting

    return(all_merged)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    # arguments for what data to process
    parser.add_argument('-D', '--Download', action = 'store_true', default=False, help='Download experiments data.')
    parser.add_argument('-t', '--Token', type=str, default=None, help='Synapse Token')
    parser.add_argument('-E', '--Experiment', action = 'store_true', default=False, help='Create experiments data.')
    parser.add_argument('-s', '--Samples', type=str, default=None, help='Path to samples file.')
    parser.add_argument('-d', '--Drugs', type=str, default=None, help='Path to drugs file')

    args = parser.parse_args()


    ###########################

    if args.Download:
        if args.Token is None:
            print("No synpase download tocken was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # download experiments data from synapse, which are split into 2 excel files
            first_experiments_path = download_experiments_data(synID="syn66401301", save_path="/tmp/", synToken=args.Token)
            rest_experiments_path = download_experiments_data(synID="syn66401302", save_path="/tmp/", synToken=args.Token)
    if args.Experiment:
        if args.Samples is None:
            print("No path to samples file detected. Cannot generate experiment data.")
            exit()
        if args.Drugs is None:
            print("No path to drugs file detected. Cannot generate experiment data.")
            exit()
        else:
            print("Parsing experiments excel sheets")
            parsed_experiments_data = parse_experiments_excel_sheets(first_experiments_path, rest_experiments_path)
            print("Generating experiments data.")
            experiments_df = merge_improve_samples_drugs(experiment_data = parsed_experiments_data, samples_data_path = args.Samples, improve_drugs_path = args.Drugs, drugs_info_path="/tmp/4_Drug_information.xlsx")
            output_path = "/tmp/liver_experiments_for_curve_fitting.tsv"
            print("Experiments data sucessfully generated.  Saving tsv to {}".format(output_path))
            experiments_df.to_csv(output_path, sep='\t')