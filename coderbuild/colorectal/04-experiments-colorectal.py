import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 

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

def create_experiments_data(experiment_data_path:str, samples_data_path:str, drugs_data_path:str):
    raw_experiment_data = pd.read_csv(experiment_data_path)
    samples_data = pd.read_csv(samples_data_path)
    drug_data = pd.read_csv(drugs_data_path, sep="\t")
    # create experiments df with these columns : 
    # DOSE: dose of drug in uM,\
    # GROWTH: percentage of cells left,\
    # study: name of study to group measurements by,\
    # source: source of the data,\
    # improve_sample_id: improve_sample_id,\
    # Drug: improve_drug_id,\
    # time: time at which measurement was taken,\
    # time_unit: unit of time')
    experiment_data = raw_experiment_data[['RESEARCH_PROJECT' ,'CELL_LINE_NAME', 'DURATION', 'CONC', 'DRUG_NAME', 'viability']]
    # remove rows where there's no drug name
    experiment_data = experiment_data[experiment_data['DRUG_NAME'].notna()]
    # merge drugs
    experiment_data['DRUG_NAME'] = experiment_data['DRUG_NAME'].str.replace("Nutlin-3a (-)","Nutlin-3a")
    experiment_data['DRUG_NAME_LOWER'] = experiment_data['DRUG_NAME'].str.lower() # make everything lower bc that's how it is in the drug data
    drug_experiment_merge = pd.merge(experiment_data, drug_data[['improve_drug_id','chem_name']], how='left', left_on='DRUG_NAME_LOWER', right_on='chem_name')
    # merge samples 
    # need to change CELL_LINE_NAME to just get patient number  (they're all tumor organoid)
    drug_experiment_merge['patient_number'] = drug_experiment_merge['CELL_LINE_NAME'].str.split("T",expand=True).iloc[:,1]+drug_experiment_merge['CELL_LINE_NAME'].str.split("T",expand=True).iloc[:,2]
    # get samples to only tumor organoid
    tumor_org_samples = samples_data[samples_data['other_id'].str.contains("Tumor-Organoid")]
    tumor_org_samples['patient_number'] = samples_data['other_id'].str.split("-",expand=True).iloc[:,0].str.replace("P","").str.replace("T","")
    sample_drug_experiment_merge = pd.merge(drug_experiment_merge,tumor_org_samples[['patient_number','improve_sample_id']], how='left', on='patient_number')

    # clean up table by dropping and renaming columns
    sample_drug_experiment_merge = sample_drug_experiment_merge.rename(columns = {'CONC':'DOSE','viability':'GROWTH','improve_drug_id':'Drug','DURATION':'time','RESEARCH_PROJECT':'study'})
    sample_drug_experiment_merge['time_unit'] = "days"
    sample_drug_experiment_merge['source'] = "vandeWetering_2015"
    sample_drug_experiment_merge = sample_drug_experiment_merge.drop(columns=['DRUG_NAME_LOWER','CELL_LINE_NAME','DRUG_NAME','chem_name','patient_number'])
    return(sample_drug_experiment_merge)




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
            # download experiments data from synapse
            experiments_data_path = download_synapse_data(synID = "syn65452842", save_path = "/tmp/", synToken = args.Token)
    if args.Experiment:
        if args.Samples is None:
            print("No path to samples file detected. Cannot generate experiment data.")
            exit()
        if args.Drugs is None:
            print("No path to drugs file detected. Cannot generate experiment data.")
            exit()
        else:
            print("Generating experiments data.")
            experiments_df = create_experiments_data(experiment_data_path = "/tmp/raw_data_GDSC_Org_restricted_11Mar25_plus_viabilities.csv", samples_data_path = args.Samples, drugs_data_path = args.Drugs)
            output_path = "/tmp/colorectal_experiments_for_curve_fitting.tsv"
            print("Experiments data sucessfully generated.  Saving tsv to {}".format(output_path))
            experiments_df.to_csv(output_path, sep='\t')
