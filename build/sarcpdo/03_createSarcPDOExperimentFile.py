import synapseclient
import pandas as pd
import numpy as np
import argparse
import os


if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of drug information for the Sarcoma PDO project into an experiments file")
    
    parser.add_argument('-t', '--token', type=str, help='Synapse Token')
    parser.add_argument('-s', '--samplesFile', nargs = "?", type=str, default = "", help = "Use this to provide previously generated sample file for this dataset to link to experiment data.")
    parser.add_argument('-d', '--drugFile', nargs = "?", type=str, default = "", help = "Use this to provide previously generated drugs file for this dataset to link with to experiment data.")

    args = parser.parse_args()
    print(args)
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)

    drug_query = synObject.tableQuery("select * from syn61892224")
    drug_data = drug_query.asDataFrame()

    # convert Drug_Name to lowercase for merge with drug info files
    drug_data['chem_name'] = drug_data['Drug_Name'].str.lower()

    sarcpdo_samples = pd.read_csv(args.samplesFile)

    sarcpdo_drugs = pd.read_csv(args.drugFile, sep="\t")
    # reformat 'other_id', specifically, alter underscores before '1' to dashes so we can split on "_" and retain appended numbers
    # rename to "Sample_ID" for merging with drug_data
    sarcpdo_samples['Sample_ID'] = sarcpdo_samples['other_id'].str.lower()
    sarcpdo_samples["Sample_ID"] = sarcpdo_samples["Sample_ID"].str.replace("_1", "-1")
    sarcpdo_samples['Sample_ID'] = sarcpdo_samples["Sample_ID"].str.split("_", expand =True)[0]
    # and change dashes back to underscores to merge with drug_data's Sample_ID
    sarcpdo_samples["Sample_ID"] = sarcpdo_samples["Sample_ID"].str.replace("-", "_")

    # inner merge with samples because there are samples without experiment info and many Sample_ID's in experiments data without sample info
    experiments = drug_data.merge(sarcpdo_drugs, how='left').merge(sarcpdo_samples, how='inner')

    final_experiment = experiments[['improve_sample_id', 'improve_drug_id', 'Viability_Score']]
    final_experiment.loc[:,['study']] = 'Landscape of Sarcoma'
    final_experiment.loc[:,['source']] = 'pharmacoGX'
    final_experiment.loc[:,['time']] = None
    final_experiment.loc[:,['time_unit']]= None
    final_experiment.loc[:,['dose_response_metric']] = 'published_auc' 
    final_experiment.loc[:,['dose_response_value']] = final_experiment['Viability_Score']

    toReturn = final_experiment[['source', 'improve_sample_id', 'improve_drug_id', 'study', 'time', 'time_unit', 'dose_response_metric', 'dose_response_value']]

    toReturn.to_csv('/tmp/sarcpdo_experiments.tsv', sep='\t', index=False)


    # to test run
    #  python3 03_createSarcPDOExperimentFile.py -t $SYNAPSE_AUTH_TOKEN -s sarcpdo_samples.csv -d sarcpdo_drugs.tsv