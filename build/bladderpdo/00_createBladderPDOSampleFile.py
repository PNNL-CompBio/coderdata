import synapseclient
import pandas as pd
import numpy as np
import argparse
import os


def get_bladder_pdo_samples(synLoginObject, maxval):
    
    # download from Synapse..
    samples_syn = synLoginObject.get('syn64765486')
    # and read the file
    samples_df = pd.read_csv(samples_syn.path, sep="\t")

    samples = samples_df[['Sample ID', 'Patient ID', 'Cancer Type Detailed', 'Sample Class']]
    samples = samples.rename({"Sample ID" : 'other_id', 'Patient ID' : 'common_name', 'Cancer Type Detailed': 'cancer_type', 'Sample Class' : 'model_type'}, axis=1)

    samples.loc[:,['species']] = 'Homo sapiens(Human)'
    samples.loc[:,['other_id_source']] = 'Synapse'
    samples.loc[:,['other_names'] ]= ''
    samples.loc[:,['cancer_type']]=samples['cancer_type'].str.lower()
    samples.loc[:, ['model_type']] = samples['model_type'].str.lower()

    samples['improve_sample_id'] = range(maxval+1, maxval+1+samples.shape[0])

    return samples


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of sample files for the Sarcoma PDO project into a single samplesheet")
    
    parser.add_argument('-t', '--token', type=str, help='Synapse Token')

    parser.add_argument("-p", '--prevSamples', nargs="?", type=str, default ="", const  = "", help = "Use this to provide previous sample file, will run sample file generation")

    args = parser.parse_args()
   
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)

    if (args.prevSamples):
        prev_max_improve_id = max(pd.read_csv(args.prevSamples).improve_sample_id)
    else: 
        prev_max_improve_id = 0

    bladder_pdo_samples = get_bladder_pdo_samples(synObject, prev_max_improve_id)

    bladder_pdo_samples.to_csv("bladderpdo_samples.csv", index=False)