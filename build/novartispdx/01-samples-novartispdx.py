import pandas as pd
import synapseclient
import numpy as np
import argparse
import os

def get_complete_novartispdx_sample_sheet(synObject):

    files = list(synObject.getChildren(parent='syn66275995', includeTypes=['file']))

    synIDs = [item['id'] for item in files]
    # leave off synIDs for drug info
    synIDs.remove('syn66276102')
    synIDs.remove('syn66276098')
    synIDs.remove("syn66477971")
    # create empty dataframe
    allsamplesheet = pd.DataFrame()
    # iterate through IDs and concatenate
    for id in synIDs:
        curr = synObject.get(id)
        currdf = pd.read_csv(curr.path)
        allsamplesheet = pd.concat([allsamplesheet, currdf], ignore_index=True)
    # rename columns and reformat cancer type from CANCER_HISTOLOGY column
    allsamplesheet['other_id'] = allsamplesheet['Sample ID']
    allsamplesheet['common_name'] = allsamplesheet['MODEL_ORIGINATOR_ID']
    allsamplesheet['cancer_type'] = allsamplesheet['CANCER_HISTOLOGY'].str.lower().str.split(pat="^[^\s]*\s", expand=True)[1]
    allsamplesheet['species'] = "Homo Sapiens(human)"
    allsamplesheet['model_type'] = 'patient derived xenograft'
    allsamplesheet['other_id_source'] = 'Synapse'
    allsamplesheet['other_names'] = ''
    finalsamplesheet = allsamplesheet[['other_id', 'common_name', 'other_id_source', 'other_names', 'cancer_type', 'species', 'model_type']]
    return finalsamplesheet

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of sample files for the Novartis PDX data into a single samplesheet")
    
    parser.add_argument('-t', '--token', type=str, help='Synapse Token')

    parser.add_argument("-p", '--prevSamples', nargs="?", type=str, default ="", const  = "", help = "Use this to provide previous sample file, will run sample file generation")

    args = parser.parse_args()
   
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)
    samplesheet = get_complete_novartispdx_sample_sheet(synObject)

    if (args.prevSamples):
        prev_max_improve_id = max(pd.read_csv(args.prevSamples).improve_sample_id)
    else: 
        prev_max_improve_id = 0

    samplesheet['improve_sample_id'] = range(prev_max_improve_id+1, prev_max_improve_id+samplesheet.shape[0]+1) 

    samplesheet.to_csv('/tmp/novartispdx_samples.csv', index=False)

        