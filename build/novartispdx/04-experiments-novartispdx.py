import synapseclient
import pandas as pd
import numpy as np
import argparse
import os

def get_novartis_pdx_experiments_file: 
    # input for the calc_pdx_metrics script

    file1 = synObject.get('syn66276102')
    rawDrugData = pd.read_csv(file1.path)
    # STILL NEED TO : link to improve ids. 
    # update a few drug ids for greater inclusion
    novartispdx_curvefile = rawDrugData[['Model', 'Days Post T0', 'Volume (mm3)', 'Treatment']]
    novartispdx_curvefile=novartispdx_curvefile.rename({'Model': 'model_id', 'Days Post T0' : 'time', 'Volume (mm3)': 'volume', 'Treatment':'treatment'}, axis=1)
    novartispdx_curvefile['treatment'] = novartispdx_curvefile['treatment'].str.lower()
    novartispdx_curvefile['treatment'] = novartispdx_curvefile['treatment'].str.replace('"', '')
    novartispdx_curvefile['treatment']=novartispdx_curvefile['treatment'].str.replace('untreated', 'control')
    novartispdx_curvefile['experiment'] = novartispdx_curvefile.groupby(['model_id']).ngroup()+1
    # remove triple combination(s)
    novartispdx_curvefile = novartispdx_curvefile[~novartispdx_curvefile['treatment'].str.contains(r'\+.*\+')]
    # remove dose information appended to some drugs in the treatment column and include in dose colum
    druganddose = novartispdx_curvefile['treatment'].str.split('-', expand=True)
    druganddose = druganddose.rename({0: 'treatment', 1:'dose'}, axis=1)
    novartispdx_curvefile['treatment']=druganddose['treatment']
    novartispdx_curvefile['dose'] = druganddose['dose']
    # remove pdxs with only one drug treatment (no control)
    unique_vals_tally = novartispdx_curvefile.groupby('experiment').nunique() 
    todiscard = unique_vals_tally[unique_vals_tally['treatment']==1].index
    novartispdx_curvefile = novartispdx_curvefile[~novartispdx_curvefile['experiment'].isin(todiscard)]
    # remove groups with no 'control' treatment
    groupeddf = test.groupby('experiment')
    no_control = groupeddf['treatment'].apply(lambda x: x.str.contains('control').any())

    missingcontrols = no_control.reset_index()[no_control.reset_index()['treatment'] ==False]['experiment']
    finaldf=test[~test['experiment'].isin(missingcontrols)]

    finalcurvefile = finaldf
    return finalcurvefile
    #finalcurvefile.to_csv('/tmp/novartispdx_doserep.tsv', sep="\t")        


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--token', help='Synapse authentication token')
    parser.add_argument('-s', '--curSampleFile', help='Sample mapping file for bladder pdo samples')
    parser.add_argument('-d', '--drugfile', help='Drug mapping file for bladder pdo samples')
    parser.add_argument('-o', '--output', default = '/tmp/novartispdx_doserep.tsv',help='Output file to be read into curve fitting code')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)
    drug_df = pd.read_csv(args.drugfile, sep='\t')
    samples_df = pd.read_csv(args.curSampleFile)

    doseresponse_data = get_novartis_pdx_experiments_file(synObject, samples_df, drug_df)
    doseresponse_data.to_csv(args.output, sep='\t')

