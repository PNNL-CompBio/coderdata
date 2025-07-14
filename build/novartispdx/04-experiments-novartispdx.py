import synapseclient
import pandas as pd
import numpy as np
import argparse
import os


# add improve IDs - for sample and drug
def get_novartis_pdx_experiments_file(synObject, samples_df): 
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
    groupeddf = novartispdx_curvefile.groupby('experiment')
    no_control = groupeddf['treatment'].apply(lambda x: x.str.contains('control').any())

    missingcontrols = no_control.reset_index()[no_control.reset_index()['treatment'] ==False]['experiment']
    nomissingcontrols=novartispdx_curvefile[~novartispdx_curvefile['experiment'].isin(missingcontrols)]
    #merge on drug names done in calc_pdx_metrics.py
    #final_w_drugIDs = finaldf.merge(drug_df, how='left',right_on='chem_name', left_on="treatment")
    final_allIDs = nomissingcontrols.merge(samples_df, how='left', right_on='common_name', left_on='model_id') 
    print(final_allIDs.head)
    final_allIDs = final_allIDs.drop('model_id', axis=1)
    finalDF = final_allIDs.rename({'improve_sample_id':'model_id'}, axis=1)
    print(finalDF.head)
    finalcurvefile = finalDF[['model_id', 'time', 'volume', 'treatment', 'experiment', 'dose']]
    print(finalcurvefile.head)
    return finalcurvefile


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--token', help='Synapse authentication token')
    parser.add_argument('-s', '--curSampleFile', default='/tmp/novartispdx_samples.csv', help='Sample mapping file for bladder pdo samples')
    parser.add_argument('-d', '--drugfile', default='/tmp/novartispdx_drugs.tsv', help='Drug mapping file for bladder pdo samples')
    parser.add_argument('-o', '--output', default = '/tmp/novartispdx_experiments.tsv',help='Output experiments file')

    args = parser.parse_args()
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)
    #drug_df = pd.read_csv(args.drugfile, sep='\t')
    samples_df = pd.read_csv(args.curSampleFile)
    
    doseresponse_data = get_novartis_pdx_experiments_file(synObject, samples_df)
    print(doseresponse_data.head)
    doseresponse_data.to_csv('/tmp/novartispdx_curvedata.tsv', columns=list({'model_id', 'time', 'volume', 'treatment','experiment', 'dose'}), sep='\t')

