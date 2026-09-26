


import os
import pandas as pd
import wget
import argparse
import synapseclient as sc
import math
import re


def get_precalc_auc():
    '''
    get pre-calculaterd AUC from supp data
    '''
    # Figshare copy of the Tiriac et al. supplementary drug table (supp_4775187).
    # The aacr.silverchair-cdn.com url previously here was a CloudFront SIGNED link
    # whose Expires lapsed on 2025-01-27, returning 403 permanently; a fresh one
    # expires about five weeks after issue. Figshare links do not expire.
    # Tiriac H, et al. Cancer Discovery (2018) 8(9):1112-1129.
    # doi:10.1158/2159-8290.CD-18-0349
    tablink = 'https://ndownloader.figstatic.com/files/39996295'

    chemo = pd.read_excel(tablink,sheet_name=1,skiprows=1)
    chemo.columns = [c.lower() for c in chemo.columns]
    targeted = res = pd.read_excel(tablink,sheet_name=2,skiprows=1)
    targeted.columns = [c.lower() for c in targeted.columns]

    cdat = chemo.melt(id_vars='sample id',value_vars=['gemcitabine','paclitaxel','sn-38','5-fu','oxaliplatin'],var_name='drug',value_name='published_auc')
    tdat = targeted.melt(id_vars='sample id',value_vars=set(targeted.columns)-set('sample id'),var_name='drug',value_name='published_auc')
    combined= pd.concat([cdat,tdat])
    combined = combined.rename(columns={'sample id':'other_id','drug':'chem_name'})
    
    return combined

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples', help='Sample mapping file for pancreatic samples')
    parser.add_argument('-d', '--drugs', help='Drug mapping file for pancreatic samples')
    parser.add_argument('-e', '--expfile', default = '/tmp/pancreatic_experiments.tsv',help='Output file to be read into curve fitting code')

    args = parser.parse_args()
    samples = pd.read_csv(args.samples,sep=',')
    drugs = pd.read_csv(args.drugs,sep='\t')

    newdat = get_precalc_auc().merge(samples).merge(drugs)
    newdat = newdat[['improve_sample_id','improve_drug_id','published_auc']].drop_duplicates()
    newdat = newdat.melt(id_vars=['improve_sample_id','improve_drug_id'],value_vars='published_auc',var_name='dose_response_metric',value_name='dose_response_value')
    newdat[['source']]='TiriacEtAl2018'
    newdat[['time']]=120
    newdat[['time_unit']]='hours'
    newdat[['study']]='pancreatic'
    olddat = pd.read_csv(args.expfile,sep='\t')
    res = pd.concat([olddat,newdat])
    res.to_csv(args.expfile)
