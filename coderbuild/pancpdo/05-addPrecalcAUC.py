


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
    tablink = 'https://aacr.silverchair-cdn.com/aacr/content_public/journal/cancerdiscovery/8/9/10.1158_2159-8290.cd-18-0349/5/21598290cd180349-sup-199398_2_supp_4775187_p95dln.xlsx?Expires=1738004990&Signature=av8XadTm9AmI20O2Y7J7aHDtPbpluKJIfI5ubsoiYJ15D0zh5p1ltF4a7-DCSWTSMs-qX5TD09shxHeqkQ2NkLWHZsXoCD5KyREGhEgcDAvWZ1V9kwXDm0bjpINipAPPtC20oeuw6c~hPooF3Mtgzp4MzMCCjcVwfn05u27a0kS0yifBi11wQj3nmHlR3ym-2fYkFuqQtnNPCzH8-yIw21y0kTvXrNodAzC5pGA8qUK4PLxBt52xUIvTEPsPiPjXwBnDCfVsLGGdDYIY25lEPKiA403q6kFYvrSQ3bsTvM4kuvltb7yS4AXjK0-tthMOKbqq8~uREmJCcueADUF91g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'

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
    parser.add_argument('-s', '--samples', help='Sample mapping file for panc pdo samples')
    parser.add_argument('-d', '--drugs', help='Drug mapping file for panc pdo samples')
    parser.add_argument('-e', '--expfile', default = '/tmp/pancpdo_experiments.tsv',help='Output file to be read into curve fitting code')

    args = parser.parse_args()
    samples = pd.read_csv(args.samples,sep=',')
    drugs = pd.read_tsv(args.drugs,sep='\t')

    newdat = get_precalc_auc().merge(samples).merge(drugs)
    newdat = newdat[['improve_sample_id','improve_drug_id','published_auc']].drop_duplicates()
    newdat = newdat.melt(id_vars=['improve_sample_id','improve_drug_id'],value_vars='published_auc',var_name='dose_response_metric',value_name='dose_response_value')
    newdat[['source']]='TiriacEtAl2018'
    newdat[['time']]=120
    newdat[['time_unit']]='hours'
    newdat[['study']]='pancpdo'
    oldat = pd.read_csv(args.expfile,sep='\t')
    res = pd.concat([olddat,newdat])
    res.to_csv(args.expfile)
