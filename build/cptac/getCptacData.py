#!/usr/bin/env python3
'''
Script to import CPTAC proteogenomic data and handle previous samples argument.
'''

import argparse
import cptac
import os
import numpy as np
import pandas as pd
import gzip
from pathvalidate.argparse import validate_filename_arg

def reclass_var(arr):
    '''
    Reclassify variant annotations to a standard set.
    '''
    variant_schema = {
        "3'UTR": ["3'UTR"],
        "5'Flank": ["5'Flank"],
        "5'UTR": ["5'UTR"],
        'Undetermined': ['COULD_NOT_DETERMINE'],
        'De_novo_Start_InFrame': ['DE_NOVO_START_IN_FRAME', 'De_novo_Start_InFrame'],
        'De_novo_Start_OutOfFrame': ['DE_NOVO_START_OUT_FRAME', 'De_novo_Start_OutOfFrame'],
        'Frame_Shift_Del': ['Frame_Shift_Del'],
        'Frame_Shift_Ins': ['Frame_Shift_Ins'],
        'IGR': ['IGR'],
        'In_Frame_Del': ['In_Frame_Del'],
        'In_Frame_Ins': ['In_Frame_Ins'],
        'Intron': ['Intron'],
        'Missense_Mutation': ['Missense_Mutation'],
        'Nonsense_Mutation': ['Nonsense_Mutation'],
        'Nonstop_Mutation': ['Nonstop_Mutation'],
        'RNA': ['RNA'],
        'Start_Codon_SNP': ['START_CODON_SNP', 'Start_Codon_SNP'],
        'Start_Codon_Del': ['Start_Codon_Del'],
        'Start_Codon_Ins': ['Start_Codon_Ins', 'START_CODON_INS'],
        'Stop_Codon_Del': ['Stop_Codon_Del'],
        'Stop_Codon_Ins': ['Stop_Codon_Ins'],
        'Silent': ['Silent'],
        'Splice_Site': ['Splice_Site'],
        'Translation_Start_Site': ['Translation_Start_Site']
    }

    newarr = []
    for a in arr:
        found = False
        for key, var_list in variant_schema.items():
            if a in var_list:
                newarr.append(key)
                found = True
                break
        if not found:
            print(f"Variant {a} not found, returning 'Undetermined'")
            newarr.append('Undetermined')
    print(f'Updated {len(arr)} old variants with {len(newarr)} new ones')
    return newarr

def getCancerObj(cancertype):
    '''
    Retrieve the CPTAC data object for a given cancer type.
    '''
    if cancertype == 'brca':
        dat = cptac.Brca()
    elif cancertype == 'ccrcc':
        dat = cptac.Ccrcc()
    elif cancertype == 'coad':
        dat = cptac.Coad()
    elif cancertype == 'gbm':
        dat = cptac.Gbm()
    elif cancertype == 'hnscc':
        dat = cptac.Hnscc()
    elif cancertype == 'lscc':
        dat = cptac.Lscc()
    elif cancertype == 'luad':
        dat = cptac.Luad()
    elif cancertype == 'ov':
        dat = cptac.Ov()
    elif cancertype == 'pdac':
        dat = cptac.Pdac()
    elif cancertype == 'ucec':
        dat = cptac.Ucec()
    else:
        print('Wrong cancer type: ' + cancertype)
        exit()
    return dat

ct_vals = {
    'brca': 'Breast carcinoma',
    'ccrcc': 'Clear cell renal cell carcinoma',
    'coad': 'Colon adenocarcinoma',
    'gbm': 'Glioblastoma multiforme',
    'hnscc': 'Head and Neck squamous cell carcinoma',
    'lscc': 'Lung squamous cell carcinoma',
    'luad': 'Lung adenocarcinoma',
    'ov': 'Ovarian carcinoma',
    'pdac': 'Pancreatic ductal adenocarcinoma',
    'ucec': 'Uterine Corpus Endometrial Carcinoma'
}

def buildTumorSampleTable(sample_names, cancer_type, samples, maxval):
    '''
    Builds or updates the samples DataFrame with new samples.

    Parameters:
    - sample_names: list of sample names to process.
    - cancer_type: string indicating the cancer type.
    - samples: existing samples DataFrame.
    - maxval: integer indicating the starting improve_sample_id.

    Returns:
    - Updated samples DataFrame.
    - Updated maxval.
    '''
    for samp in sample_names:
        if len(samples) == 0 or samp not in samples.common_name.values:
            maxval = int(maxval + 1)
            new_sample = pd.DataFrame({
                'common_name': [samp],
                'cancer_type': [ct_vals[cancer_type]],
                'other_names': [''],
                'species': ['Homo sapiens (Human)'],
                'improve_sample_id': [maxval],
                'other_id_source': ['CPTAC3'],
                'other_id': [samp],
                'model_type': ['tumor']
            })
            samples = pd.concat([samples, new_sample], ignore_index=True)
    samples = samples.reset_index(drop=True)
    return samples, maxval

def formatMutData(df, dtype, ctype, samp_names, source, samples):
    '''
    Formats mutational data.
    '''
    # Some datasets have two levels of indices, some have only one
    if df.columns.nlevels == 2:
        df.columns = df.columns.droplevel(1)
    elif df.columns.nlevels != 1:
        print("The number of column levels is not 1 or 2!\n")
        raise Exception("Invalid column levels")

    df = df.reset_index()
    subset = df[['Patient_ID', 'HGNC_Entrez_Gene_ID(supplied_by_NCBI)', 'Mutation', 'Genome_Change']]
    subset.columns = ['Patient_ID', 'entrez_gene', 'variant_classification', 'Mutation']

    improve_mapping = samples.set_index('common_name')
    blongdf = subset.join(improve_mapping, on='Patient_ID', how='inner')
    blongdf['source'] = source
    blongdf['study'] = 'CPTAC3'

    if len(blongdf) > 0:
        blongdf['variant_classification'] = reclass_var(blongdf['variant_classification'])
    else:
        blongdf['variant_classification'] = blongdf['variant_classification']

    blongdf = blongdf.rename(columns={
        'entrez_gene': 'entrez_id',
        'Mutation': 'mutation'
    })
    blongdf = blongdf[['improve_sample_id', 'entrez_id', 'mutation', 'variant_classification', 'source', 'study']]
    return blongdf


def formatData(df,dtype,ctype,samp_names,source,genes,samples):
    '''
    formats data into long form
    '''
    # some dataset has two level of indices some has only one
    if df.columns.nlevels == 2:
        df.columns = df.columns.droplevel(1)
    elif df.columns.nlevels != 1:
        print("The number of column levels is larger not 1 or 2!\n")
        raise

    gene_names = list(df.columns)
    ##pivot wide to long format
    df = df.reset_index()

    longdf = pd.melt(df,id_vars='Patient_ID',value_vars=gene_names,var_name='gene_symbol',value_name=dtype)
    ##match sample identifiers or build new ones

    #     snames=list(set(longdf.Patient_ID))
    improve_mapping = samples
    improve_mapping.reset_index(drop=True)
    improve_mapping.index=improve_mapping.common_name

    blongdf = longdf.join(improve_mapping,on='Patient_ID')
    improve_mapping.reset_index(drop=True)
    ##match gene identifiers

    genes.index=genes.gene_symbol
    mlongdf = blongdf.join(genes,on='gene_symbol',how='left',lsuffix='.e',rsuffix='.m') ##fix this
    
    mlongdf = mlongdf[['entrez_id','improve_sample_id',dtype]].drop_duplicates()

    ##if its copy number we need to include copy number call
    if dtype=='copy_number': ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
       mlongdf[['copy_call']] = mlongdf[['copy_number']].apply(copy_num)


    mlongdf[['source']] = source
    mlongdf[['study']] = 'CPTAC3'

    return mlongdf


def copy_num(arr):
    '''
    Converts copy number values to categorical calls.
    '''
    copy_call=[]
    for a in arr:
        a = 2**float(a)
        if float(a) < 0.5210507:
            b = 'deep del'
        elif float(a) < 0.7311832:
            b = 'het loss'
        elif float(a) < 1.214125:
            b = 'diploid'
        elif float(a) <1.42233:
            b = 'gain'
        else:
            b = 'amp'
        copy_call.append(b)
    return copy_call



def main():
    '''
    Main function that processes CPTAC data and handles previous samples.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--prevSampleFile', dest='sampfile', nargs='?', type=str, default=None, const=None,
                        help='Sample file to use to generate new ids. Returns sample file')
    parser.add_argument('--geneFile', dest='genefile', default='/tmp/genes.csv',
                        help='Gene file to get gene ids')
    parser.add_argument('--curSampleFile', dest='newsamps', default=None,
                        help='Sample file to use to generate data. Returns data for samples')
    opts = parser.parse_args()

    dat_files = {}

    # Decide whether to build samples, data, or both
    build_data = False
    build_samples = False
    maxval = 1  # Initialize maxval to 1

    if opts.sampfile is not None:
        if opts.sampfile != '':
            # Previous samples file provided and not empty
            print("Previous Samples File Provided: ", opts.sampfile)
            if os.path.exists(opts.sampfile):
                old_samples = pd.read_csv(opts.sampfile)
                samples = old_samples.copy()
                maxval = samples['improve_sample_id'].max()
                print("Max improve_sample_id from previous samples: ", maxval)
            else:
                print("Previous samples file not found. Exiting.")
                exit()
            build_samples = True
            print("Generating sample file")
        else:
            # opts.sampfile is an empty string
            print("No previous samples file provided. Proceeding with maxval=1")
            samples = pd.DataFrame(columns=[
                'common_name', 'cancer_type', 'other_names', 'species',
                'improve_sample_id', 'other_id_source', 'other_id', 'model_type'
            ])
            maxval = 1
            build_samples = True
            print("Generating sample file")
    elif opts.newsamps is not None and opts.newsamps != '':
        # Build data
        build_data = True
        print("Current Samples File Provided: ", opts.newsamps)
        if os.path.exists(opts.newsamps):
            samples = pd.read_csv(opts.newsamps)
            # Read in existing gene identifier table
            genes = pd.read_csv(opts.genefile)
            print("Building dataset from generated sample file")
        else:
            print("Current samples file not found. Exiting.")
            exit()
    else:
        print("No sample file provided. Exiting.")
        exit()

    # Proceed to build samples and/or data
    if build_samples:
        # Loop through the cancer types to build samples
        for cancertype in ['brca', 'coad', 'hnscc', 'lscc', 'luad', 'ov', 'gbm', 'pdac', 'ucec', 'ccrcc']:
            dat = getCancerObj(cancertype)

            # Get available data for cancer
            dat_list = dat.list_data_sources().set_index('Data type').to_dict()['Available sources']
            clinsource = dat_list.get('clinical', [])
            if 'harmonized' in clinsource:
                cs = 'harmonized'
            else:
                cs = clinsource[0] if clinsource else None
            if cs is None:
                print(f"No clinical data available for {cancertype}. Skipping.")
                continue

            tumor_samps = dat.get_clinical(cs)
            tumor_samps = list(tumor_samps.index)

            print('Building sample file for ' + cancertype)
            samples, maxval = buildTumorSampleTable(tumor_samps, cancertype, samples, maxval)

        # Save the updated samples file
        samples = samples.reset_index(drop=True)
        samples.to_csv('/tmp/cptac_samples.csv', index=False)
        print("Samples file saved to /tmp/cptac_samples.csv")

    if build_data:
        # Now proceed to build data
        # Ensure that 'samples' and 'genes' are loaded
        # Loop through the cancer types to build data
        for cancertype in ['brca', 'coad', 'hnscc', 'lscc', 'luad', 'ov', 'gbm', 'pdac', 'ucec', 'ccrcc']:
            dat = getCancerObj(cancertype)

            # Get available data for cancer
            dat_list = dat.list_data_sources().set_index('Data type').to_dict()['Available sources']
            clinsource = dat_list.get('clinical', [])
            if 'harmonized' in clinsource:
                cs = 'harmonized'
            else:
                cs = clinsource[0] if clinsource else None
            if cs is None:
                print(f"No clinical data available for {cancertype}. Skipping.")
                continue

            tumor_samps = dat.get_clinical(cs)
            tumor_samps = list(tumor_samps.index)

            print('Building data files for ' + cancertype)
            all_dfs = {}
            all_sources = {}
            for dtype in ['somatic_mutation', 'proteomics', 'transcriptomics', 'CNV']:
                if dtype not in dat_list.keys():
                    continue
                source_list = dat_list[dtype]
                if 'harmonized' in source_list:
                    source = 'harmonized'
                else:
                    source = source_list[0]
                all_sources[dtype] = source
                if dtype == 'proteomics':
                    all_dfs[dtype] = dat.get_proteomics(source)
                elif dtype == 'transcriptomics':
                    all_dfs[dtype] = dat.get_transcriptomics(source)
                elif dtype == 'somatic_mutation':
                    all_dfs[dtype] = dat.get_somatic_mutation(source)
                elif dtype == 'CNV':
                    all_dfs[dtype] = dat.get_CNV(source)

            for dtype, df in all_dfs.items():
                df = df.loc[[t for t in tumor_samps if t in df.index]]
                df.dropna(how='all', axis=0, inplace=True)
                print(cancertype + ' ' + dtype)
                if dtype == 'somatic_mutation':
                    fdf = formatMutData(df, 'mutation', cancertype, tumor_samps, all_sources[dtype], samples)
                    fdf = fdf.reset_index(drop=True)
                    dtype_key = 'mutations'
                elif dtype == 'CNV':
                    fdf = formatData(df, 'copy_number', cancertype, tumor_samps, all_sources[dtype], genes, samples)
                    fdf = fdf.reset_index(drop=True)
                    dtype_key = 'copy_number'
                else:
                    fdf = formatData(df, dtype, cancertype, tumor_samps, all_sources[dtype], genes, samples)
                    fdf = fdf.reset_index(drop=True)
                    dtype_key = dtype

                if dtype_key in dat_files.keys():
                    of = dat_files[dtype_key].dropna()
                    fdf2 = pd.concat([of, fdf])
                    dat_files[dtype_key] = fdf2
                else:
                    dat_files[dtype_key] = fdf.dropna()
                
                print(dtype_key)

        # Now concatenate all the cancers into a single file
        for dtype_key, df in dat_files.items():
            print('Saving ' + "cptac_" + dtype_key + '.csv.gz' + ' file')
            print(df.to_string())
            df['entrez_id'] = df['entrez_id'].fillna(0)
            df['entrez_id'] = df['entrez_id'].astype(int)
            df.to_csv("/tmp/" + "cptac_" + dtype_key + '.csv.gz', sep=',', index=False, compression='gzip')

if __name__ == '__main__':
    main()