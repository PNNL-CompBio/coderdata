#!/usr/local/bin/python
'''
Basic CLI to import CPTAC protegenomic data
'''
import argparse
import cptac
import os
import numpy as np
import pandas as pd
import gzip
from pathvalidate.argparse import validate_filename_arg


def reclass_var(arr):
    # the dictionary started with
    # CCLE data
    newarr=[]
    variant_schema = {"3'UTR":["3'UTR"],\
                      "5'Flank":["5'Flank"],\
                      "5'UTR": ["5'UTR"], \
                      'Undetermined':['COULD_NOT_DETERMINE'],\
                      'De_novo_Start_InFrame':['DE_NOVO_START_IN_FRAME','De_novo_Start_InFrame'],\
                      'De_novo_Start_OutOfFrame':['DE_NOVO_START_OUT_FRAME','De_novo_Start_OutOfFrame'],\
                      'Frame_Shift_Del':['Frame_Shift_Del'],\
                      'Frame_Shift_Ins':['Frame_Shift_Ins'],\
                      'IGR':['IGR'],\
                      'In_Frame_Del':['In_Frame_Del'],\
                      'In_Frame_Ins':['In_Frame_Ins'],\
                      'Intron':['Intron'],\
                      'Missense_Mutation': ['Missense_Mutation'],\
                      'Nonsense_Mutation': ['Nonsense_Mutation'],\
                      'Nonstop_Mutation': ['Nonstop_Mutation'],\
                      'RNA':['RNA'],\
                      'Start_Codon_SNP':['START_CODON_SNP','Start_Codon_SNP'],\
                      'Start_Codon_Del':['Start_Codon_Del'],\
                      'Start_Codon_Ins':['Start_Codon_Ins','START_CODON_INS'],\
                      'Stop_Codon_Del' :['Stop_Codon_Del'],\
                      'Stop_Codon_Ins' :['Stop_Codon_Ins'],\
                      'Silent':['Silent'],\
                      'Splice_Site':['Splice_Site'],\
                      'Translation_Start_Site':['Translation_Start_Site']}
    #there has to be a better way to loop through each element of the dictionary
    #but here we are
    for a in arr:
        found = False
        for key,var in variant_schema.items():
            if a in var:
                newarr.append(key)
                found = True
        if not found:
            print("Variant "+a+' not found, returning undetermined')
            newarr.append('Undetermined')
    ##double check here
    print('udpdated '+str(len(arr))+' old variants with '+str(len(newarr))+' new ones')
    return newarr

def getCancerObj(cancertype):
   # cptac.download(dataset=cancertype,source='harmonized',)
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
    elif cancertype =='pdac':
        dat = cptac.Pdac()
    elif cancertype =='ucec':
        dat = cptac.Ucec()
    else:
        print('Wrong cancer type: '+cancertype)
        exit()
    return dat


ct_vals={'brca':'Breast carcinoma','ccrcc':'Clear cell renal cell carcinoma',\
         'coad':'Colon adenocarcinoma',\
         'gbm':'Glioblastoma multiforme',\
         'hnscc':'Head and Neck squamous cell carcinoma',\
         'lscc':'Lung squamous cell carcinoma','luad':'Lung adenocarcinoma',\
         'ov':'Ovarian carcinoma',\
         'pdac': 'Pancreatic ductal adenocarcinoma','ucec':'Uterine Corpus Endometrial Carcinoma'}


def buildTumorSampleTable(sample_names,cancer_type,samples,maxval):
    '''
    we have to read in existing sample information first to see if it is part of the samples
    '''
    for samp in sample_names:
       if len(samples)==0 or samp not in samples.common_name:
           maxval = int(maxval+1)
           samples = pd.concat([samples,
                                pd.DataFrame({'common_name':[samp],\
                                              'cancer_type':[ct_vals[cancer_type]],\
                                              'other_names':[''],\
                                              'species':['Homo Sapiens (Human)'],\
                                              'improve_sample_id':[maxval],\
                                              'other_id_source':['CPTAC3'],\
                                              'other_id':[samp],\
                                              'model_type':['tumor']}).reset_index(drop=True)])


    samples = samples.reset_index(drop=True)
    #print(samples)
    samples.to_csv('/tmp/cptac_samples.csv',index=False)
    return samples
   
def formatMutData(df,dtype,ctype,samp_names,source,samples):
    '''
    formats mutational data, which is different as its not a matrix!
    '''
    # some dataset has two level of indices some has only one
    if df.columns.nlevels == 2:
        df.columns = df.columns.droplevel(1)
    elif df.columns.nlevels != 1:
        print("The number of column levels is larger not 1 or 2!\n")
        raise

    ##firstr get relevant columnsz
    df = df.reset_index()
    subset = df[['Patient_ID','HGNC_Entrez_Gene_ID(supplied_by_NCBI)','Mutation','Genome_Change']]
    subset.columns=['Patient_ID','entrez_gene','variant_classification','Mutation']

    improve_mapping = samples
    improve_mapping.reset_index(drop=True)
    improve_mapping.index=improve_mapping.common_name

    blongdf = subset.join(improve_mapping,on='Patient_ID')
    blongdf.reset_index(drop=True)


    blongdf[['source']]=source
    blongdf[['study']]='CPTAC3'
    if(len(blongdf)) > 0:
        blongdf[['new_class']] =  blongdf[['variant_classification']].apply(reclass_var)
    else:
        blongdf[['new_class']] = blongdf[['variant_classification']]
    blongdf = blongdf.rename(columns={'variant_classification':'old_class',\
                                      'new_class':'variant_classification',\
                                      'entrez_gene':'entrez_id',
                                      'Mutation':'mutation'})
    blongdf = blongdf[['improve_sample_id','entrez_id','mutation','variant_classification','source','study']]
    print(blongdf)

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
      main function that takes a sample file and gene file as arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--prevSampleFile', dest='sampfile', 
                        default=None, help='Sample file to use to generate new ids. Returns sample file')
    parser.add_argument('--geneFile', dest='genefile',default='./genes.csv',
                       help='gene file to get gene ids')
    parser.add_argument('--curSampleFile', dest='newsamps',default=None, 
                        help='Sample file to use to generate data. Returns data for samples')
    opts = parser.parse_args()
    dat_files = {}
    ##read in existing gene identifier table from cell line data
    genes = pd.read_csv(opts.genefile)

    ####here is where we decide to build the data or not
    ## if there is an old sample file, we build only build sample file
    build_data=False
    if(opts.sampfile is not None):
        print(opts.sampfile)
        old_samples = pd.read_csv(opts.sampfile)
        samples = pd.DataFrame()
        print("generating sample file")
    elif(opts.newsamps is not None): ##otherwise we build the data
        build_data=True
        print(opts.newsamps)
        samples = pd.read_csv(opts.newsamps)
        print("Building dataset from generated sample file")
    else:
        print("Need a sample file to continue")
        exit()
    ##this loops through the 10 main cancer types and collects samples and then data
    for cancertype in ['brca','coad','hnscc','lscc','luad','ov','gbm','pdac','ucec','ccrcc']:
        dat = getCancerObj(cancertype)

        ##get available data for cancer
        #this call changed in recent version
        dat_list = dat.list_data_sources().set_index('Data type').to_dict()['Available sources']
        clinsource = dat_list['clinical']
        if 'harmonized' in clinsource:
            cs = 'harmonized'
    #    elif 'umich' in clinsource:
    #        cs = 'umich'
        else:
            cs = clinsource[0]
        tumor_samps = dat.get_clinical(cs)#.Sample_Tumor_Normal=='Tumor'
        tumor_samps = list(tumor_samps.index)
        if not build_data:
            print('Building sample file only for '+cancertype)
            if len(samples)==0:
                maxval = max(old_samples.improve_sample_id)
            else:
                maxval = max(samples.improve_sample_id)
            samples = buildTumorSampleTable(tumor_samps,cancertype,samples,maxval)
        if build_data:
            print('Building data file')
            all_dfs = {}
            all_sources = {} ##keep track of sources for long table
            ##all the data types we're collecting so far
            ## these are specifically the names of the file data keys in the `list_data_sources` function
            ## so have to match EXACTLY
            for dtype in ['somatic_mutation','proteomics','transcriptomics','CNV']: #'miRNA doesnt work
                if dtype not in dat_list.keys():
                    continue
                ###figure out whic source, prioritize harmonized when available
                source_list = dat_list[dtype]
                if 'harmonized' in source_list:
                    source = 'harmonized'
                elif 'umich' in dat_list[dtype]:
                    source = 'umich'
                else:
                    source = source_list[0]
                all_sources[dtype] = source #we can keep the source for future analysis/batch
                    #print(source)
                if dtype=='proteomics':
                    all_dfs[dtype] = dat.get_proteomics(source)
                if dtype=='transcriptomics':
                    all_dfs[dtype] = dat.get_transcriptomics(source)
                if dtype=='somatic_mutation':
                    all_dfs[dtype] = dat.get_somatic_mutation(source)
                if dtype=='miRNA':
                    all_dfs[dtype] = dat.get_miRNA(source)
                if dtype=='CNV':
                    all_dfs[dtype] = dat.get_CNV(source)

            for dtype,df in all_dfs.items():
                #tumor_samps = [t for t in df.index]  ##clinical sample info broke, so used this as a bypass
                df = df.loc[[t for t in tumor_samps if t in df.index]] ##get the data for those samples
                dfU = df
                dfU.dropna(how='all', axis=0, inplace=True)
                print(cancertype+' '+dtype)
                ##now we move to long form and match identifiers
                if dtype=='somatic_mutation': ## we actually want the column to be mutation
                    fdf = formatMutData(dfU,'mutation',cancertype,tumor_samps,all_sources[dtype],samples)
                    fdf = fdf.reset_index(drop=True)
                    dtype='mutations' ##but we want the file name to be mutations
                elif dtype=='CNV':##column should be copy_number
                    fdf = formatData(dfU,'copy_number',cancertype,tumor_samps,all_sources[dtype],genes,samples)
                    fdf = fdf.reset_index(drop=True)
                    dtype='copy_number'##so should file name
                else:
                    fdf = formatData(dfU,dtype,cancertype,tumor_samps,all_sources[dtype],genes,samples)
                    fdf = fdf.reset_index(drop=True)
                    
                if dtype in dat_files.keys():
                    of = dat_files[dtype]
                    fdf = pd.concat([of,fdf])
                    dat_files[dtype] = fdf
                else:
                    dat_files[dtype] = fdf
    print(build_data)
    if build_data:
        ##now concatenate all the cancers into a single file
        for dtype,df in dat_files.items():
            print('saving '+"cptac_"+dtype+'.csv.gz'+' file')
            df.to_csv("/tmp/"+"cptac_"+dtype+'.csv.gz',sep=',',index=False, compression='gzip')

if __name__ == '__main__':
    main()

