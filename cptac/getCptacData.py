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


##read in existing gene identifier table from cell line data
genes = pd.read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/genes.csv')

##read in samples, get max value
if os.path.exists('./samples.csv'):
    samples = pd.read_csv('./samples.csv')
    maxval = max(samples.improve_sample_id)
else:
    samples = pd.DataFrame.from_dict({'common_name':[],'cancer_type':[],'other_names':[],'species':[],\
                            'improve_sample_id':[],'id_source':[],'other_id':[],'model_type':[]})
    maxval = max(pd.read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/samples.csv').improve_sample_id)



def reclass_var(variant):
    # the dictionary started with
    # CCLE data
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
                      'Start_Codon_Ins':['Start_Codon_Ins'],\
                      'Stop_Codon_Del' :['Stop_Codon_Del'],\
                      'Stop_Codon_Ins' :['Stop_Codon_Ins'],\
                      'Silent':['Silent'],\
                      'Splice_Site':['Splice_Site'],\
                      'Translation_Start_Site':['Translation_Start_Site']}
    for key,var in variant_schema.items():
        if variant in var:
            return(key)
    print("Variant "+variant+' not found, returning undetermined')
    return('Undetermined')

def getCancerObj(cancertype):
   # cptac.download(dataset=cancertype,source='harmonized',)
    if cancertype == 'brca':
        dat = cptac.Brca()
    elif cancertype == 'ccrcc':
        dat = cptac.Ccrcc()
    elif cancertype == 'coad':
        dat = cptac.Coad()
    elif cancertype == 'endometrial':
        dat = cptac.Endometrial()
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


def buildTumorSampleTable(sample_names,cancer_type):
    '''
    we have to read in existing sample information first to see if it is part of the samples
    '''
    global maxval,samples
    
    for samp in sample_names:
       if samp not in samples.common_name:
           maxval = int(maxval+1)
           samples = pd.concat([samples,
                                pd.DataFrame({'common_name':[samp],\
                                              'cancer_type':[ct_vals[cancer_type]],\
                                              'other_names':[''],\
                                              'species':['Homo Sapiens (Human)'],\
                                              'improve_sample_id':[maxval],\
                                              'id_source':['CPTAC3'],\
                                              'other_id':[samp],\
                                              'model_type':['Tumor']}).reset_index(drop=True)])


    samples = samples.reset_index(drop=True)
    #print(samples)
    samples.to_csv('samples.csv',index=False)
    return samples
   
def formatMutData(df,dtype,ctype,samp_names,source):
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
    subset.columns=['Patient_ID','entrez_gene','variant_classification','mutation']

    improve_mapping = buildTumorSampleTable(samp_names,ctype)
    improve_mapping.reset_index(drop=True)
    improve_mapping.index=improve_mapping.common_name

    blongdf = subset.join(improve_mapping,on='Patient_ID')
    blongdf.reset_index(drop=True)
    
    blongdf[['source']]=source
    blongdf[['study']]='CPTAC3'
    print(blongdf)
    blongdf = blongdf[['improve_sample_id','entrez_gene','mutation','variant_classification','source','study']]
    return blongdf

def formatData(df,dtype,ctype,samp_names,source):
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
    improve_mapping = buildTumorSampleTable(samp_names,ctype)
    improve_mapping.reset_index(drop=True)
    improve_mapping.index=improve_mapping.common_name

    blongdf = longdf.join(improve_mapping,on='Patient_ID')
    improve_mapping.reset_index(drop=True)
    ##match gene identifiers

    genes.index=genes.gene_symbol
    mlongdf = blongdf.join(genes,on='gene_symbol',how='left',lsuffix='.e',rsuffix='.m') ##fix this
    
    mlongdf = mlongdf[['entrez_id','improve_sample_id',dtype]].drop_duplicates()

    ##if its copy number we need to include copy number call
    if dtype=='CNV': ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
       mlongdf[['copy_call']] = mlongdf[['CNV']].apply(copy_num)


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
   # parser = argparse.ArgumentParser()
   # parser.add_argument('--cancerType', dest='type',
   #                     help='Cancer type to be collected')
   # parser.add_argument('--sampleType', dest='sample', default='all',
   #                     help='Sample type, tumor vs normal vs all (default), \
   #                     to be collected')
   # opts = parser.parse_args()
    dat_files = {}

    for cancertype in ['brca','coad','hnscc','lscc','luad','ov','gbm','pdac','ucec','ccrcc']:
        dat = getCancerObj(cancertype)
        #dat = dat['mssm']
        ##get the tumor sample identifiers
        tumor_samps = dat.get_clinical()#.Sample_Tumor_Normal=='Tumor'
        tumor_samps = list(tumor_samps.index)

        ##get available data for cancer
        #this call changed in recent version
        dat_list = dat.list_data_sources().set_index('Data type').to_dict()['Available sources']
        #print(dat_list)
        all_dfs = {}
        all_sources = {} ##keep track of sources for long table
        ##all the data types we're collecting so far
        for dtype in ['somatic_mutation','proteomics','transcriptomics','CNV']: #'miRNA doesnt work
            if dtype not in dat_list.keys():
                continue
            ###figure out whic source, prioritize harmonized when available
            source_list = dat_list[dtype].split(',')
            if 'harmonized' in source_list:
                source = 'harmonized'
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
            df = df.loc[[t for t in tumor_samps if t in df.index]] ##get the data for those samples
            #dfE = np.exp(df)
            #dfU = np.log(dfE.sum(axis=1, level=0, min_count=1))
            dfU = df
            dfU.dropna(how='all', axis=0, inplace=True)
            print(cancertype+' '+dtype)
            ##now we move to long form and match identifiers
            if dtype=='somatic_mutation':
                fdf = formatMutData(dfU,dtype,cancertype,tumor_samps,all_sources[dtype]).reset_index(drop=True)
            else:
                fdf = formatData(dfU,dtype,cancertype,tumor_samps,all_sources[dtype]).reset_index(drop=True)
            if dtype in dat_files.keys():
                of = dat_files[dtype]
                fdf = pd.concat([of,fdf])
            dat_files[dtype] = fdf
    ##now concatenate all the cancers into a single file
    for dtype,df in dat_files.items():
        df.to_csv(dtype+'.csv.gz',sep=',',index=False, compression='gzip')
        
    #fdf.to_csv(path_or_buf=fname, sep=',',index=False)
   

if __name__ == '__main__':
    main()

