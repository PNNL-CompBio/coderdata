#!/usr/local/bin/python
'''
Basic CLI to import CPTAC protegenomic data
'''
import argparse
import cptac
import os

import pandas as pd


##read in existing gene identifier table from cell line data
genes = pd.read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/genes.csv')

##read in samples, get max value
if os.path.exists('./samples.csv'):
    samples = pd.read_csv('./samples.csv')
    maxval = max(samples.improve_sample_id)
else:
    samples = pd.DataFrame({'common_name':[],'cancer_type':[],'other_names':[],'species':[],\
                            'improve_sample_id':[],'id_source':[],'other_id':[],'model_type':[]})
    maxval = max(pd.read_csv('https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/data/samples.csv').improve_sample_id)
    

def getCancerObj(cancertype):
    cptac.download(dataset=cancertype)
    if cancertype == 'brca':
        dat = cptac.Brca()
    elif cancertype == 'ccrcc':
        dat = cptac.Ccrcc()
    elif cancertype == 'colon':
        dat = cptac.Colon()
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
    elif cancertype == 'ovarian':
        dat = cptac.Ovarian()
    else:
        exit()
    return dat


ct_vals={'brca':'Breast carcinoma','ccrcc':'Clear cell renal cell carcinoma','colon':'Colon adenocarcinoma',\
         'endometrial':'Endometrial carcinoma','gbm':'Glioblastoma','hnscc':'Head and neck squamous cell carcinoma',\
         'lscc':'Lung squamous cell carcinoma','luad':'Lung adenocarcinoma','ovarian':'Ovarian carcinoma'}
def buildTumorSampleTable(sample_names,cancer_type):
    '''
    we have to read in existing sample information first to see if it is part of the samples
    '''
    for samp in sample_names:
       if samp not in samples.common_name:
           maxval = maxval+1
           samples.append({'common_name':samp,'cancer_type':ct_vals[cancer_type],'other_names':'',\
                           species:'Homo Sapiens (Human)',\
                           'improve_sample_id':maxval,'id_source':'CPTAC3',\
                           'other_id':samp,'model_type':'Tumor'})
    samples.write_csv('samples.csv')
    return samples
   
   

def formatData(df,dtype,ctype,samp_names):
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
    improve_mapping.index=improve_mapping.common_name
    blongdf = longdf.join(improve_mapping,on='Patient_ID'))
     
    ##match gene identifiers
    genes.index=genes.gene_symbol
    mlongdf = nlongdf.join(genes,on='gene_symbol',how='left',lsuffix='.e',rsuffix='.m')) ##fix this

     

     
def main():
   # parser = argparse.ArgumentParser()
   # parser.add_argument('--cancerType', dest='type',
   #                     help='Cancer type to be collected')
   # parser.add_argument('--sampleType', dest='sample', default='all',
   #                     help='Sample type, tumor vs normal vs all (default), \
   #                     to be collected')
   # opts = parser.parse_args()

    for cancertype in ['brca','endometrial','ccrcc','hnscc','lscc','luad','ovarian']:
        dat = getCancerObj(cancertype)

        ##get the tumors
        tumor_samps = dat.get_clinical().Sample_Tumor_Normal=='Tumor'
        tumor_samps = list(tumor_samps.index)
        dat_list=dat.get_data_list()
        all_dfs = {}
        for dtype in ['proteomics','transcriptomics','somatic_mutations','miRNA','CNV']:
            if dtype not in dat_list.keys():
                continue
            if dtype=='proteomics':
                all_dfs[dtype] = dat.get_proteomics()
            if dtype=='transcriptomics':
                all_dfs[dtype] = dat.get_transcriptomics()
            if dtype=='somatic_mutations':
                all_dfs[dtype] = dat.get_somatic_mutation()
            if dtype=='miRNA':
                all_dfs[dtype] = dat.get_miRNA()
            if dtype=='CNV':
                all_dfs[dtype] = dat.get_CNV()

        for dtype,df in all_dfs.items():
#       df = df.loc[ind]
            dfE = np.exp(df)
            dfU = np.log(dfE.sum(axis=1, level=0, min_count=1))
            dfU.dropna(how='all', axis=0, inplace=True)

            ##now we move to long form and correct
            fdf = formatData(dfU,dtype,cancertype,samples)
            fname=cancertype+'_'+dtype+'.csv'
            dfU.transpose().to_csv(path_or_buf=fname, sep=',')

        
#    sample_type='all'
    # Get the sample type specific dataframe
#    if opts.sample.lower() != 'all':
#        if opts.sample.lower() == 'tumor':
#            ind = meta[meta["Sample_Tumor_Normal"] == "Tumor"].index
#            ind = [i for i in ind if i in df.index]
#            df = df.loc[ind]
#        elif opts.sample.lower() == 'normal':
#            nIDs = list(meta[meta["Sample_Tumor_Normal"] == "Normal"].index)
#            nIDs = list(set(nIDs) & set(df.index))
#            df = df.loc[nIDs]
            #df.index = [nID[:-2] if nID[-2:] ==
            #            ".N" else nID for nID in nIDs]
#        else:
#            exit("The sample type, tumor vs normal vs all (default), \
#            is not correctly set!")

   

if __name__ == '__main__':
    main()
