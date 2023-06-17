#!/usr/local/bin/python
'''
Basic CLI to import CPTAC protegenomic data
'''
import argparse
import cptac
import numpy as np

genes = ''
samples = ''

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


def buildTumorSampleTable(samplenames):
   '''
   we have to read in existing sample information first
   '''
   

def formatData(df): 
   '''
   formats data into long form
   '''
       # some dataset has two level of indices some has only one
    if df.columns.nlevels == 2:
        df.columns = df.columns.droplevel(1)
    elif df.columns.nlevels != 1:
        print("The number of column levels is larger not 1 or 2!\n")
        raise

     ##pivot wide to long format

     ##match sample identifiers or build new ones
     improve_mapping = buildTumorSampleTable(snames)
     
     ##match gene identifiers

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cancerType', dest='type',
                        help='Cancer type to be collected')
    parser.add_argument('--sampleType', dest='sample', default='all',
                        help='Sample type, tumor vs normal vs all (default), \
                        to be collected')
    opts = parser.parse_args()
    dat = getCancerObj(opts.type.lower())
  
    prot = dat.get_proteomics()
    mrna = dat.get_transcriptomics()
    mut = dat.get_somatic_mutation()
    mirna = dat.get_mirna()
    cnv = dat.get_cnv()
    clin = dat.get_clinical()

    
    # Get the sample type specific dataframe
    if opts.sample.lower() != 'all':
        meta = dat.get_clinical()
        if opts.sample.lower() == 'tumor':
            ind = meta[meta["Sample_Tumor_Normal"] == "Tumor"].index
            ind = [i for i in ind if i in df.index]
            df = df.loc[ind]
        elif opts.sample.lower() == 'normal':
            nIDs = list(meta[meta["Sample_Tumor_Normal"] == "Normal"].index)
            nIDs = list(set(nIDs) & set(df.index))
            df = df.loc[nIDs]
            df.index = [nID[:-2] if nID[-2:] ==
                        ".N" else nID for nID in nIDs]
        else:
            exit("The sample type, tumor vs normal vs all (default), \
            is not correctly set!")


    dfE = np.exp(df)
    dfU = np.log(dfE.sum(axis=1, level=0, min_count=1))
    dfU.dropna(how='all', axis=0, inplace=True)
    dfU.transpose().to_csv(path_or_buf="file.tsv", sep='\t')


if __name__ == '__main__':
    main()
