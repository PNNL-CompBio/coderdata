'''
python script to read in counts matrix and gene lengths to calculate tpm

'''

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser("Quick script to get TPM from counts matrix")
    parser.add_argument('--counts')


    args = parser.parse_args()
    print('Creating TPM from '+args.counts)
    counts = pd.read_csv(args.counts,sep='\t')
    counts.index=counts.stable_id

    ##get list of patients
    pats = set(counts.columns)-set(['stable_id','display_label','description','biotype'])

    ##transcript info from grc37
    # gtf = pd.read_csv("https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz",sep='\t',comment='#') # the "current" dir no longer exists...
    gtf = pd.read_csv("https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz",sep='\t',comment='#')
    gtf.index = [a.split(';')[0].split(' ')[1].strip('"') for a in gtf[gtf.columns[8]]]
    ##first select only exons
    gtf = gtf[gtf.gene=='exon']
    ##compute length and convert to kilobases
    length = abs(gtf[gtf.columns[4]]-gtf[gtf.columns[3]])/1000 ##get difference between start and end, then divide by kb

    #gtf = gtf[gtf.gene=='exon'].groupby(level=0).sum() ##sum over all the exons for a particular gene. yes this doesn't really tell you about which exon
    length = length.groupby(level=0).sum() ##sum over all exon lengths for each gene
    
    
    ##set counts matrix
    X = counts[list(pats)]
    tg = [g for g in X. index if g in length.index]
    
    X = X.loc[tg].transpose()
    length = pd.Series(length)[X.columns]
   # length =length.loc[tg]
    
    ##    df = pd.DataFrame(lengths=lengths,Genes=gnames)
    C = X.values
    L = length.values
    N = X.sum(axis=1).values.reshape(-1,1)
    rpk = C/L
    per_million_scaling_factor = (rpk.sum(axis=1)/1e6).reshape(-1,1)
    tpm = pd.DataFrame( rpk/per_million_scaling_factor, index=X.index, columns=X.columns).transpose()
    tpm.to_csv('tpm_'+args.counts,sep='\t')

if __name__=='__main__':
    main()
