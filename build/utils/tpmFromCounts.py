'''
python script to read in counts matrix and gene lengths to calculate tpm

'''

import argparse
import pandas as pd

def main(counts_data, genome_link, gene_column, exclude_columns, out_file):
    """
    Converts RNA count matrix to tpm matrix (transcripts per million).  
    
    Parameters
    ----------
    counts_data : string
        Path to RNA sequencing counts data. No default
        
    genome_link : string
        Link to human genome build. Defaults to "https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"

    gene_column : string
        Column name of column with gene name information.  Defaults to "stable_id".
    
    exclude_columns : string
        Column names of columns to exclude from patient list.  NO SPACES.  Defaults to "stable_id,display_label,description,biotype".

    out_file : string
        Path to output csv. No default.

    Returns
    -------
    None
    
    """
    # read in counts data
    counts = pd.read_csv(counts_data,sep='\t')
    counts.index=counts[gene_column]

    # parse list of columns to exclude from patients list
    exclude_cols_array = exclude_columns.split(",") # split long string by commas

    ##get list of patients
    pats = set(counts.columns)-set(counts.select_dtypes(include='object') + exclude_columns) # get patient names from column names, excluding columns were the datatype is a string and any columns in the exclude_columns arg


    ##transcript info from grc37
    # gtf = pd.read_csv("https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz",sep='\t',comment='#') # the "current" dir no longer exists...
    gtf = pd.read_csv(genome_link,sep='\t',comment='#')
    gtf.index = [a.split(';')[0].split(' ')[1].strip('"') for a in gtf[gtf.columns[8]]]
    ##first select only exons
    gtf = gtf[gtf.gene=='exon']
    ##compute length and convert to kilobases
    length = abs(gtf[gtf.columns[4]]-gtf[gtf.columns[3]])/1000 ##get difference between start and end, then divide by kb
    #gtf = gtf[gtf.gene=='exon'].groupby(level=0).sum() ##sum over all the exons for a particular gene. yes this doesn't really tell you about which exon
    length = length.groupby(level=0).sum() ##sum over all exon lengths for each gene

    
    ##set counts matrix
    X = counts[list(pats)]
    tg = [g for g in X.index if g in length.index]

    
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
    tpm.to_csv(out_file, sep='\t')


if __name__=='__main__':
    parser = argparse.ArgumentParser("Quick script to get TPM from counts matrix")

    parser.add_argument('--counts', default=None, help='Transcriptomics counts matrix')
    parser.add_argument('--genome_build', default="https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz", help='Link to human genome build')
    parser.add_argument('--gene_col', default="stable_id", help='Name of column with gene names')
    parser.add_argument('--exclude_col', default="stable_id,display_label,description,biotype", help='Name of column with gene names')
    parser.add_argument('--out_file', default=None, help='Output csv name.')


    args = parser.parse_args()
    print('Creating TPM from '+args.counts)
    main(counts_data = args.counts, genome_link = args.genome_build, gene_column = args.gene_col, exclude_columns = args.exclude_col, out_file = args.out_file)
