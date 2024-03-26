import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample',dest='samplefile',default=None,help='DepMap sample file')
    parser.add_argument('--gene',dest='genefile',default=None,help='DepMap sample file')

    opts = parser.parse_args()

    samplefile = opts.samplefile
    gfile = opts.genefile

    samps = pd.read_csv(samplefile)
    print(samps)
    genes = pd.read_csv(gfile)[['gene_symbol','entrez_id']]
    genes = genes.drop_duplicates()

    print(genes)
    protfile='https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx'

    prots = pd.read_excel(protfile,'Normalized Protein Expression')
    print(prots)


    vvars = [col for col in prots.columns if '_Ten' in col]

    prot2 = prots[['Gene_Symbol']+vvars]

    ##we can just do the metlt here
    plong = pd.melt(prot2,id_vars='Gene_Symbol',value_vars=vvars,var_name='cellline',value_name='proteomics')

    ##rename gene symbol column
    plong = plong.rename({'Gene_Symbol':'gene_symbol'},axis=1)
    print(plong)
    
    ##split cell lin
    plong['other_id'] = [a.split('_Ten')[0] for a in plong.cellline]

    full = plong.merge(genes,on='gene_symbol')
    full = full.merge(samps,on='other_id')

    full = full[['entrez_id','proteomics','improve_sample_id']].drop_duplicates().dropna()

    full[['study']] = 'DepMap'
    full[['source']] = 'Broad'
    full.to_csv('/tmp/depmap_proteomics.csv',index=False)
    
main()
