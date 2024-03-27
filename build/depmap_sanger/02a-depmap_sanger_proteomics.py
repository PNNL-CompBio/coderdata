import pandas as pd
import argparse
from zipfile import ZipFile
import requests

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


    ##now get sanger
    sanger_protfile='https://cog.sanger.ac.uk/cmp/download/Proteomics_20221214.zip'
    r = requests.get(sanger_protfile)
    sanger_loc ='/tmp/sp.zip'
    open(sanger_loc , 'wb').write(r.content)

    zf = ZipFile(sanger_loc,'r')
    zf.extractall(path='/tmp/')
    pdat = pd.read_csv('/tmp/Protein_matrix_averaged_zscore_20221214.tsv',sep='\t',skiprows=[0])
    vv=pdat.columns[2:]
    plong = pd.melt(pdat,id_vars='symbol',value_vars=vv)
    pres = plong.rename({'symbol':'other_names','variable':'gene_symbol','value':'proteomics'},axis=1)
    pres = pres.merge(genes,on='gene_symbol')
    pres = pres.merge(samps,on='other_names')

    full2 = pres[['entrez_id','improve_sample_id','proteomics']]
    full2[['study']] = 'Sanger'
    full2[['source']] = 'Sanger'

    full3 = pd.concat([full,full2])
    print(full3)
    full3.dropna(axis=0)
    full3.to_csv('/tmp/depmap_sanger_proteomics.csv',index=False)
    
main()
