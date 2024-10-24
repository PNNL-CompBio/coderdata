

'''
gets nci60 drug information
'''

import polars as pl
import os
import argparse
import pubchem_retrieval as pr
import random as rand
from urllib import request

##drug files
smi_strings='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_smiles.csv?version=1&modificationDate=1710381820000&api=v2&download=true'
#oct 2024
smi_strings = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_smiles.csv?version=3&modificationDate=1727924130457&api=v2&download=true'

pc_ids='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_sid_cid.csv?version=2&modificationDate=1712766341112&api=v2&download=true'
pc_ids = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_sid_cid.csv?version=4&modificationDate=1727924129121&api=v2&download=true'

#oct 2024
chemnames='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_chemcal_name.csv?version=1&modificationDate=1710382716000&api=v2&download=true'

chemnames='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_chemical_name.csv?version=1&modificationDate=1727924127004&api=v2'
#oct 2024
cas='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_cas.csv?version=1&modificationDate=1710381783000&api=v2&download=true'
#oct 2024
cas = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_cas.csv?version=3&modificationDate=1727924126194&api=v2&download=true'
conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=11&modificationDate=1712351454136&api=v2'
##OCT 2024
conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=13&modificationDate=1727922354561&api=v2'


def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument('--test',action='store_true',default=False,help='Test script by sampling 100 chemicals')
    parser.add_argument('--output',default='/tmp/broad_sanger_drugs.tsv')
    opts = parser.parse_args()

    ###primary DF
    df = {'improve_drug_id':[],'chem_name':[],'canSMILES':[],'isoSMILES':[],\
          'InChIKey':[],'formula':[],'weight':[],'pubchem_id':[]}

    print('Downloading NSC identifiers for nci60 data')
    names = pl.read_csv(chemnames,ignore_errors=True)
    #castab = pd.read_csv(cas)
    pubchems = pl.read_csv(pc_ids)
    smiles = pl.read_csv(smi_strings)

    print('Getting experimental data to filter drugs')
    if not os.path.exists('DOSERESP.csv'):
        resp = request.urlretrieve(conc_data,'doseresp.zip')
        os.system('unzip doseresp.zip')
    dose_resp = pl.read_csv("DOSERESP.csv",quote_char='"',infer_schema_length=10000000,ignore_errors=True)
    pubchems = pubchems.filter(pl.col('NSC').is_in(dose_resp['NSC']))
    smiles = smiles.filter(pl.col("NSC").is_in(dose_resp['NSC']))
    ##first retreive pubchem data
    if opts.test:
        arr = rand.sample(list(pubchems['CID']),100)
    else:
        arr = set(pubchems['CID'])

    ##first filter to see if there are structures/drugs in teh data already. i dont think this does much.
    if os.path.exists(opts.output):
        curdrugs = pl.read_csv(opts.output,separator='\t')
       # cs = set(curdrugs['isoSMILES'])
        smiles = smiles.filter(pl.col('SMILES').is_not_null())
        upper=[a.upper() for a in smiles['SMILES']]
        smiles= pl.DataFrame({'NSC':smiles['NSC'],'upper':upper})#smiles.with_columns(upper=upper)
        ##reduce to smiels only in current drugs
        ssmiles = smiles.filter(~pl.col('upper').is_in(curdrugs['isoSMILES']))
        ssmiles = ssmiles.filter(~pl.col('upper').is_in(curdrugs['canSMILES']))
        pubchems = pubchems.filter(pl.col('NSC').is_in(ssmiles['NSC']))
        arr = set(pubchems['CID'])
        
    print("Querying pubchem from CIDs")
    pr.update_dataframe_and_write_tsv(arr,opts.output,'/tmp/ignore_chems.txt',batch_size=400,isname=False,time_limit=10*60*60)
    
    ##then make sure to paste `nsc` in front of all nsc idds
    res = pl.read_csv(opts.output,separator='\t')


    nsc = list(pubchems.filter(pl.col('CID').is_in(list(res['pubchem_id'])))['NSC'])

    print('Checking NSCs to see what we missed')
    missing = [n for n in nsc if 'nsc'+str(n) not in res['chem_name'] and 'nsc-'+str(n) not in res['chem_name']]
    
    ##check ignore_chems.txt
    print('missing '+str(len(missing))+' nsc ids')

    msmi = smiles.filter(pl.col('NSC').is_in(missing))
    print('Found SMILE strings for '+str(msmi.shape[1])+' NSCs')
    
    ##add in improve ids, nsc name and structure for all.
    mdf = msmi.join(names,on='NSC',how='left').join(pubchems,on='NSC',how='left')

    max_imp = max(int(a.split('_')[1]) for a in res['improve_drug_id'])

    smicount=len(set(mdf['SMILES'])) ## unique smiles in our missing data frame
    newdf = pl.DataFrame(
        {
            "improve_drug_id": ["SMI_"+str(a) for a in range(max_imp+1,max_imp+1+smicount,1)],
            'canSMILES': [a for a in set(mdf['SMILES'])],
            'isoSMILES': [a for a in set(mdf['SMILES'])],
            'InChIKey': [None for a in range(smicount)],
            'formula': [None for a in range(smicount)],
            'weight': [None for a in range(smicount)]
        }
    )

    #create updated nsc ids and names
    namedf = pl.DataFrame(
        {
            "nscid": ['nsc-'+str(a) for a in mdf['NSC']],
            'lower_name': [a if a is None else str(a).lower() for a in mdf['NAME']],
            'canSMILES': list(mdf['SMILES']),
            'pubchem_id': list(mdf['CID'])
        }
    )
    #merge and melt
    merged = pl.concat([mdf,namedf],how='horizontal').select(['SMILES','pubchem_id','nscid','lower_name'])
    melted = merged.melt(id_vars=['SMILES','pubchem_id'],value_vars=['nscid','lower_name']).select(['SMILES','pubchem_id','value']).unique()
    melted.columns = ['canSMILES','pubchem_id','chem_name']

    if newdf.shape[0] > 0:
        res = res.with_columns([
            pl.col("InChIKey").cast(pl.Utf8),
            pl.col("formula").cast(pl.Utf8)
        ])
        newdf = newdf.with_columns([
            pl.col("InChIKey").cast(pl.Utf8),
            pl.col("formula").cast(pl.Utf8)
        ])
    
        newdf = newdf.join(melted, on='canSMILES', how='inner').select(res.columns)
        res = pl.concat([res, newdf], how='vertical')
    res.write_csv(opts.output,separator='\t')
    
if __name__=='__main__':
    main()
