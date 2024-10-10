'''
gets nci60 data from 10/2023 release

'''

import polars as pl
import argparse
#from zipfile import ZipFile
import os
#from io import BytesIO
import re
from urllib import request

conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=11&modificationDate=1712351454136&api=v2'
cancelled = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP_Cancelled.csv?version=1&modificationDate=1660871847000&api=v2&download=true'

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleFile',dest='samplefile',default=None,help='DepMap sample file') 
    parser.add_argument('--drugFile',dest='dfile',default=None,help='Drug database')
    
    opts = parser.parse_args()
    
    samplefile = opts.samplefile
    drugfile = opts.dfile
    if not os.path.exists('DOSERESP.csv'):
        resp = request.urlretrieve(conc_data,'doseresp.zip')
        os.system('unzip doseresp.zip')
    
    samples = pl.read_csv(samplefile,quote_char='"')
    drugs = pl.read_csv(drugfile,separator='\t',quote_char='"')
    
    dose_resp = pl.read_csv("DOSERESP.csv",quote_char='"',infer_schema_length=10000000)

    ##update drug mapping
    drugmapping = pl.DataFrame(
        {
            'chem_name' : ['nsc-'+str(nsc) for nsc in set(dose_resp['NSC'])],
            'NSC' : [a for a in set(dose_resp['NSC'])]
        }
    )

    drugmapping = drugmapping.join(drugs,on='chem_name')[['NSC','improve_drug_id']]
    drugmapping = drugmapping.unique()

    ###update sample mapping
    on = samples[['other_names','improve_sample_id']]
    on.columns=['common_name','improve_sample_id']
    
    sampmapping = pl.concat([on[['common_name','improve_sample_id']],samples[['common_name','improve_sample_id']]])
                            
    sampmapping = sampmapping.unique()
    sampmapping.columns = ['CELL_NAME','improve_sample_id']

    ###create a time mapping tabel
    timemapping = pl.DataFrame(
        {
            'EXPID':dose_resp['EXPID'],
            'time':[72 if int(a[0:2])>22 and int(a[0:2])<50 and int(a[2:4])>0 else 48 for a in dose_resp['EXPID']],
            'time_unit':['hours' for a in dose_resp['EXPID']]
        }
        ).unique()


    ##now we can merge all the data into the dose response data frame
    merged = dose_resp[['AVERAGE_PTC','CONCENTRATION','CELL_NAME','EXPID','NSC']].join(sampmapping,on='CELL_NAME',how='left')
    merged = merged.join(timemapping,on='EXPID',how='left')
    
    ##clean up mssing samples
    nonulls = merged.filter(pl.col('improve_sample_id').is_not_null())

    nulls = merged.filter(pl.col('improve_sample_id').is_null())
    
    newnames = pl.DataFrame(
        {
            'new_name': [re.split(r' |\(|/', a)[0] for a in nulls['CELL_NAME']],
            'CELL_NAME':nulls['CELL_NAME']
        }
    )
    newnames = newnames.unique()
    
    fixed = nulls[['AVERAGE_PTC','CONCENTRATION','CELL_NAME','EXPID','NSC','time','time_unit']].join(newnames,on='CELL_NAME',how='left')
    fixed.columns = ['AVERAGE_PTC','CONCENTRATION','old_CELL_NAME','EXPID','NSC','time','time_unit','CELL_NAME']
    fixed = fixed.join(sampmapping,on='CELL_NAME',how='left')[['AVERAGE_PTC','CONCENTRATION','old_CELL_NAME','EXPID','NSC','improve_sample_id','time','time_unit']]
    fixed.columns = ['AVERAGE_PTC','CONCENTRATION','CELL_NAME','EXPID','NSC','improve_sample_id','time','time_unit']
    fixed = fixed.filter(pl.col('improve_sample_id').is_not_null())

    merged = pl.concat([nonulls,fixed])
    
    ###we get a few more results added, but still missing a bunch    
    merged = merged.join(drugmapping,on='NSC',how='left')
    nulldrugs = merged.filter(pl.col('improve_drug_id').is_null())
    nonulls =  merged.filter(pl.col('improve_drug_id').is_not_null())
    finaldf = pl.DataFrame(
        {
            'source':['NCI60_24' for a in nonulls['improve_drug_id']], ##2024 build
            'improve_sample_id':nonulls['improve_sample_id'],
            'Drug':nonulls['improve_drug_id'],
            'study':['NCI60' for a in nonulls['improve_drug_id']],
            'time':nonulls['time'],
            'time_unit':nonulls['time_unit'],
            'DOSE': [10**a for a in nonulls['CONCENTRATION']],
            'GROWTH':nonulls['AVERAGE_PTC']
        }
    )
    ##write to file
    finaldf.write_csv('nci60DoseResponse',separator='\t')
    

if __name__=='__main__':
    main()
