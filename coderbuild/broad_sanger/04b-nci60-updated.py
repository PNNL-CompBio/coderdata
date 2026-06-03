'''
gets nci60 data from 10/2024 release

'''

import polars as pl
import argparse
import os
import re
from urllib import request

_BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)

def _retrieve_url(url: str, dest: str) -> None:
    """Download url with a browser User-Agent (NCI wiki blocks Python's default)."""
    req = request.Request(url, headers={"User-Agent": _BROWSER_UA})
    with request.urlopen(req) as resp, open(dest, "wb") as f:
        f.write(resp.read())

##oct 2024
#conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=13&modificationDate=1727922354561&api=v2'
#jan 2025
#conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=14&modificationDate=1735932462303&api=v2'
#may 2025
conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=19&modificationDate=1775183341937&api=v2'

#may 2025
oneconc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/ONECONC.zip?version=19&modificationDate=1775183401634&api=v2'

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleFile',dest='samplefile',default=None,help='DepMap sample file')
    parser.add_argument('--drugFile',dest='dfile',default=None,help='Drug database')


    opts = parser.parse_args()

    samplefile = opts.samplefile
    drugfile = opts.dfile
    if not os.path.exists('DOSERESP.csv'):
        _retrieve_url(conc_data, 'doseresp.zip')
        os.system('unzip doseresp.zip')

    samples = pl.read_csv(samplefile,quote_char='"')
    drugs = pl.read_csv(drugfile,separator='\t',quote_char='"')

    dose_resp = pl.read_csv("DOSERESP.csv",quote_char='"',infer_schema_length=10000000,ignore_errors=True)

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
    on = samples[['other_names','improve_sample_id']].rename({'other_names': 'common_name'})

    #there should be 71 cell lines, but there are 163.
    # 82 map to the 'other_names'
    # 81 map to neither
    sampmapping = pl.concat([on[['common_name','improve_sample_id']],samples[['common_name','improve_sample_id']]])

    sampmapping = sampmapping.unique().rename({'common_name': 'CELL_NAME'})

    ###create a time mapping tabel
    timemapping = pl.DataFrame(
        {
            'EXPID':dose_resp['EXPID'],
            'time':[72 if int(a[0:2])>22 and int(a[0:2])<50 and int(a[2:4])>0 else 48 for a in dose_resp['EXPID']],
            'time_unit':['hours' for a in dose_resp['EXPID']]
        }
        ).unique()


    ##now we can merge all the data into the dose response data frame
    merged = dose_resp[['AVERAGE_PTC','CONCENTRATION_UNIT','CONCENTRATION','CELL_NAME','EXPID','NSC']].join(sampmapping,on='CELL_NAME',how='left')
    merged = merged.join(timemapping,on='EXPID',how='left')

    ##clean up mssing samples
    nonulls = merged.filter(pl.col('improve_sample_id').is_not_null())

    nulls = merged.filter(pl.col('improve_sample_id').is_null())

    merged = nonulls

    ###we get a few more results added, but still missing a bunch
    merged = merged.join(drugmapping,on='NSC',how='left')
    nulldrugs = merged.filter(pl.col('improve_drug_id').is_null())
    nonulls =  merged.filter(pl.col('improve_drug_id').is_not_null())

    ###now update all the concentrations to be in Moles (some are in uM, all are log10)
    ##some are provided as molecular weights ('v') or other ('s') and we can't compare
    molar = merged.filter(pl.col('CONCENTRATION_UNIT')=='M')

    finaldf = pl.DataFrame(
        {
            'source':['NCI60_24' for a in molar['improve_drug_id']], ##2024 build
            'improve_sample_id':molar['improve_sample_id'],
            'Drug':molar['improve_drug_id'],
            # 'study': molar['EXPID'],#['NCI60' for a in nonulls['improve_drug_id']],
            'study': "NCI60",
            'time':molar['time'],
            'time_unit':molar['time_unit'],
            'DOSE': [(10**a)*1000000 for a in molar['CONCENTRATION']], ##move from molar to uM to match pharmacoDB
            'GROWTH':molar['AVERAGE_PTC']
        }
    )
    ##write to file
    finaldf.write_csv('nci60DoseResponse',separator='\t')


if __name__=='__main__':
    main()
