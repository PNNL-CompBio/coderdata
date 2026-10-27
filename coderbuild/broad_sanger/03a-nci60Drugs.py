

'''
gets nci60 drug information
'''

import polars as pl
import os
import argparse
import pubchem_retrieval as pr
import random as rand
import download_utils

_BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)

def _retrieve_url(url: str, dest: str) -> None:
    """Download url to dest with retries, resuming an interrupted transfer.

    NCI's wiki truncates downloads now and then (DOSERESP.zip in build v32,
    nsc_chemical_names.csv in the first 2.4 rebuild), so every file goes
    through download_utils rather than a single urlopen().read(). A browser
    User-Agent is sent because the wiki blocks Python's default.
    """
    download_utils.retrieve_url(url, dest, headers={"User-Agent": _BROWSER_UA})

def _read_csv_url(url: str, dest: str, **kwargs) -> pl.DataFrame:
    _retrieve_url(url, dest)
    return pl.read_csv(dest, **kwargs)

##drug files
smi_strings='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_smiles.csv?version=1&modificationDate=1710381820000&api=v2&download=true'
#oct 2024
smi_strings = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_smiles.csv?version=3&modificationDate=1727924130457&api=v2&download=true'

#jan 2025
smi_strings ='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_smiles.csv?version=5&modificationDate=1735932227387&api=v2'

pc_ids='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_sid_cid.csv?version=2&modificationDate=1712766341112&api=v2&download=true'
pc_ids = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_sid_cid.csv?version=4&modificationDate=1727924129121&api=v2&download=true'

#jan 2025
pc_ids = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_sid_cid.csv?version=5&modificationDate=1735931099618&api=v2'


chemnames = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_chemcal_name.csv?version=1&modificationDate=1710382716000&api=v2&download=true'
#oct 2024
chemnames = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_chemical_name.csv?version=1&modificationDate=1727924127004&api=v2'

#jan 2025
chemnames = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_chemical_names.csv?version=1&modificationDate=1735931084378&api=v2'

#oct 2024
cas='https://wiki.nci.nih.gov/download/attachments/155844992/nsc_cas.csv?version=1&modificationDate=1710381783000&api=v2&download=true'
#oct 2024
cas = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_cas.csv?version=3&modificationDate=1727924126194&api=v2&download=true'

#jan 2025
cas = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_cas.csv?version=4&modificationDate=1735931080720&api=v2'

##oct 2024
#conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=13&modificationDate=1727922354561&api=v2'
#jan 2025
#conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=14&modificationDate=1735932462303&api=v2'
#may 2025
conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=19&modificationDate=1775183341937&api=v2'

def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument('--test',action='store_true',default=False,help='Test script by sampling 100 chemicals')
    parser.add_argument('--output',default='/tmp/broad_sanger_drugs.tsv')
    opts = parser.parse_args()

    ###primary DF
    df = {'improve_drug_id':[],'chem_name':[],'canSMILES':[],\
          'InChIKey':[],'formula':[],'weight':[],'pubchem_id':[]}

    print('Downloading NSC identifiers for nci60 data')
    names = _read_csv_url(chemnames, 'nsc_chemical_names.csv', ignore_errors=True)
    #castab = pd.read_csv(cas)
    pubchems = _read_csv_url(pc_ids, 'nsc_sid_cid.csv')
    smiles = _read_csv_url(smi_strings, 'nsc_smiles.csv')

    print('Getting experimental data to filter drugs')
    if not os.path.exists('DOSERESP.csv'):
        _retrieve_url(conc_data, 'doseresp.zip')
        # -o overwrites without prompting; the exit status is checked so a
        # corrupt or truncated archive fails here rather than surfacing later
        # as a confusing parse error on a stale DOSERESP.csv.
        _rc = os.system('unzip -o doseresp.zip')
        if _rc != 0:
            raise RuntimeError(f"unzip of doseresp.zip failed (exit status {_rc}).")
    dose_resp = pl.read_csv("DOSERESP.csv",quote_char='"',infer_schema_length=10000000,ignore_errors=True)
    pubchems = pubchems.filter(pl.col('NSC').is_in(dose_resp['NSC']))
    nsc_cid = pubchems  # every NSC with dose data
    smiles = smiles.filter(pl.col("NSC").is_in(dose_resp['NSC'])
                           )
    ##first retreive pubchem data
    
    if opts.test:
        arr = rand.sample(list(pubchems['CID']),100)
    else:
        arr = set(pubchems['CID'])

    # A rerun after a failure finds this script's earlier output. Every CID is
    # still requested: update_dataframe_and_write_tsv skips CIDs already in the
    # output, and it also limits the output to the requested CIDs, so asking
    # for fewer (as this used to on a rerun) deleted drugs it had already
    # written -- about 560 NCI60 drugs in the second 2.4 rebuild.
    print("Querying pubchem from CIDs")
    pr.update_dataframe_and_write_tsv(
        unique_names=arr,
        output_filename=opts.output,
        ignore_chems="/tmp/ignore_chems.txt",
        batch_size=400,
        isname=False,
        time_limit=10*60*60,
        prev_drug_filepaths=None,
        restrict_to_raw_names=arr
    )
    
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

    res = add_missing_nscs(res, mdf)
    res = add_nsc_names(res, nsc_cid)
    res.write_csv(opts.output,separator='\t')


def _cid_str(column):
    """A CID column as text without a float suffix ('359323.0' -> '359323')."""
    return pl.col(column).cast(pl.Utf8).str.replace(r'\.0$', '')


def _drug_per_cid(res):
    """Return res with text columns, and one drug row (id + structure) per CID."""
    str_cols = ['canSMILES', 'InChIKey', 'formula', 'weight']
    res = res.with_columns([pl.col(c).cast(pl.Utf8) for c in str_cols] + [_cid_str('pubchem_id')])
    # One structure per CID: the lowest SMI id if a CID somehow has several.
    known = (
        res.with_columns(pl.col('improve_drug_id').str.extract(r'SMI_(\d+)', 1)
                         .cast(pl.Int64).alias('_n'))
           .sort('_n')
           .unique(subset=['pubchem_id'], keep='first', maintain_order=True)
           .select(['pubchem_id', 'improve_drug_id', 'canSMILES', 'InChIKey', 'formula', 'weight'])
    )
    return res, known


def add_nsc_names(res, nsc_cid):
    """
    Give every NSC an 'nsc-<id>' name under its CID's drug id.

    04b-nci60-updated.py joins dose responses to drugs on 'nsc-<id>' only.
    PubChem mostly lists NSCs as 'nsc123', and the drug cache no longer
    carries 'nsc-' names between builds, so without this most NCI60 drugs
    had no name 04b could match. NCI's NSC->CID table is used rather than
    PubChem synonyms, which can list an NSC under a related compound.
    """
    res, known = _drug_per_cid(res)
    names = (nsc_cid.select([_cid_str('CID').alias('pubchem_id'),
                             pl.format('nsc-{}', _cid_str('NSC')).alias('chem_name')])
                    .drop_nulls().unique())
    added = names.join(known, on='pubchem_id', how='inner').select(res.columns)
    out = pl.concat([res, added], how='vertical').unique(maintain_order=True)
    print(f'NSC names: {added.height} of {names.height} NSCs attached to a drug id')
    return out


def add_missing_nscs(res, mdf):
    """
    Attach NSC ids that PubChem does not list as synonyms to their drug.

    Every NSC in mdf maps (via NCI's NSC->CID table) to a CID already in res,
    so its 'nsc-<id>' name and NCI name are added under that CID's existing
    improve_drug_id and structure. Minting a new id here instead gave one
    chemical two ids; when both NSCs had dose data the build fit two curves
    that the stable-id step later collapsed onto one id (issue #512).
    """
    res, known = _drug_per_cid(res)
    nsc_names = mdf.select([_cid_str('CID').alias('pubchem_id'),
                            pl.format('nsc-{}', pl.col('NSC')).alias('chem_name')])
    # NCI's own compound names are only added when no other drug already uses
    # the name, so a PSet drug of that name cannot be matched to this CID.
    taken = res.select(['chem_name', 'pubchem_id']).unique()
    nci_names = (mdf.select([_cid_str('CID').alias('pubchem_id'),
                             pl.col('NAME').cast(pl.Utf8).str.to_lowercase().alias('chem_name')])
                    .join(taken.filter(pl.col('pubchem_id').is_not_null()).rename({'pubchem_id': '_other'}),
                          on='chem_name', how='left')
                    .filter(pl.col('_other').is_null() | (pl.col('_other') == pl.col('pubchem_id')))
                    .drop('_other'))
    names = pl.concat([nsc_names, nci_names]).drop_nulls().unique()
    added = names.join(known, on='pubchem_id', how='inner').select(res.columns)
    unmatched = names.join(known, on='pubchem_id', how='anti')
    if unmatched.height:
        print(f'{unmatched.height} NSC names had no CID in the drug table and were skipped')
    print(f'Attached {added.height} NSC names to existing drug ids')
    return pl.concat([res, added], how='vertical').unique(maintain_order=True)
    
if __name__=='__main__':
    main()
