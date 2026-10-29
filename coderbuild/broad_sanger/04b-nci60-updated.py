'''
gets nci60 data from 10/2024 release

'''

import polars as pl
import argparse
import os
import re
import json
from urllib import request
import download_utils

_BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)

def _retrieve_url(url: str, dest: str) -> None:
    """Download url, resuming an interrupted transfer.

    DOSERESP.zip is ~346MB. The previous single urlopen().read() truncated at
    32MB in build v32 with IncompleteRead, which cost the entire nci60 dataset:
    this script exited 1, its caller only warned, nci60DoseResponse was never
    written, and nci60 shipped as eight empty files that still validated.
    """
    download_utils.retrieve_url(url, dest, headers={"User-Agent": _BROWSER_UA})

##oct 2024
#conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=13&modificationDate=1727922354561&api=v2'
#jan 2025
#conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=14&modificationDate=1735932462303&api=v2'
#may 2025
conc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/DOSERESP.zip?version=19&modificationDate=1775183341937&api=v2'

#may 2025
oneconc_data = 'https://wiki.nci.nih.gov/download/attachments/147193864/ONECONC.zip?version=19&modificationDate=1775183401634&api=v2'

#NCI's NSC -> PubChem CID table, the same release 03a-nci60Drugs.py uses
pc_ids = 'https://wiki.nci.nih.gov/download/attachments/155844992/nsc_sid_cid.csv?version=5&modificationDate=1735931099618&api=v2'


def _id_str(column):
    """An id column as text without a float suffix ('359323.0' -> '359323')."""
    return pl.col(column).cast(pl.Utf8).str.replace(r'\.0$', '')


def _load_nsc_cids():
    """Return NCI's NSC -> CID table as a DataFrame with str columns."""
    download_utils.retrieve_url(pc_ids, 'nsc_sid_cid.csv', headers={"User-Agent": _BROWSER_UA})
    tab = pl.read_csv('nsc_sid_cid.csv', ignore_errors=True)
    return tab.select([_id_str('NSC'), _id_str('CID')]).drop_nulls().unique()


def _pubchem_inchikeys(cids):
    """Return {CID: InChIKey} for the given CIDs from PubChem."""
    out = {}
    for i in range(0, len(cids), 100):
        batch = ','.join(cids[i:i + 100])
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{batch}/property/InChIKey/JSON'
        with request.urlopen(request.Request(url, headers={"User-Agent": _BROWSER_UA})) as resp:
            for prop in json.load(resp)['PropertyTable']['Properties']:
                out[str(prop['CID'])] = prop.get('InChIKey')
    return out


def resolve_nsc_ids(drugmapping, drugs, load_nsc_cids, fetch_inchikeys=_pubchem_inchikeys):
    """
    Keep one improve_drug_id per NSC.

    An 'nsc-<id>' name can sit under more than one drug id when PubChem lists
    the NSC as a synonym of a related compound. Joining on it copied that
    NSC's dose response to every such drug. Ambiguous NSCs keep the id whose
    PubChem CID matches NCI's own NSC->CID table (the lowest SMI id if several
    share that CID). If NCI's CID is none of the
    candidates (e.g. it is the non-stereo record of a salt), the candidate
    with the same InChIKey skeleton (first block) is kept. Any NSC still
    ambiguous stops the build.
    """
    counts = drugmapping.group_by('NSC').agg(pl.col('improve_drug_id').n_unique().alias('n'))
    ambiguous = counts.filter(pl.col('n') > 1)['NSC'].to_list()
    if len(ambiguous) == 0:
        return drugmapping
    nsc_cid = load_nsc_cids().select([_id_str('NSC'), _id_str('CID')]).drop_nulls().unique()
    id_cid = drugs.select([pl.col('improve_drug_id'), _id_str('pubchem_id')]).drop_nulls().unique()
    amb = (drugmapping.filter(pl.col('NSC').is_in(ambiguous))
           .with_columns(_id_str('NSC'))
           .join(id_cid, on='improve_drug_id', how='left')
           .join(nsc_cid.rename({'CID': 'nci_cid'}), on='NSC', how='left'))
    # Ids sharing NCI's CID are one chemical listed twice (the stable-id step
    # later merges them by structure), so the lowest SMI id is kept.
    keep = (amb.filter(pl.col('pubchem_id') == pl.col('nci_cid'))
               .select(['NSC', 'improve_drug_id']).unique()
               .with_columns(pl.col('improve_drug_id').str.extract(r'SMI_(\d+)', 1)
                             .cast(pl.Int64).alias('_n'))
               .sort(['NSC', '_n']).unique(subset=['NSC'], keep='first').drop('_n'))
    per = keep.group_by('NSC').agg(pl.col('improve_drug_id').n_unique().alias('n'))
    unresolved = {str(a) for a in ambiguous} - set(per.filter(pl.col('n') == 1)['NSC'])

    if unresolved:
        keep = keep.filter(~pl.col('NSC').is_in(list(unresolved)))
        rest = amb.filter(pl.col('NSC').is_in(list(unresolved)))
        nci_keys = fetch_inchikeys(sorted(set(rest['nci_cid'].drop_nulls().to_list())))
        id_keys = {}
        if 'InChIKey' in drugs.columns:
            for drug_id, key in drugs.select(['improve_drug_id', 'InChIKey']).drop_nulls().unique().rows():
                id_keys.setdefault(drug_id, key)
        picked = []
        for nsc in sorted(unresolved):
            cand = rest.filter(pl.col('NSC') == nsc)
            nci_key = nci_keys.get(cand['nci_cid'][0]) if cand['nci_cid'][0] is not None else None
            skeleton = nci_key.split('-')[0] if nci_key else None
            matches = {d for d in cand['improve_drug_id'].to_list()
                       if skeleton and (id_keys.get(d) or '').split('-')[0] == skeleton}
            if len(matches) == 1:
                picked.append((nsc, matches.pop()))
        if picked:
            keep = pl.concat([keep, pl.DataFrame(picked, schema=['NSC', 'improve_drug_id'], orient='row')])
        unresolved -= {nsc for nsc, _ in picked}
    if unresolved:
        raise ValueError(
            f"{len(unresolved)} NSC id(s) map to more than one drug and NCI's NSC->CID table "
            f"does not pick one: {sorted(unresolved)[:20]}")
    print(f'Resolved {len(ambiguous)} NSC ids that matched more than one drug id')
    keep = keep.with_columns(pl.col('NSC').cast(drugmapping['NSC'].dtype))
    return pl.concat([drugmapping.filter(~pl.col('NSC').is_in(ambiguous)), keep])


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleFile',dest='samplefile',default=None,help='DepMap sample file')
    parser.add_argument('--drugFile',dest='dfile',default=None,help='Drug database')


    opts = parser.parse_args()

    samplefile = opts.samplefile
    drugfile = opts.dfile
    if not os.path.exists('DOSERESP.csv'):
        _retrieve_url(conc_data, 'doseresp.zip')
        # -o overwrites without prompting; the exit status is checked so a
        # corrupt or truncated archive fails here rather than surfacing later
        # as a confusing parse error on a stale DOSERESP.csv.
        _rc = os.system('unzip -o doseresp.zip')
        if _rc != 0:
            raise RuntimeError(f"unzip of doseresp.zip failed (exit status {_rc}).")

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
    drugmapping = resolve_nsc_ids(drugmapping, drugs, _load_nsc_cids)

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
