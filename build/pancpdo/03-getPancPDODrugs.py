import pandas as pd
import os
import argparse
import synapseclient as sc



###figshare link:

filelink='https://aacr.figshare.com/ndownloader/files/39996295'
synid = 'syn64333325'
##get third tab and drugsa re listeda cross top



def getDrugNames(token=""):
    if token !="":
        syn = sc.login(token)
    else:
        syn = sc.login()
    fpath = syn.get(synid).path
    print(fpath)
    tab = pd.read_excel(fpath,sheet_name='concentrations')
    drugs = [a.lower() for a in tab.columns]
    return drugs



def main():
    parser = argparse.ArgumentParser(description='Download and match pancpdodrugs')
    parser.add_argument('-p', '--pat',help='Synapse authentication token with permission to syn64333325')
    parser.add_argument('-d', '--prevDrugFile',help='Comma-delimited list of previous drug files')
    parser.add_argument('-o', '--output', default = '/tmp/pancpdo_drugs.tsv.gz')

    args = parser.parse_args()
    newdrugs = getDrugNames(args.pat)

    prevdrugs = [pd.read_csv(t,sep='\t') for t in args.prevDrugFile.split(',')]
    alldrugs = pd.concat(prevdrugs).drop_duplicates()

    imps = alldrugs[alldrugs.chem_name.isin(newdrugs)]
    newdrugs = alldrugs[alldrugs.improve_drug_id.isin(imps.improve_drug_id)]

    ##write drugs
    newdrugs.to_csv(args.output, sep='\t', compression='gzip', index=False)

    ##calculate drug descriptors
    
    
if __name__=='__main__':
    main()
