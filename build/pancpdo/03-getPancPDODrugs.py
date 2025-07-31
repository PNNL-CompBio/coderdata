import pandas as pd
import os
import argparse
import synapseclient as sc
import pubchem_retrieval as pr

###figshare link:

filelink='https://aacr.figshare.com/ndownloader/files/39996295'
#synid = 'syn64333325'
##get third tab and drugsa re listeda cross top

##sup table drug list (in column names)
tablink = 'https://aacr.silverchair-cdn.com/aacr/content_public/journal/cancerdiscovery/8/9/10.1158_2159-8290.cd-18-0349/5/21598290cd180349-sup-199398_2_supp_4775187_p95dln.xlsx?Expires=1738004990&Signature=av8XadTm9AmI20O2Y7J7aHDtPbpluKJIfI5ubsoiYJ15D0zh5p1ltF4a7-DCSWTSMs-qX5TD09shxHeqkQ2NkLWHZsXoCD5KyREGhEgcDAvWZ1V9kwXDm0bjpINipAPPtC20oeuw6c~hPooF3Mtgzp4MzMCCjcVwfn05u27a0kS0yifBi11wQj3nmHlR3ym-2fYkFuqQtnNPCzH8-yIw21y0kTvXrNodAzC5pGA8qUK4PLxBt52xUIvTEPsPiPjXwBnDCfVsLGGdDYIY25lEPKiA403q6kFYvrSQ3bsTvM4kuvltb7yS4AXjK0-tthMOKbqq8~uREmJCcueADUF91g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'

def getDrugNames(token=""):
    #chemo drugs
    ctab = pd.read_excel(tablink,sheet_name=1,skiprows=1)
    #targeted drugs
    ttab = pd.read_excel(tablink,sheet_name=2,skiprows=1)
    drugs = [a.lower() for a in ctab.columns]+[a.lower() for a in ttab.columns]
    drugs = set(drugs)-set(['sample id','insensitive'])
    return drugs


def main():
    parser = argparse.ArgumentParser(description='Download and match pancpdodrugs')
#    parser.add_argument('-p', '--pat',help='Synapse authentication token with permission to syn64333325')
    parser.add_argument('-d', '--prevDrugFile', default=None, help='Comma-delimited list of previous drug files')
    parser.add_argument('-o', '--output', default = '/tmp/pancpdo_drugs.tsv')

    args = parser.parse_args()
    newdrugnames = getDrugNames()
    print(f"Raw pancpdo drug names ({len(newdrugnames)}): {sorted(newdrugnames)}")

    final_df = pr.update_dataframe_and_write_tsv(
        unique_names=newdrugnames,
        output_filename=args.output,
        batch_size=50,
        isname=True,
        prev_drug_filepaths=args.prevDrugFile if args.prevDrugFile and args.prevDrugFile.strip() else None,
        restrict_to_raw_names=newdrugnames
    )

    if final_df.empty:
        print("Warning: no pancpdo drugs were found.")
    else:
        kept_ids = set(final_df.get('improve_drug_id', []))
        print(f"Retained {len(final_df)} rows across {len(kept_ids)} improve_drug_id(s).")
    
    
if __name__=='__main__':
    main()
