import pandas as pd
import numpy as np
import os
import math
import argparse
import synapseclient 
from pubchem_retrieval import update_dataframe_and_write_tsv
import warnings
warnings.filterwarnings("ignore")


def _get_max_old_smi(prev_df: pd.DataFrame) -> int:
    """
    Extract the numeric part from the improve_drug_id column like 'SMI_999' and return the max val.
    Returns 0 if nothing parseable is found.
    """
    if "improve_drug_id" not in prev_df.columns:
        return 0
    # pull digits out of SMI_###, ignore malformed
    extracted = (
        prev_df["improve_drug_id"]
        .astype(str)
        .str.extract(r"SMI_(\d+)", expand=False)
        .dropna()
    )
    if extracted.empty:
        return 0
    try:
        nums = extracted.astype(int)
        return int(nums.max())
    except ValueError:
        return 0
    

def _load_prev_drugs(prevDrugFilepath: str) -> pd.DataFrame:
    """
    Accepts a comma-separated list of previous drug file paths, reads each if it exists,
    and returns a deduplicated concatenated DataFrame. If none are readable, returns an
    empty DataFrame with a 'chem_name' column to avoid key errors.
    """
    if not prevDrugFilepath or prevDrugFilepath.strip() == "":
        return pd.DataFrame(columns=["chem_name"])

    paths = [p.strip() for p in prevDrugFilepath.split(",") if p.strip()]
    dfs = []
    for p in paths:
        if not os.path.exists(p):
            print(f"Warning: previous drug file '{p}' not found; skipping.")
            continue
        try:
            if p.lower().endswith(".tsv"):
                df = pd.read_csv(p, sep="\t")
            else:
                df = pd.read_csv(p)
            dfs.append(df)
        except Exception as e:
            print(f"Warning: failed to read previous drug file '{p}': {e}; skipping.")

    if not dfs:
        return pd.DataFrame(columns=["chem_name"])

    combined = pd.concat(dfs, ignore_index=True)
    # Drop exact duplicate rows to avoid inflated counts
    combined = combined.drop_duplicates()
    # Report what was loaded for debugging
    if "chem_name" in combined.columns:
        unique_prev = combined["chem_name"].nunique()
        print(f"Loaded {len(paths)} previous file(s); {unique_prev} unique previous drug names found.")
    else:
        print(f"Loaded {len(paths)} previous file(s), but 'chem_name' column missing in combined. Proceeding with empty previous set.")
        combined = pd.DataFrame(columns=["chem_name"])

    return combined



# function for loading data
def download_parse_drug_data(synID:str , save_path:str = None, synToken:str = None):
    """ 
    Download omics data from Synapse at synapseID syn66401303. Requires a synapse token, which requires you to make a Synapse account
    and create a Personal Access Token.  More information here: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens 
    
    Parameters
    ----------
    synID : string
        SynapseID of dataset to download. Default is synapseID of the omics dataset.
        
    save_path : string
        atal path where the downloaded file will be saved.

    synToken : string
        Synapse Personal Access Token of user.  Requires a Synapse account. More information at: https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens
        
    Returns
    -------
    mutations_data : pd.Dataframe
        A pandas dataframe containing mutations data

    copynum_data : pd.Dataframe
        A pandas dataframe containing copy number data

    proteomics_data : pd.Dataframe
        A pandas dataframe containing proteomics data
    """
    
    syn = synapseclient.Synapse() 
    syn.login(authToken=synToken) 
    
    # Obtain a pointer and download the data 
    downloaded_data = syn.get(entity=synID, downloadLocation = save_path) 

    # Get the path to the local copy of the data file 
    drugs_filepath = downloaded_data.path
    
    return(drugs_filepath)

def create_liverpdo_drug_data(drug_info_path: str, prevDrugFilepath: str, output_drug_data_path: str):
    """
    Using one or more previous drug files and the current liverpdo drug information,
    keep all prior entries, append only new drugs (with incremented SMI IDs), and
    write the full union to output_drug_data_path.
    """
    # Parse current liverpdo drug names (normalize to lowercase for matching) 
    drug_info_df = pd.read_csv(drug_info_path)
    liverpdo_drugs_df = pd.DataFrame({
        "chem_name": drug_info_df["Drug"].astype(str).str.lower().unique()
    })
    raw_names = set(liverpdo_drugs_df["chem_name"])

    # Load prior drug data union and determine max existing SMI index 
    prev_df = _load_prev_drugs(prevDrugFilepath)
    max_old = _get_max_old_smi(prev_df)
    print(f"Current max existing SMI index: {max_old}")

    # Determine which current names are already seen (case-insensitive) 
    seen_names = set()
    if not prev_df.empty and "chem_name" in prev_df.columns:
        seen_names = set(prev_df["chem_name"].astype(str).str.lower())
    new_names = [n for n in raw_names if n not in seen_names]

    # Prime temp file with full previous union 
    prime_file = output_drug_data_path.replace(".tsv", "_prime.tsv")
    # write previous entries (if any) into prime file
    to_write = prev_df.copy()
    to_write.to_csv(prime_file, sep="\t", index=False)

    # Append only new drugs 
    if new_names:
        update_dataframe_and_write_tsv(
            unique_names=new_names,
            output_filename=prime_file
        )
        print(f"Searched PubChem for {len(new_names)} new liverPDO drugs")
    else:
        print("No new liverPDO drugs to retrieve; only prior entries will be kept.")

    # Load back combined (prior + appended) and emit final union 
    combined = pd.read_csv(prime_file, sep="\t")
    combined.to_csv(output_drug_data_path, sep="\t", index=False)

    # cleanup
    try:
        os.remove(prime_file)
    except OSError:
        pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='###')

    # arguments for what data to process
    parser.add_argument('-d', '--Download', action = 'store_true', default=False, help='Download drug data.')
    parser.add_argument('-t', '--Token', type=str, default=None, help='Synapse Token')
    parser.add_argument('-D', '--Drug', action = 'store_true', default=False, help='Generate drug data.')
    parser.add_argument('-p', '--PrevDrugs', nargs='?', type=str, default='', const='', help='Previous drug file')

    args = parser.parse_args()


    ###########################

    if args.Download:
        if args.Token is None:
            print("No synpase download token was provided. Cannot download data.")
            exit()
        else:
            print("Downloading Files from Synapse.")
            # download fitted and raw drug data from synapse
            fitted_drug_data_path = download_parse_drug_data(synID = "syn66401300", save_path = "/tmp/", synToken = args.Token)
            drug_excel = pd.ExcelFile(open(fitted_drug_data_path, 'rb'))
            druginfo_df = pd.read_excel(drug_excel)
            druginfo_df.to_csv("/tmp/raw_druginfo.csv")
    if args.Drug:
        if args.PrevDrugs is None or args.PrevDrugs=='':
            print("No previous drugs file provided.  Starting improve_drug_id from SMI_1. Running drug file generation")
            create_liverpdo_drug_data(drug_info_path = "/tmp/raw_druginfo.csv", output_drug_data_path = "/tmp/liverpdo_drugs.tsv", prevDrugFilepath = "")
        else:
            print("Previous drugs file {} detected. Running drugs file generation and checking for duplicate IDs.".format(args.PrevDrugs))
            create_liverpdo_drug_data(drug_info_path = "/tmp/raw_druginfo.csv", prevDrugFilepath = args.PrevDrugs, output_drug_data_path = "/tmp/liverpdo_drugs.tsv")