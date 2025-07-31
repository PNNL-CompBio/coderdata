import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import threading
import time
import signal
import sys

# Global variables
request_counter = 0
last_request_time = time.time()
lock = threading.Lock()
should_continue = True
improve_drug_id = 0
existing_synonyms = set()
existing_structures = dict()
existing_pubchemids = set()


def fetch_url(url, retries=3, backoff_factor=1):
    """
    Fetches a URL with retry mechanism and backoff.
    
    Parameters:
    - url (str): The URL to fetch.
    - retries (int): Number of retry attempts.
    - backoff_factor (float): Factor to calculate backoff time.
    
    Returns:
    - dict: JSON response if successful.
    
    Raises:
    - Exception: If all retry attempts fail.
    """
    global last_request_time, lock, request_counter
    with lock:
        current_time = time.time()
        # Reset counter if more than 1 second has passed
        if current_time - last_request_time >= 1:
            request_counter = 0
            last_request_time = current_time
        
        # Wait if the limit is reached
        while request_counter >= 4:
            time.sleep(0.2)  # Sleep a bit to check again
            current_time = time.time()
            if current_time - last_request_time >= 1:
                request_counter = 0
                last_request_time = current_time
        
        request_counter += 1
    
    for attempt in range(retries + 1):  # Total attempts = retries + 1
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                return response.json()
            else:
                raise Exception(f"Failed to fetch {url}, Status Code: {response.status_code}")
        except Exception as exc:
            if attempt < retries:
                wait = backoff_factor * (2 ** attempt)
                print(f"Attempt {attempt + 1} for URL {url} failed with error: {exc}. Retrying in {wait} seconds...")
                time.sleep(wait)
            else:
                print(f"All {retries + 1} attempts failed for URL {url}.")
                raise


def retrieve_drug_info(compound, ignore_chems, isname=True):
    """
    Retrieves information for a given compound from PubChem.

    Parameters:
    - compound (str or int): Name or CID of the compound.
    - ignore_chems (str): File path to log ignored compounds.
    - isname (bool): True if the compound is a name, False if it's a CID.
    
    Returns:
    - list: List of dictionaries containing drug information, or None if unsuccessful.
    """
    global improve_drug_id, existing_synonyms, existing_structures
    if pd.isna(compound):
        return None

    if isname:
        urls = {
            "properties": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/property/SMILES,InChIKey,MolecularFormula,MolecularWeight/JSON",
            "synonyms": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/synonyms/JSON"
        }
    else:
        urls = {
            "properties": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{compound}/property/SMILES,InChIKey,MolecularFormula,MolecularWeight/JSON",
            "synonyms": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{compound}/synonyms/JSON"
        }

    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_url = {executor.submit(fetch_url, url): key for key, url in urls.items()}
        results = {}

        for future in as_completed(future_to_url):
            key = future_to_url[future]
            try:
                data = future.result()
                results[key] = data
            except Exception as exc:
                print(f'{compound} generated an exception: {exc}')
                with open(ignore_chems, "a") as f:
                    f.write(f"{compound}\n")
                return None

    if all(key in results for key in ["properties", "synonyms"]):
        properties = results["properties"]['PropertyTable']['Properties'][0]
        synonyms_list = results["synonyms"]['InformationList']['Information'][0]['Synonym']

        # Check if this compound or any of its synonyms already has an assigned improve_drug_id
        new_syns = set()
        if isname:
            sl = synonyms_list + [compound]
        else:
            sl = synonyms_list
        for synonym in sl:
            synonym_lower = str(synonym).lower()
            if synonym_lower not in existing_synonyms:
                new_syns.add(synonym_lower)
        if len(new_syns) == 0:  # Ensure there are new synonyms before proceeding
            return None
        for synonym in new_syns:
            synonym_lower = str(synonym).lower()
            existing_synonyms.add(synonym_lower)

        # Check for structure
        if properties['SMILES'] in existing_structures.keys():
            print(f'Found structure for {compound}')
            SMI_assignment = existing_structures[properties['SMILES']]
        else:
            improve_drug_id += 1
            SMI_assignment = f"SMI_{improve_drug_id}"
            existing_structures[properties['SMILES']] = SMI_assignment
        
        #print(new_syns)
        data_for_tsv = [{
            'improve_drug_id': SMI_assignment,
            'name': str(synonym).lower(),
            **properties
        } for synonym in new_syns]

        return data_for_tsv
    else:
        return None


def fetch_data_for_batch(batch, ignore_chems, isname):
    """
    Fetches drug information for a batch of compounds.

    Parameters:
    - batch (list): List of compound names or CIDs.
    - ignore_chems (str): File path to log ignored compounds.
    - isname (bool): True if compounds are names, False if they're CIDs.
    
    Returns:
    - list: Combined list of drug information for the batch.
    """
    all_data = []
    for compound_name in batch:
        data = retrieve_drug_info(compound_name, ignore_chems, isname)
        if data:
            all_data.extend(data)
    return all_data


def read_existing_data(output_filename):
    """
    Reads existing data from the output file to prevent duplication.

    Parameters:
    - output_filename (str): File path to the output file.
    
    Returns:
    - None
    """
    global improve_drug_id, existing_synonyms, existing_structures, existing_pubchemids
    try:
        df = pd.read_csv(output_filename, sep='\t', quoting=3)
        existing_synonyms = set([str(a).lower() for a in set(df.chem_name)])
        existing_pubchemids = set([str(a) for a in df['pubchem_id']])
        max_id = df['improve_drug_id'].str.extract(r'SMI_(\d+)').astype(float).max()
        improve_drug_id = int(max_id[0]) + 1 if pd.notna(max_id[0]) else 1
        existing_structures = {row['canSMILES']: row['improve_drug_id'] for _, row in df.iterrows()}
        print(f'Read in {len(existing_synonyms)} drug names and {len(existing_pubchemids)} pubchem IDs')
    except FileNotFoundError:
        return {}


def timeout_handler(signum, frame):
    """
    Handles timeouts by setting the global `should_continue` flag to False.
    """
    global should_continue
    print("Time limit reached, exiting gracefully...")
    should_continue = False


# def update_dataframe_and_write_tsv(unique_names, output_filename="drugs.tsv", ignore_chems="ignore_chems.txt",
#                                    batch_size=1, isname=True, time_limit=48 * 60 * 60):
#     """
#     Updates the data frame with drug information and writes it to a TSV file.

#     Parameters:
#     - unique_names (iterable): List of unique compound names or CIDs.
#     - output_filename (str): File path to the output TSV file.
#     - ignore_chems (str): File path to log ignored compounds.
#     - batch_size (int): Number of compounds to process in each batch.
#     - isname (bool): True if unique_names are names, False if they're CIDs.
#     - time_limit (int): Time limit for the script in seconds. This is a remnant of the GitHub Action CI.
    
#     Returns:
#     - None
#     """
#     global should_continue, existing_synonyms, existing_pubchemids
#     signal.signal(signal.SIGALRM, timeout_handler)
#     signal.alarm(time_limit)
#     print(f'Starting with {len(unique_names)} unique drug names/IDs')

#     try:
#         print(f'Reading existing data from {output_filename}')
#         read_existing_data(output_filename)
#         if isname:
#             unique_names = set([str(name).lower() for name in unique_names if not pd.isna(name)])
#             unique_names = set(unique_names) - set(existing_synonyms)
#             print(f'Looking at {len(unique_names)} names')
#         else:
#             unique_names = set([str(name) for name in unique_names if not pd.isna(name)])
#             unique_names = set(unique_names) - set(existing_pubchemids)
#             print(f'Looking at {len(unique_names)} IDs')
#         ignore_chem_set = set()
#         if os.path.exists(ignore_chems):
#             with open(ignore_chems, 'r') as file:
#                 for line in file:
#                     ignore_chem_set.add(line.strip())
#         unique_names = list(set(unique_names) - ignore_chem_set)

#         print(f"{len(unique_names)} Drugs to search")
#         for i in range(0, len(unique_names), batch_size):
#             if not should_continue:
#                 break
#             if unique_names[i] in existing_synonyms or unique_names[i] in existing_pubchemids:
#                 continue

#             batch = unique_names[i:i + batch_size]
#             data = fetch_data_for_batch(batch, ignore_chems, isname)
#             if data:
#                 file_exists = os.path.isfile(output_filename)
#                 mode = 'a' if file_exists else 'w'
#                 with open(output_filename, mode) as f:
#                     if not file_exists:
#                         f.write("improve_drug_id\tchem_name\tpubchem_id\tcanSMILES\tInChIKey\tformula\tweight\n")
#                     for entry in data:
#                         f.write(f"{entry['improve_drug_id']}\t{entry['name']}\t{entry.get('CID', '')}\t"
#                                 f"{entry['SMILES']}\t{entry['InChIKey']}\t"
#                                 f"{entry['MolecularFormula']}\t{entry['MolecularWeight']}\n")

#                 with open(ignore_chems, "a") as ig_f:
#                     for entry in data:
#                         if isname:
#                             ig_f.write(f"{entry['name']}\n")
#                         else:
#                             ig_f.write(f"{entry.get('CID', '')}\n")

#     except Exception as e:
#         print(f"An unexpected error occurred: {e}")
#     finally:
#         signal.alarm(0)





def _load_prev_drugs_union(prevDrugFilepath: str) -> pd.DataFrame:
    """
    Load and concatenate comma-separated prior drug TSVs, deduplicate, and return.
    """
    if not prevDrugFilepath or str(prevDrugFilepath).strip() == "":
        return pd.DataFrame(columns=["improve_drug_id", "chem_name", "pubchem_id", "canSMILES", "InChIKey", "formula", "weight"])

    paths = [p.strip() for p in str(prevDrugFilepath).split(",") if p.strip()]
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
        return pd.DataFrame(columns=["improve_drug_id", "chem_name", "pubchem_id", "canSMILES", "InChIKey", "formula", "weight"])

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.drop_duplicates()
    return combined


def _max_smi_in_df(df: pd.DataFrame) -> int:
    """
    Extract max numeric part of improve_drug_id like SMI_123 from a dataframe.
    """
    if "improve_drug_id" not in df.columns:
        return 0
    extracted = df["improve_drug_id"].astype(str).str.extract(r"SMI_(\d+)", expand=False)
    nums = pd.to_numeric(extracted, errors="coerce")
    if nums.empty or nums.dropna().empty:
        return 0
    return int(nums.max())


# --- revised main function --- #

def update_dataframe_and_write_tsv(unique_names,
                                   output_filename="drugs.tsv",
                                   ignore_chems="ignore_chems.txt",
                                   batch_size=1,
                                   isname=True,
                                   time_limit=48 * 60 * 60,
                                   prev_drug_filepaths=None,
                                   restrict_to_raw_names=None):
    """
    Updates the data frame with drug information and writes it to a TSV file.

    New features:
    - Accepts previous drug file(s) via `prev_drug_filepaths` (comma-separated) to temp existing entries and
      continue SMI numbering from the global max across those and the existing output.
    - Only retains drugs relevant to `restrict_to_raw_names` (e.g., liverpdo/bladderpdo raw drug names).
    - Avoids re-querying names/IDs already present in either previous files or existing output.
    
    Parameters:
    - unique_names (iterable): Current raw compound names or CIDs to consider for this dataset.
    - output_filename (str): Final filtered output TSV path.
    - ignore_chems (str): File path to log ignored compounds.
    - batch_size (int): Number of compounds to process in each batch.
    - isname (bool): True if unique_names are names, False if they are CIDs.
    - time_limit (int): Timeout in seconds.
    - prev_drug_filepaths (str or None): Comma-separated prior drug TSV file paths.
    - restrict_to_raw_names (iterable or None): If provided, final output is filtered to only these names (lowercased for names, raw for CIDs).
    
    Returns:
    - pd.DataFrame: The final written DataFrame (subset of relevant drugs).
    """
    global should_continue, existing_synonyms, existing_pubchemids, improve_drug_id
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(time_limit)

    # Normalize input raw names
    if isname:
        raw_names = {str(n).strip().lower() for n in unique_names if not pd.isna(n)}
    else:
        raw_names = {str(n).strip() for n in unique_names if not pd.isna(n)}
    if restrict_to_raw_names is not None:
        if isname:
            restrict_set = {str(n).strip().lower() for n in restrict_to_raw_names if not pd.isna(n)}
        else:
            restrict_set = {str(n).strip() for n in restrict_to_raw_names if not pd.isna(n)}
    else:
        restrict_set = raw_names  # default filtering

    print(f"Starting with {len(raw_names)} provided {'names' if isname else 'IDs'}; restricting output to {len(restrict_set)} of them.")

    # --- 1) read existing output to bootstrap state ---
    print(f"Reading existing data from {output_filename}")
    # capture existing output file (if any) to include in base
    existing_output_df = pd.DataFrame()
    if os.path.exists(output_filename):
        try:
            existing_output_df = pd.read_csv(output_filename, sep="\t", quoting=3)
        except Exception:
            existing_output_df = pd.read_csv(output_filename, sep="\t")
    # read_existing_data populates globals (synonyms, pubchemids, and sets improve_drug_id based on output)
    read_existing_data(output_filename)
    existing_output_max = improve_drug_id - 1  # because improve_drug_id was set to last+1

    # --- 2) load previous union and incorporate its names/IDs into seen sets ---
    prev_union_df = _load_prev_drugs_union(prev_drug_filepaths)
    prev_union_max = _max_smi_in_df(prev_union_df)

    # adjust improve_drug_id to be max of existing output and previous union, so new IDs start after both
    desired_start = max(existing_output_max, prev_union_max) + 1
    if improve_drug_id < desired_start:
        improve_drug_id = desired_start
    print(f"SMI numbering will start from {improve_drug_id} (max prior was {desired_start - 1})")

    # build seen names/IDs (to avoid re-query)
    seen_names = set(existing_synonyms)
    seen_pubchemids = set(existing_pubchemids)
    if not prev_union_df.empty:
        if "chem_name" in prev_union_df.columns:
            seen_names.update({str(n).strip().lower() for n in prev_union_df["chem_name"].astype(str)})
        if "pubchem_id" in prev_union_df.columns:
            seen_pubchemids.update({str(n).strip() for n in prev_union_df["pubchem_id"].astype(str)})

    # --- 3) determine new candidates to query ---
    if isname:
        candidates = raw_names - seen_names
        print(f"{len(raw_names)} raw names provided; {len(seen_names)} already seen; {len(candidates)} new to fetch.")
    else:
        candidates = raw_names - seen_pubchemids
        print(f"{len(raw_names)} raw IDs provided; {len(seen_pubchemids)} already seen; {len(candidates)} new to fetch.")

    # apply ignore_chems filtering
    ignore_chem_set = set()
    if os.path.exists(ignore_chems):
        with open(ignore_chems, "r") as file:
            for line in file:
                ignore_chem_set.add(line.strip())
    candidates = set(candidates) - ignore_chem_set
    print(f"{len(candidates)} candidates remain after removing ignored.")

    # --- 4) make a temp union file with previous union + existing output ---
    temp_file = output_filename.replace(".tsv", "_temp.tsv")
    base_dfs = []
    if not prev_union_df.empty:
        base_dfs.append(prev_union_df)
    if not existing_output_df.empty:
        base_dfs.append(existing_output_df)
    if base_dfs:
        base_union = pd.concat(base_dfs, ignore_index=True).drop_duplicates()
    else:
        base_union = pd.DataFrame()
    # Write that tempd base for appending
    if not base_union.empty:
        with open(temp_file, "w") as f:
            header_written = False
            for _, row in base_union.iterrows():
                if not header_written:
                    cols = row.index.tolist()
                    f.write("\t".join(cols) + "\n")
                    header_written = True
                f.write("\t".join(str(row[col]) if pd.notna(row[col]) else "" for col in row.index) + "\n")
    else:
        # create empty temp file so fetch logic can append headers
        open(temp_file, "a").close()

    # --- 5) fetch new ones in batches and append to temp_file ---
    candidates_list = list(candidates)
    for i in range(0, len(candidates_list), batch_size):
        if not should_continue:
            break
        batch = candidates_list[i : i + batch_size]
        data = fetch_data_for_batch(batch, ignore_chems, isname)
        if data:
            file_exists = os.path.isfile(temp_file)
            mode = "a" if file_exists else "w"
            with open(temp_file, mode) as f:
                if os.path.getsize(temp_file) == 0:
                    f.write("improve_drug_id\tchem_name\tpubchem_id\tcanSMILES\tInChIKey\tformula\tweight\n")
                for entry in data:
                    f.write(
                        f"{entry['improve_drug_id']}\t{entry['name']}\t{entry.get('CID', '')}\t"
                        f"{entry['SMILES']}\t{entry['InChIKey']}\t"
                        f"{entry['MolecularFormula']}\t{entry['MolecularWeight']}\n"
                    )
            with open(ignore_chems, "a") as ig_f:
                for entry in data:
                    if isname:
                        ig_f.write(f"{entry['name']}\n")
                    else:
                        ig_f.write(f"{entry.get('CID', '')}\n")

    # --- 6) load combined temp results ---
    combined = pd.read_csv(temp_file, sep="\t")

    # Determine previous max (before new fetches) to identify newly assigned SMI IDs
    previous_max = desired_start - 1  # desired_start was max(existing_output_max, prev_union_max) + 1

    # --- 7) compute hit improve_drug_id(s) from restrict_set (preserves all synonyms) ---
    hit_ids = set()
    if isname:
        mask_hit = combined["chem_name"].astype(str).str.lower().isin(restrict_set)
        hit_ids = set(combined.loc[mask_hit, "improve_drug_id"])
    else:
        if "pubchem_id" in combined.columns:
            mask_hit = combined["pubchem_id"].astype(str).isin(restrict_set)
            hit_ids = set(combined.loc[mask_hit, "improve_drug_id"])

    # --- 8) identify newly assigned improve_drug_id(s) ---
    new_ids = set()
    if "improve_drug_id" in combined.columns:
        extracted_comb = combined["improve_drug_id"].astype(str).str.extract(r"SMI_(\d+)", expand=False)
        nums_comb = pd.to_numeric(extracted_comb, errors="coerce")
        if not nums_comb.empty:
            new_ids = set(combined.loc[nums_comb > previous_max, "improve_drug_id"])
            if new_ids:
                print(f"Newly assigned improve_drug_id(s): {new_ids}")

    # --- 9) union and filter final DataFrame by improve_drug_id(s) ---
    keep_ids = hit_ids.union(new_ids)
    if keep_ids:
        final_df = combined[combined["improve_drug_id"].isin(keep_ids)].copy()
    else:
        print("Warning: no relevant drugs were retained/fetched for the restriction set.")
        final_df = pd.DataFrame(columns=combined.columns)

    # --- 10) write final filtered output ---
    final_df.to_csv(output_filename, sep="\t", index=False)

    if os.path.exists(temp_file):
        try:
            os.remove(temp_file)
        except OSError as e:
            print(f"Warning: failed to delete temp file {temp_file}: {e}")


    return final_df
