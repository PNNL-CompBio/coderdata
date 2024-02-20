import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
import argparse

def fetch_url(url):
    """Helper function to fetch a URL."""
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to fetch {url}")

def retrieve_drug_info2(compound_name, improve_drug_id, existing_synonyms):
    lower_compound_name = compound_name.lower()
    if pd.isna(compound_name) or lower_compound_name in existing_synonyms:
        return None

    urls = {
        "properties": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight/JSON",
        "synonyms": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/synonyms/JSON"
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
                print(f'{compound_name} generated an exception: {exc}')
                return None

    if all(key in results for key in ["properties", "synonyms"]):
        properties = results["properties"]['PropertyTable']['Properties'][0]
        synonyms_list = results["synonyms"]['InformationList']['Information'][0]['Synonym']
        # Ensure all synonyms for a compound have the same SMI_ number
        improve_drug_id = existing_synonyms.get(lower_compound_name, improve_drug_id)
        data_for_tsv = [{
            'improve_drug_id': improve_drug_id,
            'name': synonym.lower(),
            **properties
        } for synonym in synonyms_list if synonym.lower() not in existing_synonyms]
        return data_for_tsv
    else:
        return None

def fetch_data_for_batch(batch, improve_drug_id_start, existing_synonyms):
    all_data = []
    for compound_name in batch:
        data = retrieve_drug_info2(compound_name, f"SMI_{improve_drug_id_start}", existing_synonyms)
        if data:
            # Update the existing_synonyms with the new compound_name and its SMI_ number
            existing_synonyms[compound_name.lower()] = f"SMI_{improve_drug_id_start}"
            all_data.extend(data)
            # Increment improve_drug_id_start for the next unique compound
            improve_drug_id_start += 1
#         sleep(0.25)  # Prevent exceeding PubChem's rate limit
    return all_data, improve_drug_id_start

def read_existing_data(output_filename):
    try:
        df = pd.read_csv(output_filename, sep='\t')
        existing_synonyms = {row['chem_name'].lower(): "" for index, row in df.iterrows()}
        max_id = df['improve_drug_id'].str.extract(r'SMI_(\d+)').astype(float).max()
        if pd.isna(max_id[0]):
            improve_drug_id_start = 1
        else:
            improve_drug_id_start = int(max_id[0]) + 1
        return df, existing_synonyms, improve_drug_id_start
    except FileNotFoundError:
        return pd.DataFrame(), {}, 1
    
    
def update_dataframe_and_write_tsv(chem_name_file=None, unique_names=None, output_filename="drugs_new.tsv", batch_size=50):
    if bool(chem_name_file) == bool(unique_names):
        raise ValueError("You must provide either input_file or unique_names, but not both or neither.")
    processed_df, existing_synonyms, improve_drug_id_start = read_existing_data(output_filename)

    if chem_name_file:
        chem_name_file = open(chem_name_file, "r") 
        data = chem_name_file.read() 
        unique_names = data.split("\n") 
        chem_name_file.close() 
       
    unique_names = [name.lower() for name in unique_names if not pd.isna(name)]
    unique_names = list(set(unique_names) - set(existing_synonyms.keys()))

    print(f"Drugs to search: {len(unique_names)}")

    mode = 'w' if processed_df.empty else 'a'
    for i in range(0, len(unique_names), batch_size):
        batch = unique_names[i:i+batch_size]
        data, improve_drug_id_start = fetch_data_for_batch(batch, improve_drug_id_start, existing_synonyms)
        with open(output_filename, mode) as f:
            if mode == 'w':
                f.write("improve_drug_id\tchem_name\tpubchem_id\tcanSMILES\tisoSMILES\tInChIKey\tformula\tweight\n")
                mode = 'a'  # Change mode to 'append' after first write
            for entry in data:
                f.write(f"{entry['improve_drug_id']}\t{entry['name']}\t{entry['CID']}\t{entry['CanonicalSMILES']}\t{entry['IsomericSMILES']}\t{entry['InChIKey']}\t{entry['MolecularFormula']}\t{entry['MolecularWeight']}\n")

