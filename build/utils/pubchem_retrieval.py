import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import threading
import time
import signal
import sys

request_counter = 0
last_request_time = time.time()
lock = threading.Lock()
should_continue = True
improve_drug_id = 0
existing_synonyms = set()
existing_structures = dict()
existing_pubchemids = set()

def fetch_url(url):
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
    
    # Proceed with the request
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to fetch {url}")
    
def retrieve_drug_info(compound,ignore_chems,isname=True):
    '''
    compound_name: name of compound or CID
    ignore_chems: list of chemicals to ignore
    isname: true if compound is name, false if is CID
    '''
    global improve_drug_id, existing_synonyms, existing_structures
    if pd.isna(compound):
        return None

    if isname:
        urls = {
            "properties": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/property/CanonicalSMILES,InChIKey,MolecularFormula,MolecularWeight/JSON",
            "synonyms": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/synonyms/JSON"
        }
    else:
        urls = {
            "properties": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{compound}/property/CanonicalSMILES,InChIKey,MolecularFormula,MolecularWeight/JSON",
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
                with open(ignore_chems,"a") as f:
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
           # print(synonym)
            synonym_lower = str(synonym).lower()
#            if synonym_lower in existing_synonyms: ### THIS IS CAUSING THE LOOP TO END BEFORE IT GETS TO THE COMPOUND NAME
#                return None
            if synonym_lower not in existing_synonyms:
                new_syns.add(str(synonym).lower())
        if len(new_syns) == 0: #JUST BE SURE WE HAVE NO NEW SYNONYMS BEFORE RETURNING
            return None
        for synonym in new_syns:#synonyms_list + [compound]: ##NOW JUST ADD THOSE
            synonym_lower = str(synonym).lower()
            existing_synonyms.add(synonym_lower)
            

        ###now check for structure
        if properties['CanonicalSMILES'] in existing_structures.keys():
            print('found structure for '+str(compound))
            SMI_assignment = existing_structures[properties['CanonicalSMILES']]
        else:
            improve_drug_id += 1
            SMI_assignment = f"SMI_{improve_drug_id}"
            existing_structures[properties['CanonicalSMILES']] = SMI_assignment

        #print(new_syns)
        data_for_tsv = [{
            'improve_drug_id': SMI_assignment,
            'name': str(synonym).lower(),
            **properties
        } for synonym in new_syns]##synonyms_list]

        return data_for_tsv
    else:
        return None

def fetch_data_for_batch(batch,ignore_chems,isname):
    all_data = []
    for compound_name in batch:
        data = retrieve_drug_info(compound_name,ignore_chems,isname)
        if data:
            all_data.extend(data)
    return all_data

def read_existing_data(output_filename):
    global improve_drug_id,existing_synonyms,existing_structures
    try:
        df = pd.read_csv(output_filename, sep='\t', quoting=3)
        existing_synonyms = set([str(a).lower() for a in set(df.chem_name)])
        existing_pubchemids = set([str(a) for a in  df['pubchem_id']])
        max_id = df['improve_drug_id'].str.extract(r'SMI_(\d+)').astype(float).max()
        improve_drug_id = int(max_id[0]) + 1 if pd.notna(max_id[0]) else 1
        existing_structures = {row['canSMILES']:row['improve_drug_id'] for index,row in df.iterrows()}
        print('Read in '+str(len(existing_synonyms))+' drug names and '+str(len(existing_pubchemids))+' pubchem ids')
    except FileNotFoundError:
        return {}


def timeout_handler(signum, frame):
    global should_continue
    print("Time limit reached, exiting gracefully...")
    should_continue = False

# Call this function from other scripts. 
def update_dataframe_and_write_tsv(unique_names, output_filename="drugs.tsv",ignore_chems="ignore_chems.txt", batch_size=1,isname=True,time_limit=5*60*60):
    global should_continue, existing_synonyms, existing_pubchemids
    #time_limit=5*60*60 # 5 hours
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(time_limit)
    print('starting with '+str(len(unique_names))+' drug names/ids')
    try:
        print('reading in '+output_filename)
        read_existing_data(output_filename)
        if isname:
            unique_names = set([str(name).lower() for name in unique_names if not pd.isna(name)])
            unique_names = set(unique_names) - set(existing_synonyms)
            print('looking at '+str(len(unique_names))+' names')
        else:
            unique_names = set([str(name) for name in unique_names if not pd.isna(name)])
            unique_names = set(unique_names) - set(existing_pubchemids)
            print('looking at '+str(len(unique_names))+' ids')
        ignore_chem_set = set()
        if os.path.exists(ignore_chems):
            with open(ignore_chems, 'r') as file:
                for line in file:
                    ignore_chem_set.add(line.strip())
        unique_names = list(set(unique_names) - ignore_chem_set)
                
        print(f"{len(unique_names)} Drugs to search")#: {(unique_names)}")
        for i in range(0, len(unique_names), batch_size):
            if not should_continue:
                break
            if unique_names[i] in existing_synonyms or unique_names[i] in existing_pubchemids:
                continue
            
            batch = unique_names[i:i+batch_size]
            data = fetch_data_for_batch(batch,ignore_chems,isname)
            if data:
                file_exists = os.path.isfile(output_filename)
                mode = 'a' if file_exists else 'w' 
                with open(output_filename, mode) as f:
                    if not file_exists:
                        f.write("improve_drug_id\tchem_name\tpubchem_id\tcanSMILES\tInChIKey\tformula\tweight\n")
                    for entry in data:
                        f.write(f"{entry['improve_drug_id']}\t{entry['name']}\t{entry.get('CID', '')}\t{entry['CanonicalSMILES']}\t{entry['InChIKey']}\t{entry['MolecularFormula']}\t{entry['MolecularWeight']}\n")
                
                with open(ignore_chems,"a") as ig_f:
                    for entry in data:
                        if isname:
                            ig_f.write(f"{entry['name']}\n")
                        else:
                            ig_f.write(f"{entry.get('CID','')}\n")
                    
                    
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    finally:
        # Cancel the alarm
        signal.alarm(0)
