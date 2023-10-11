import pandas as pd
import requests

def download_from_github(raw_url, save_path):
    """ 
    Download a file from github.
    
    This will use requests to pull a file from github and save it locally.
    
    Parameters
    ----------
    raw_url : string
        github url to download from
        
    save_path : string
        path of location to save file to
        
    Returns
    -------
    None
    """
    response = requests.get(raw_url)
    with open(save_path, 'wb') as f:
        f.write(response.content)
    return

def extract_uuids_from_manifest(manifest_data):
    """
    Extract UUIDs from the provided manifest data. 
    
    Takes a manifests file generated from GDC portal (or manually) and parses through while collecting UUIDs.
    
    Parameters
    ----------
    manifest_data : string
        file path to manifests file

    Returns
    -------
    List of UUIDs
    """
    with open(manifest_data, 'r') as f:
        lines = f.readlines()[1:]  # Skip header
        return [line.split("\t")[0] for line in lines]
    
    
def fetch_metadata_for_samples(uuids):
    """
    Fetch metadata for given UUIDs.
    
    Takes a list of UUIDs and fetches metadata from the GDC Portal
    
    Parameters
    ----------
    uuids : list
        list of UUIDs

    Returns
    -------
    JSON Request Data
    """
    
    endpoint = "https://api.gdc.cancer.gov/files"
    
    filters_content = {
        "field": "files.file_id",
        "value": uuids
    }
    
    payload = {
        "filters": {
            "op": "in",
            "content": filters_content
        },
        "fields": (
            "cases.sample_ids,"
            "cases.case_id,"
            "cases.samples.sample_id,"
            "cases.samples.portions.analytes.aliquots.aliquot_id,"
            "cases.samples.sample_type,"
            "cases.diagnoses.tissue_or_organ_of_origin,"
            "cases.diagnoses.primary_diagnosis,"
            "cases.samples.tumor_descriptor,"
            "cases.samples.composition"
        ),
        "format": "JSON",
        "size": str(len(uuids))
    }
    
    response = requests.post(endpoint, json=payload)
    return response.json()


def extract_data(data):
    """
    Write API returned JSON Data to Pandas Table
        
    Parameters
    ----------
    data : json data
        json data from GDC Portal

    Returns
    -------
    Pandas Dataframe
    """
    extracted = []
    for hit in data['data']['hits']:
        for case in hit['cases']:
            for idx, sample in enumerate(case['samples']):
                for portion in sample['portions']:
                    for analyte in portion['analytes']:
                        for aliquot in analyte['aliquots']:
                            if idx < len(case['diagnoses']):
                                diagnosis = case['diagnoses'][idx]
                                extracted.append({
                                    'id': hit['id'],
                                    'case_id': case['case_id'],
                                    'tissue_or_organ_of_origin': diagnosis['tissue_or_organ_of_origin'],
                                    'primary_diagnosis': diagnosis['primary_diagnosis'],
                                    'sample_id': sample['sample_id'],
                                    'sample_type': sample['sample_type'],
                                    'tumor_descriptor': sample.get('tumor_descriptor', None),
                                    'composition': sample.get('composition', None),
                                    'aliquot_id': aliquot['aliquot_id']
                                })
    return pd.DataFrame(extracted)

def filter_and_subset_data(df):
    """
    Filter and subset the data.
    
    Taking a pandas dataframe containing all sample information, filter it to desired columns and rename them to match schema.
        
    Parameters
    ----------
    df : pandas dataframe
        full samples table

    Returns
    -------
    Pandas Dataframe
    """
    duplicates_mask = df.drop('id', axis=1).duplicated(keep='first')
    filt = df[~duplicates_mask]
    filt = filt[(filt.tumor_descriptor != "Not Applicable")]
    filt = filt.rename(
        columns={"tissue_or_organ_of_origin":"common_name",
                 "primary_diagnosis": "cancer_type",
                 "composition": "model_type",
                 "case_id": "other_names",
                 "aliquot_id": "other_id"}
            )
    filt = filt[["cancer_type","common_name","other_names","other_id","model_type"]]
    filt["other_id_source"] = "HCMI"
    return filt


def main():
    manifest_path = "full_manifest.txt"
    manifest_url = "https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/hcmi_update/hcmi/full_manifest.txt"
    download_from_github(manifest_url, manifest_path)
    uuids = extract_uuids_from_manifest(manifest_path)
    metadata = fetch_metadata_for_samples(uuids)
    df = extract_data(metadata)
    output = filter_and_subset_data(df)
    output.to_csv("samples_new_version.csv",index=False)

    
main()