import pandas as pd
import requests
import os

def download_from_github(raw_url, save_path):
    """ 
    Download a file from a raw GitHub URL and save it to a local path.
    
    Parameters
    ----------
    raw_url : string
        The raw GitHub URL to download the file from.
        
    save_path : string
        Local path where the downloaded file will be saved.
        
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
    
    This function makes a POST request to the GDC API endpoint to fetch relevant metadata for the provided UUIDs.
    
    Parameters
    ----------
    uuids : list
        list of UUIDs

    Returns
    -------
    dict 
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
    filt= filt.drop_duplicates(subset='aliquot_id', keep=False)
    filt = filt.rename(
        columns={"tissue_or_organ_of_origin":"common_name",
                 "primary_diagnosis": "cancer_type",
                 "composition": "model_type",
                 "case_id": "other_names",
                 "aliquot_id": "other_id"}
            )
    filt = filt[["cancer_type","common_name","other_names","other_id","model_type"]]
    filt["other_id_source"] = "HCMI"
    # Create new improve sample IDs
    
    #Non-docker:
    # maxval = max(pd.read_csv('../cptac/cptac_samples.csv').improve_sample_id)
    # Docker:
    maxval = max(pd.read_csv('cptac_samples.csv').improve_sample_id)
    mapping = {other_id: i for i, other_id in enumerate(filt['other_id'].unique(), start=(int(maxval)+1))}
    # Use the map method to create the new column based on the lab-id column
    filt['improve_sample_id'] = filt['other_id'].map(mapping)
    return filt

def align_to_linkml_schema(input_df):
    """
    Maps the 'model_type' column of the input DataFrame to a set of predefined categories 
    according to a specified mapping dictionary. This alignment is intended to ensure 
    the DataFrame's 'model_type' values conform to a schema compatible with the LinkML model.
    
    Parameters
    ----------
    input_df : pd.DataFrame
        The input DataFrame containing a 'model_type' column with values to be mapped 
        according to the predefined categories.
    
    Returns
    -------
    pd.DataFrame
        A copy of the input DataFrame with the 'model_type' column values mapped to 
        a set of predefined categories ('tumor', 'organoid', 'cell line'). 
        The mapping is designed to align the DataFrame with the LinkML schema requirements.
    """
    
    mapping_dict = {
    'Solid Tissue': 'tumor',
    '3D Organoid': 'organoid',
    'Peripheral Blood Components NOS': 'tumor',
    'Buffy Coat': np.nan,
     None: np.nan,
    'Peripheral Whole Blood': 'tumor',
    'Adherent Cell Line': 'cell line',
    '3D Neurosphere': 'organoid',
    '2D Modified Conditionally Reprogrammed Cells': 'cell line',
    'Pleural Effusion': np.nan,
    'Human Original Cells': 'cell line',
    'Not Reported': np.nan, 
    'Mixed Adherent Suspension': 'cell line',
    'Cell': 'cell line',
    'Saliva': np.nan
    }

    # Apply mapping
    input_df['model_type'] = input_df['model_type'].map(mapping_dict)
    input_df.dropna(subset=['model_type'], inplace=True)
    
    return input_df

def main():
    """
    Retrieve and process HCMI (Human Cancer Models Initiative) samples metadata from GDC (Genomic Data Commons).
    Create samples.csv file for schema.

    This function automates the workflow of:
    1. Downloading a manifest file from the GitHub repository.
    2. Extracting UUIDs (Unique Universal Identifiers) from the manifest.
    3. Fetching the metadata for the samples corresponding to the UUIDs from GDC API via POST request.
    4. Structuring the fetched metadata into a pandas dataframe.
    5. Filtering and subsetting the dataframe to align with the schema.
    6. Writing the processed dataframe to a CSV file.

    Notes:
    ------
    The GDC API is publicly accessible, so no authentication is required.

    To Run:
    --------
    python createHCMISamplesFile.py

    Output:
    -------
    A local CSV file named 'samples.csv' containing the processed metadata.
    """
    
    manifest_path = "full_manifest.txt"
    #manifest_url = "https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/hcmi_update/hcmi/full_manifest.txt"
    #download_from_github(manifest_url, manifest_path)
    uuids = extract_uuids_from_manifest(manifest_path)
    metadata = fetch_metadata_for_samples(uuids)
    df = extract_data(metadata)
    output = filter_and_subset_data(df)
    aligned = align_to_linkml_schema(output)
    aligned.to_csv("hcmi_samples.csv",index=False)

 
main()


