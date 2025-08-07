import pandas as pd
import requests
import os
import argparse
import numpy as np


#this is table S1 - it has a mapping from patient number to organoid
sample_mapping='https://aacr.silverchair-cdn.com/aacr/content_public/journal/cancerdiscovery/8/9/10.1158_2159-8290.cd-18-0349/5/21598290cd180349-sup-199398_2_supp_4775186_p95dln.xlsx?Expires=1738004990&Signature=yngaaKNaXfIPCr-xLS2bDjX49n9py8JC7NBwi3q7m7ARYnK573eZwavFYmJOZVanL555vUWAr5x5k9b7IKj4VWHtZ-dts7BDzHd14AZh15LbsorJh-r3gjPliF7v1PIoAcGnEXjma2~kosmoDmyK0EDWXQCOE48tAaG5hFtaWAMMAINRMeBNgtDYk937Npc3Wb0IcGAdlgD2TJd8KJW2jQmcRspY1hfYssiS3BcWzuJrP-DVJeb-1V7-BnVNL6cVCkr7zHhau50H6aVgMVzk33F0gjCphl4r90OIx9UwE59hyNHbN9rFeeW26kDQpgCQKCj98Ol6CNQfLDsb2Zc5dQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'


def get_organoid_samples(sample_tab):
    '''
    takes as input a processed list of samples from HCMI and appends it with the 'organoid' identifier from the papers table S1 described above
    '''
    map = pd.read_excel(sample_mapping, sheet_name='Patient-Derived Organoid Cohort', skiprows=1)
    pmap = map[['Patient number','Organoid']]
    pmap = pmap.rename(columns={'Patient number':'common_name','Organoid':'experimentId'})
    
    #join with sampletab
    sample_tab.common_name=[str(a) for a in sample_tab.common_name]
    pmap.common_name = [str(a) for a in pmap.common_name]
    ocols=['common_name','other_names','model_type','cancer_type','improve_sample_id','species']
    red_tab = sample_tab[ocols].merge(pmap)
    
    #then add in organoid number
    newsamp = red_tab.melt(id_vars=ocols,value_vars='experimentId',var_name='other_id_source',value_name='other_id').drop_duplicates()
    res = pd.concat([sample_tab,newsamp])
    return res

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
    '3D Organoid': 'patient derived organoid',
    'Peripheral Blood Components NOS': 'tumor',
    'Buffy Coat': np.nan,
     None: np.nan,
    'Peripheral Whole Blood': 'tumor',
    'Adherent Cell Line': 'cell line',
    '3D Neurosphere': 'patient derived organoid',
    '2D Modified Conditionally Reprogrammed Cells': 'cell line',
    'Pleural Effusion': np.nan,
    'Human Original Cells': 'cell line',
    'Not Reported': np.nan, 
    'Mixed Adherent Suspension': 'cell line',
    'Cell': 'cell line',
    'Saliva': np.nan
    }

    # Apply mapping
    input_df['species'] = 'Homo sapiens (Human)' ##i assume they're lal human? 
    input_df['model_type'] = input_df['model_type'].map(mapping_dict)
    input_df.dropna(subset=['model_type'], inplace=True)
    input_df = input_df.sort_values(by='improve_sample_id')
    
    return input_df

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
            "cases.submitter_id,"
            "cases.annotations.case_submitter_id,"
            "cases.samples.sample_id,"
            "cases.samples.portions.analytes.aliquots.aliquot_id,"
            "cases.samples.sample_type,"
            "cases.diagnoses.submitter_id,"
            "cases.diagnoses.diagnosis_id,"
            "cases.diagnoses.classification_of_tumor,"
            "cases.diagnoses.tissue_or_organ_of_origin,"
            "cases.diagnoses.primary_diagnosis,"
            "cases.diagnoses.treatments.treatment_id,"##getting these but ignoring for now
            "cases.diagnoses.treatments.submitter_id," ##getting these but ignoring for now
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
                                    'entry_id': hit['id'],
                                    'case_uuid': case['case_id'],
                                    'case_id': case['submitter_id'],
                                    'tissue_or_organ_of_origin': diagnosis['tissue_or_organ_of_origin'],
                                    'primary_diagnosis': diagnosis['primary_diagnosis'],
                                    'diagnosis_id':diagnosis['submitter_id'],
                                    'tumor_classification':diagnosis['classification_of_tumor'],
                                    'sample_id': sample['sample_id'],
                                    'sample_type': sample['sample_type'],
                                    #'tumor_descriptor': sample.get('tumor_descriptor', None),
                                    'composition': sample.get('composition', None),
                                    'id': aliquot['aliquot_id']
                                })
    return pd.DataFrame(extracted)

def filter_and_subset_data(df, maxval, mapfile):
    """
    Filter and subset the data, then assign improve_sample_id at the end.

    Parameters
    ----------
    df : pd.DataFrame
        A tidied pandas DataFrame containing the full samples table.
    maxval : int
        The maximum value of improve_sample_id from previous samples, used to continue numbering.
    mapfile : str
        File path to the mapping file that maps primary diagnosis and tissue of origin to common cancer types.

    Returns
    -------
    pd.DataFrame
        The processed DataFrame ready for further use.
    """
    # Remove duplicates based on all columns except 'id'
    duplicates_mask = df.drop('id', axis=1).duplicated(keep='first')
    cmap = pd.read_csv(mapfile, encoding='ISO-8859-1')
    filt = df[~duplicates_mask]
    filt = filt.drop_duplicates()

    # Merge with the cancer type mapping file
    filt = pd.merge(
        filt,
        cmap,
        right_on=['tissue_or_organ_of_origin', 'primary_diagnosis'],
        left_on=['tissue_or_organ_of_origin', 'primary_diagnosis'],
        how='left'
    )

    # Rename columns to match the schema
    filt = filt.rename(
        columns={
            "composition": "model_type",
            "case_id": "common_name",
            "id": "other_names"
        }
    )

    # Melt the dataframe to create 'other_id' and 'other_id_source'
    longtab = pd.melt(
        filt,
        id_vars=['common_name', 'other_names', 'model_type', 'cancer_type'],
        value_vars=['diagnosis_id', 'tumor_classification', 'sample_type']
    )
    longtab = longtab.rename(columns={'variable': 'other_id_source', 'value': 'other_id'}).drop_duplicates()

    # Handle missing 'other_names'
    missing_other_names = longtab[longtab['other_names'].isnull()]
    if not missing_other_names.empty:
        print("Warning: Some samples have missing 'other_names' (aliquot_id). These samples will be excluded.")
        print(missing_other_names)
    longtab = longtab.dropna(subset=['other_names'])

    # Convert 'other_names' to string to ensure consistency
    longtab['other_names'] = longtab['other_names'].astype(str)

    #print(longtab)
    # Reassign 'improve_sample_id's at the end
    unique_other_names = longtab['other_names'].unique()
    print("Number of unique 'other_names' after filtering:", len(unique_other_names))

    ##UPDATE: assign them to common_names instead!
    unique_common_names = longtab['common_name'].unique()
    print("Number of unique 'common_names' after filtering:", len(unique_common_names))        
    # Create a new mapping
    #mapping = pd.DataFrame({
    #    'other_names': unique_other_names,
    #    'improve_sample_id': range(int(maxval) + 1, int(maxval) + len(unique_other_names) + 1)
        #})
    mapping = pd.DataFrame({
        'common_name':unique_common_names,
        'improve_sample_id': range(int(maxval) +1, int(maxval) + len(unique_common_names)+1)
        })

    # Merge the mapping back into 'longtab'
    longtab = pd.merge(longtab, mapping, on='common_name', how='left')

    # Debugging: Check longtab after reassigning IDs
    print("\nlongtab columns after reassigning 'improve_sample_id':", longtab.columns)
    print("longtab head after reassigning IDs:")
    print(longtab.head())

    # Verify that all 'improve_sample_id's are assigned
    missing_ids = longtab[longtab['improve_sample_id'].isnull()]
    if not missing_ids.empty:
        print("\nWarning: Some samples could not be assigned an 'improve_sample_id'.")
        print(missing_ids)
        
    # Missing cancer type indicates that it is normal tissue.
    longtab['cancer_type'] = longtab['cancer_type'].replace('', np.nan)
    longtab['cancer_type'] = longtab['cancer_type'].fillna('Normal Tissue')
    return longtab

def main():
    """
    Retrieve and process PANCPDO (Human Cancer Models Initiative) samples metadata from GDC (Genomic Data Commons).
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
    python createPANCPDOSamplesFile.py

    Output:
    -------
    A local CSV file named '/tmp/pancpdo_samples.csv' containing the processed metadata.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--prevSamples',dest='prev_samps', nargs='?',type=str, default='', const='', help='Previous sample file')
    parser.add_argument('--mapfile',dest='map',help='Mapping to common_cancer from primary_diagnosis and tissue_or_organ_of_origin',default='pancpdo_cancer_types.csv')

    args = parser.parse_args()
    manifest_path = "full_manifest.txt"
    #manifest_url = "https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/hcmi_update/pancpdo/full_manifest.txt"
    #download_from_github(manifest_url, manifest_path)
    uuids = extract_uuids_from_manifest(manifest_path)
    metadata = fetch_metadata_for_samples(uuids)
    df = extract_data(metadata)

    if args.prev_samps is None or args.prev_samps=='':
        print("No Previous Samples file was found. PANCPDO Data will not align with other datasets. Use ONLY for testing purposes.")
        maxval = 0
    else:
        print("Previous Samples File Provided. Running PANCPDO Sample File Generation")
        maxval = max(pd.read_csv(args.prev_samps).improve_sample_id)
    
    output = filter_and_subset_data(df,maxval,args.map)
    aligned = align_to_linkml_schema(output)
    print(aligned)
    aligned = get_organoid_samples(aligned)
    aligned.to_csv("/tmp/pancpdo_samples.csv",index=False)
 
main()


