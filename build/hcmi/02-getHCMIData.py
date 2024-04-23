import os
import subprocess
import pyarrow
import pandas as pd
from shutil import which, unpack_archive, rmtree
import shutil
import platform
import stat
import gzip
import requests
import wget
import argparse
import math
import time
import hashlib
import json
import polars as pl
import gc

    
def download_tool(url):
    """
    Download, extract, and make a tool (GDC Client) executable from the provided URL.

    Parameters
    ----------
    url : str
        The URL from where the tool needs to be downloaded.

    Returns
    -------
    str
        Name of the downloaded file.
    """
    
    filename = wget.download(url)
    files_before = os.listdir()
    shutil.unpack_archive(filename)
    files_after = os.listdir()
    new_file = str(next(iter((set(files_after) - set(files_before)))))
    st = os.stat(new_file)
    os.chmod(new_file, st.st_mode | stat.S_IEXEC)
    return filename

def is_tool(name):
    """
    Check if a specific tool is available on the system or in the current directory.

    Parameters
    ----------
    name : str
        The name of the tool to check.

    Returns
    -------
    bool
        True if the tool is found, otherwise False.
    """
    
    return which(name) is not None or name in os.listdir()

def ensure_gdc_client():
    """
    Ensure that the gdc-client is available on the system.
    
    If the gdc-client tool isn't found, this function will automatically
    download the appropriate version based on the operating system.
    """
    
    tool_name = "gdc-client"
    if not is_tool(tool_name):
        print("Downloading gdc-client")
        urls = {
            "Darwin": 'https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_OSX_x64.zip',
            "Windows": 'https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Windows_x64.zip',
            "Linux": 'https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip'
        }
        download_tool(urls.get(platform.system()))
    else:
        print("gdc-client already installed")

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

def fetch_metadata(uuids):
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
    payload = {
        "filters": {
            "op": "in",
            "content": {
                "field": "files.file_id",
                "value": uuids
            }
        },
        "fields": "cases.sample_ids,cases.case_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id",
        "format": "JSON",
        "size": str(len(uuids))
    }
    response = requests.post(endpoint, json=payload)
    return response.json()

def use_gdc_tool(manifest_data, data_type, download_data):
    """
    Use the gdc-client tool to download data based on the provided manifest and data type.
    
    The function first checks if the manifest location exists and cleans it if required. Then, it uses the gdc-client to download the data. Finally, it extracts UUIDs and fetches associated metadata.

    Parameters
    ----------
    manifest_data : str
        File path to the manifests file.
    data_type : str
        Type of the data to be downloaded.
    download_data : bool
        Flag to determine whether data should be downloaded or not.

    Returns
    -------
    dict
        Metadata associated with the UUIDs extracted from the manifest.
    """

    ##first, let's filter by type
    tdict={'transcriptomics':'rna_seq','copy_number':'copy_number','mutations':'ensemble_masked'}
    fm = pd.read_csv(manifest_data,sep='\t')
    fm['include'] = [tdict[data_type] in a for a in fm.filename]
    newfm = fm[fm.include]
#    newfm.reset_index(drop=True,inplace=True)
    newfm.to_csv('new_manifest.txt',sep='\t',index=False)
    
    manifest_loc = "full_manifest_files"

    if download_data:
        if os.path.isdir(manifest_loc):
            shutil.rmtree(manifest_loc)
        os.makedirs(manifest_loc)

        # Download files using gdc-client
        subprocess.run(['./gdc-client', 'download', '-d', manifest_loc, '-m','new_manifest.txt'])

    # Extract UUIDs and fetch metadata
    uuids = extract_uuids_from_manifest('new_manifest.txt')
    metadata = fetch_metadata(uuids)

    return metadata


def get_clean_files(data_type):
    """
    Extract clean files of a specified data type from manifest folders.
    
    Given a specific data type, this function looks through manifest folders to find 
    matching files and process them accordingly.
    
    Parameters
    ----------
    data_type : string
        The type of data being processed, e.g., "transcriptomics", "copy_number", or "mutations".
    
    Returns
    -------
    list of pl.DataFrame
        A list of polars dataframes containing cleaned data extracted from the manifest folders.
    """
    
   
    data_suffixes = {
        "transcriptomics": "rna_seq.augmented_star_gene_counts.tsv",
        "copy_number": "copy_number_variation.tsv",
        "mutations": "ensemble_masked.maf.gz"
    }

    suffix = data_suffixes.get(data_type)
    manifest = 'full_manifest_files'
    manifest_folders = [folder for folder in os.listdir(manifest) if folder != '.DS_Store']
    all_dataframes = []

    for folder_name in manifest_folders:
        folder_path = os.path.join(manifest, folder_name)
        folder_files = os.listdir(folder_path)

        sample_filenames = [x for x in folder_files if suffix in x and '.ipynb_checkpoints' not in x]

        for sample in sample_filenames:
            filepath = os.path.join(manifest, folder_name, sample)
            #gzipped data is mutation data
            if ".gz" in filepath:
                with gzip.open(filepath, 'rt') as f:
                    # Read into pandas DataFrame then convert. This is the only time pandas is used.
                    dataframe_pd = pd.read_csv(f, sep='\t', skiprows=7,low_memory=False)
                    dataframe = pl.DataFrame(dataframe_pd)
            else:
                if data_type == "transcriptomics":
                    dataframe = pl.read_csv(filepath, separator='\t',skip_rows=1)
                else: 
                    dataframe = pl.read_csv(filepath, separator='\t')

            dataframe = dataframe.with_columns(pl.lit(folder_name).alias('file_id'))

            if data_type == "transcriptomics":
                dataframe = dataframe[4:]
                if 'tpm_unstranded' in dataframe.columns:
                    new_columns = ['gene_id', 'gene_name', 'gene_type', 'tpm_unstranded', 'file_id']
                    dataframe = dataframe.select(new_columns)
                    dataframe = dataframe.filter(dataframe['gene_type'] == 'protein_coding')

            all_dataframes.append(dataframe)
        
    return all_dataframes

#old
def map_and_combine(dataframe_list, data_type, metadata, entrez_map_file):
    """
    Map and combine dataframes based on their data type, and merge with provided metadata
    using Polars.
    
    Parameters
    ----------
    dataframe_list : list of pl.DataFrame
        List of dataframes containing data to be mapped and combined.
        
    data_type : string
        The type of data being processed.
        
    metadata : dict
        Metadata to be merged with the processed data.
        
    entrez_map_file : string
        File path to the CSV file mapping gene names to Entrez IDs.
    
    Returns
    -------
    pl.DataFrame
        A dataframe containing the combined, mapped, and merged data.
    """
    
    # Initialize the list to hold mapped dataframes
    final_dataframe = pl.DataFrame()  # Initialize an empty DataFrame

    # Load mapping files using Polars
    genes = pl.read_csv(entrez_map_file)  # Map gene_name to entrez_id

    # Process each dataframe based on its data_type
    while dataframe_list:
        df = dataframe_list.pop()
        if data_type == "transcriptomics":
            mapped_df = df.join(genes, left_on='gene_name', right_on='gene_symbol', how='left')
            mapped_df = mapped_df.select(['gene_id', 'gene_name', 'gene_type', 'tpm_unstranded', 'entrez_id', 'file_id'])
            mapped_df = mapped_df.rename({'tpm_unstranded': 'transcriptomics'})
            mapped_df = mapped_df.with_columns([pl.lit('GDC').alias('source'),
                                               pl.lit('HCMI').alias('study')])
            
        elif data_type == "copy_number":
            joined_df = df.join(genes, left_on='gene_name', right_on='gene_symbol', how='left')
            selected_df = joined_df.select(['entrez_id', 'copy_number', 'file_id'])
            copy_call_series = copy_num(selected_df['copy_number'])
            mapped_df = selected_df.with_columns([
                copy_call_series.alias('copy_call'), 
                pl.lit('GDC').alias('source'),
                pl.lit('HCMI').alias('study')
            ])
            
        elif data_type == "mutations":
            mapped_df = df.rename({'Entrez_Gene_Id': 'entrez_id', 'HGVSc': 'mutation'})
            mapped_df = mapped_df.select(['entrez_id', 'mutation', 'Variant_Classification', 'file_id'])
            mapped_df = mapped_df.with_columns([pl.lit('GDC').alias('source'),
                                               pl.lit('HCMI').alias('study')])
            mapped_df = mapped_df.with_columns(mapped_df["entrez_id"].cast(str))

        final_dataframe = pl.concat([final_dataframe, mapped_df])
        del df, mapped_df
        gc.collect()

    
    # Convert the metadata into a DataFrame
    metadata_dict = {
    'file_id': [item['id'] for item in metadata['data']['hits']],
    'case_id': [item['cases'][0]['case_id'] for item in metadata['data']['hits']],
    'sample_id': [item['cases'][0]['samples'][0]['sample_id'] for item in metadata['data']['hits']],
    'aliquot_id': [item['cases'][0]['samples'][0]["portions"][0]["analytes"][0]['aliquots'][0]['aliquot_id'] for item in metadata['data']['hits']]
    }
    df_metadata = pl.DataFrame(metadata_dict)

#    print(final_dataframe)
#    print(df_metadata)
    # Merge the metadata DataFrame with the final dataframe based on 'file_id'
    final_dataframe = final_dataframe.join(df_metadata, on='file_id', how='left')
    
    return final_dataframe



def retrieve_figshare_data(url):
    """
    Download data from a given Figshare URL.
    
    Parameters
    ----------
    url : string
        The Figshare URL to download data from.
    
    Returns
    -------
    string
        Name of the downloaded file.
    """
    
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    new_file = str(next(iter(set(files_1) - set(files_0))))
    return new_file

def copy_num(arr):
    """
    Determine copy number variations for a given array of values.
    
    The function maps numerical copy number values to their respective categorical descriptions.
    
    Parameters
    ----------
    arr : list of float
        List of copy number values.
    
    Returns
    -------
    list of string
        List of copy number descriptions corresponding to the input values.
    """
        
    def get_copy_call(a):
        """
        Helper Function - Determine copy call for a value.
        """

        if a is None:
            return float('nan')

        if math.isnan(a):
            return float('nan')
        
        a_val = math.log2(float(a)+0.000001) ###this should not be exponent, should be log!!! 2**float(a)
        if a_val < 0.0: #0.5210507:
            return 'deep del'
        elif a_val < 0.7311832:
            return 'het loss'
        elif a_val < 1.214125:
            return 'diploid'
        elif a_val < 1.731183:
            return 'gain'
        else:
            return 'amp'

    return pl.Series([get_copy_call(a) for a in arr])


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


def align_to_schema(data, data_type, chunksize=7500,samples_path='/tmp/hcmi_samples.csv'):
    """
    Modify the data match the CANDLE schema based on its type, using Polars for processing.
    Essentially just adding improve_sample_id
    
    Parameters
    ----------
    data : pl.DataFrame
        The data to be aligned.
        
    data_type : string
        The type of data being processed.
    
    Returns
    -------
    pl.DataFrame
        The final form of the dataframe.
    """
#    samples_path = "/tmp/hcmi_samples.csv"
    samples = pl.read_csv(samples_path)
    samples = samples.drop(["cancer_type", "common_name", "other_names", "model_type", "other_id_source"])

    # Determine columns to select based on data_type
    columns = {
        "transcriptomics": ["entrez_id", "transcriptomics", "source", "study", "aliquot_id"],
        "copy_number": ["entrez_id", "copy_number", "copy_call", "source", "study", "aliquot_id"],
        "mutations": ["entrez_id", "mutation", "variant_classification", "source", "study", "aliquot_id"]
    }
    selected_columns = columns.get(data_type, [])

    # Process in chunks
    merged_data = pl.DataFrame()
    for i in range(0, len(data), chunksize):
        chunk = data[i:i + chunksize]
        if data_type == "mutations":
            chunk = chunk.rename({"Variant_Classification": "variant_classification"})
        chunk = chunk.select(selected_columns)
        
        merged_chunk = samples.join(chunk, left_on='other_id', right_on='aliquot_id', how='inner')
        merged_chunk = merged_chunk.drop(["aliquot_id", "other_id"])

        # Append the processed chunk
        merged_data = pl.concat([merged_data, merged_chunk])
        gc.collect()

    return merged_data


def write_dataframe_to_csv(dataframe, outname):
    """
    Wrapper function for pandas "to_csv" function. 
    Write the provided dataframe to a CSV file.
    
    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataframe to be written to a CSV file.
        
    outname : string
        The desired name for the output CSV file.
    
    Returns
    -------
    None
    """
    if('gz' in outname):
        dataframe.to_pandas().to_csv(outname,compression='gzip')
    else:
        dataframe.to_pandas().to_csv(outname)
    return

def upload_to_figshare(token, title, filepath):
    """
    Uploads a file to Figshare and publishes the article.

    This function automates the process of uploading a file to Figshare and
    subsequently publishing the associated article.

    Parameters
    ----------
    token : str
        The authentication token for Figshare API access.
    
    title : str
        The title of the article to be created on Figshare.
    
    filepath : str
        The path to the file that is to be uploaded to Figshare.

    Notes
    -----
    The function uses various helper functions for each specific tasks.
    - Sending requests to the Figshare API.
    - Creating a new article on Figshare.
    - Computing the MD5 checksum and size of a given file.
    - Handling the multipart upload for large files.
    - Publishing the article once the upload is complete.

    The Figshare API endpoint is defined by the `BASE_URL` constant and the size
    of each chunk in the multipart upload is defined by the `CHUNK_SIZE` constant.

    Raises
    ------
    HTTPError
        If there is any issue with the Figshare API requests.
    
    ValueError
        If there's an issue parsing the response from the Figshare API.

    """
    
    BASE_URL = 'https://api.figshare.com/v2/{endpoint}'
    CHUNK_SIZE = 1048576
    def raw_issue_request(method, url, data=None, binary=False):
        """
        Sends an HTTP request and returns the response.
        """
        headers = {'Authorization': 'token ' + token}
        if data is not None and not binary:
            data = json.dumps(data)
        response = requests.request(method, url, headers=headers, data=data)
        response.raise_for_status()
        try:
            data = json.loads(response.content)
        except ValueError:
            data = response.content
        return data

    def issue_request(method, endpoint, *args, **kwargs):
        """
        Sends a request to a specific Figshare API endpoint.
        """
        return raw_issue_request(method, BASE_URL.format(endpoint=endpoint), *args, **kwargs)

    def create_article(title):
        """
        Creates a new article on Figshare with a given title.
        """
        data = {
            'title': title,
            'description': "Cancer Organoid Data",
            'keywords': ["cancer","organoid","hcmi"],
            "categories_by_source_id": [
            "321101",
            "400207"
          ]
        }
        result = issue_request('POST', 'account/articles', data=data)
        result = raw_issue_request('GET', result['location'])
        return result['id']
    
    

    def get_file_check_data(file_name):
        """
        Check the MD5 checksum and size of a given file.
        """
        with open(file_name, 'rb') as fin:
            md5 = hashlib.md5()
            size = 0
            data = fin.read(CHUNK_SIZE)
            while data:
                size += len(data)
                md5.update(data)
                data = fin.read(CHUNK_SIZE)
            return md5.hexdigest(), size

    def initiate_new_upload(article_id, file_name):
        """
        Initiates the file upload process for a specific article.
        """
        endpoint = 'account/articles/{}/files'.format(article_id)
        md5, size = get_file_check_data(file_name)
        data = {'name': os.path.basename(file_name),
                'md5': md5,
                'size': size}
        result = issue_request('POST', endpoint, data=data)
        result = raw_issue_request('GET', result['location'])
        return result

    def complete_upload(article_id, file_id):
        """
        Marks the file upload as complete for a specific article.
        """
        issue_request('POST', 'account/articles/{}/files/{}'.format(article_id, file_id))

    def upload_parts(file_info):
        """
        Handles the multipart upload for large files.
        """
        url = '{upload_url}'.format(**file_info)
        print("url: ", url)
        result = raw_issue_request('GET', url)
        print("result: ", result)
        with open(filepath, 'rb') as fin:
            for part in result['parts']:
                upload_part(file_info, fin, part)

    def upload_part(file_info, stream, part):
        """
        Uploads a specific part/chunk of the file.
        """
        udata = file_info.copy()
        udata.update(part)
        url = '{upload_url}/{partNo}'.format(**udata)
        stream.seek(part['startOffset'])
        data = stream.read(part['endOffset'] - part['startOffset'] + 1)
        raw_issue_request('PUT', url, data=data, binary=True)
        

    def change_article_status(article_id, status='public'):
        """
        Changes the visibility status of a specific article on Figshare.
        """
        data = {'status': status}
        response = issue_request('PUT', f'account/articles/{article_id}', data=data)
        
    def publish_article(token, article_id):
        """
        Publishes a specific article on Figshare.
        """
        headers = {
            'Authorization': 'token ' + token,
            'Content-Type': 'application/json'
        }
    
        url = BASE_URL.format(endpoint=f'account/articles/{article_id}/publish')
        response = requests.post(url, headers=headers)

        # Handle the response
        if response.status_code == 201:
            print("Article published successfully!")
        else:
            print("Error:", response.status_code)
            print(response.text)
            
    article_id = create_article(title)
    file_info = initiate_new_upload(article_id, filepath)
    upload_parts(file_info)
    complete_upload(article_id, file_info['id'])
    change_article_status(article_id, 'public')
    print("Make the file public")
    publish_article(token,article_id)


def main():
    """
    Automates the process of retrieving and processing HCMI (Human Cancer Models Initiative)
    data from Genomic Data Commons Data Portal.

    This function handles the entire workflow of the data processing, including:
    - Parsing command-line arguments.
    - Downloading and processing the data based on the manifest file or folder.
    - Using the GDC-client tool to retrieve metadata.
    - Extracting relevant data files based on the manifest.
    - Retrieving mapping gene/samples data from Figshare/Github.
    - Combining and mapping data.
    - Formatting the combined data the schema.
    - Writing the final dataset to a CSV file.
    - Uploading and Publishing the data on figshare.

    Parameters (Command-line Arguments)
    -----------------------------------
    -m, --manifest : str
        Path to the manifest file to specify which datasets to download and process.
    
    -M, --manifestfolder : str, optional
        Path to a folder containing manifest files. If provided, data download is skipped.
    
    -t, --type : {'transcriptomics', 'copy_number', 'mutations'}
        Type of data to process. This determines how the data is parsed and combined.
    
    -o, --outname : str
        Name of the output CSV file where the processed data will be saved.
    
    -z, --token : str
        Authentication token for accessing private Figshare datasets.
    -g, --genes : str
        File containing list of genes

    Example Usage
    -------------
    - Without manifest files already downloaded:
        python getHCMIData.py -m path_to_manifest_file -t transcriptomics -o transcriptomics.csv
    
    - With manifest files already downloaded:
        python getHCMIData.py -M path_to_manifest_folder -t copy_number -o copy_number.csv
    
    Notes
    -----
    The function checks if a manifest folder is provided. If so, it skips the data download step
    and proceeds with the provided data. Otherwise, it ensures the GDC client is present and 
    downloads the data using the provided manifest file.
    
    This will publish the data to figshare. Don't include token if this is not desired.
    
    Raises
    ------
    ValueError, HTTPError, and other exceptions based on underlying functions.
    
    """
    
    parser = argparse.ArgumentParser(description='Process data.')
    parser.add_argument('-m', '--manifest', help='Path to manifest file', required=True)
    parser.add_argument('-M', '--manifestfolder', help='Path to manifest folder', required=False)
    tc = ['transcriptomics', 'copy_number', 'mutations']
    parser.add_argument('-t', '--type', help='Type of data (e.g., transcriptomics, copy_number)',choices = tc, required=True)
    parser.add_argument('-o', '--outname', help='Output CSV Name', required=True)
    parser.add_argument('-z', '--token', help='figshare token ID', required=False)
    parser.add_argument('-s', '--samples',help='Samples file', required=False,default='/tmp/hcmi_samples.csv')
    parser.add_argument('-g', '--genes',help='File containing valid gene ids',required=False,default='/tmp/genes.csv')
    args = parser.parse_args()
    
        
    if bool(args.manifestfolder):
        download_option = False
        print("Using provided manifest folder without downloading data...")
    else:
        download_option = True
        ensure_gdc_client()
        print("Using provided manifest and downloading data...")
        
          
    # Use gdc tool to get metadata
    print("Using gdc tool and retrieving get metadata...")
    metadata = use_gdc_tool(args.manifest, args.type, download_data=download_option)
    #print(metadata)
    # Extract data files
    print("Running 'get_clean_files' function")
    data_files = get_clean_files(args.type)

    # Retrieve figshare gene data for entrez map
    print("Running 'retrieve_figshare_data' function")
    gene_url = "https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47"
    entrez_map_file = args.genes #retrieve_figshare_data(gene_url)
    gc.collect()
    
    # Combine the data
    print("Running 'map_and_combine' function")
    combined_data = map_and_combine(data_files, args.type, metadata, entrez_map_file)
    gc.collect()
    data_files = None
    metadata = None
    
    # Final formatting
    print("Aligning to Schema")
    final_data = align_to_schema(combined_data,args.type,7500,args.samples)
    gc.collect()
    combined_data = None
    
    print(f"final data:\n{final_data}")

    # Save to CSV
    print("Running 'write_dataframe_to_csv' function")
    write_dataframe_to_csv(final_data, args.outname)
    
    print("Data processing complete!")
    if args.token:
        token = args.token
        print("Running 'upload_to_figshare' function")
        upload_to_figshare(token, args.outname, args.outname)
    else:
        print("No token provided. Data not uploaded to Figshare")


if __name__ == "__main__":
    main()
    
    
