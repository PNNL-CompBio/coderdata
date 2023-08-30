# from __future__ import print_statement
import os
import subprocess
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
# import time
# import swagger_client
# from swagger_client.rest import ApiException
# from pprint import pprint

    
def download_tool(url):
    """Download and extract tool from a given URL, and then make it executable."""
    filename = wget.download(url)
    files_before = os.listdir()
    shutil.unpack_archive(filename)
    files_after = os.listdir()
    new_file = str(next(iter((set(files_after) - set(files_before)))))
    st = os.stat(new_file)
    os.chmod(new_file, st.st_mode | stat.S_IEXEC)
    return filename

def is_tool(name):
    """Check whether the given tool is available on the system or in the directory."""
    return which(name) is not None or name in os.listdir()


def ensure_gdc_client():
    """Ensure the gdc-client is available, if not, download it."""
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
    """Extract UUIDs from the provided manifest data."""
    with open(manifest_data, 'r') as f:
        lines = f.readlines()[1:]  # Skip header
        return [line.split("\t")[0] for line in lines]


def fetch_metadata(uuids):
    """Fetch metadata for given UUIDs."""
    endpoint = "https://api.gdc.cancer.gov/files"
    payload = {
        "filters": {
            "op": "in",
            "content": {
                "field": "files.file_id",
                "value": uuids
            }
        },
        "fields": "cases.sample_ids,cases.case_id,cases.samples.sample_id",
        "format": "JSON",
        "size": str(len(uuids))
    }
    response = requests.post(endpoint, json=payload)
    return response.json()


def use_gdc_tool(manifest_data, data_type, download_data):
    """Use gdc-client tool to download data for given manifest and data type."""
    manifest_loc = data_type + "_manifest_files"

    if download_data:
        if os.path.isdir(manifest_loc):
            shutil.rmtree(manifest_loc)
        os.makedirs(manifest_loc)

        # Download files using gdc-client
        subprocess.run(['./gdc-client', 'download', '-d', manifest_loc, '-m', manifest_data])

    # Extract UUIDs and fetch metadata
    uuids = extract_uuids_from_manifest(manifest_data)
    metadata = fetch_metadata(uuids)

    return metadata


def get_clean_files(data_type):
    
    data_suffixes = {
        "transcriptomics": "rna_seq.augmented_star_gene_counts.tsv",
        "copy_number": "copy_number_variation.tsv",
        "mutations": "ensemble_masked.maf.gz"
    }

    suffix = data_suffixes.get(data_type)
    
    manifest = f'{data_type}_manifest_files'
    
    # Filter out unwanted folders like .DS_Store
    manifest_folders = [folder for folder in os.listdir(manifest) if folder != '.DS_Store']
    
    # Extract filenames that match the given suffix, exclude any with '.ipynb_checkpoints'
    sample_filenames = [
        x for folder_name in manifest_folders 
        for x in os.listdir(os.path.join(manifest, folder_name))
        if suffix in x and '.ipynb_checkpoints' not in x
    ]

    all_dataframes = []
    first = True

    # Create file paths, read files into dataframes and process them
    for folder_n, sample in zip(manifest_folders, sample_filenames):
        filepath = os.path.join(manifest, folder_n, sample)
        if ".gz" in filepath:
            with gzip.open(filepath) as f:
                dataframe = pd.read_csv(f, delimiter='\t', skiprows=7)
        else:
            dataframe = pd.read_csv(filepath, delimiter='\t')

        dataframe['file_id'] = folder_n
        dataframe.reset_index(inplace=True)

        # Special handling for transcriptomics data type
        if data_type == "transcriptomics":
            dataframe.columns = dataframe.iloc[0]
            if 'tpm_unstranded' in dataframe.columns:
                dataframe = dataframe[5:]
                new_index = ['gene_id', 'gene_name', 'gene_type', 'tpm_unstranded', 'file_id']
                dataframe = dataframe.reindex(columns=new_index)
                dataframe['file_id'] = folder_n
                
        # Print the first dataframe for sanity
        if first:
            print(dataframe)
            first = False

        all_dataframes.append(dataframe)
        
    return all_dataframes


def retrieve_figshare_data(url):
    """
    Retrieve data from a given figshare URL.
    
    Return: File name (genes.csv)
    """
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    new_file = str(next(iter(set(files_1) - set(files_0))))
    return new_file

def copy_num(arr):
    """Determine copy number variations."""
    def get_copy_call(a):
        if math.isnan(a):
            return float('nan')
        
        a_val = 2**float(a)
        if a_val < 0.5210507:
            return 'deep del'
        elif a_val < 0.7311832:
            return 'het loss'
        elif a_val < 1.214125:
            return 'normal'
        elif a_val < 1.731183:
            return 'low gain'
        else:
            return 'high gain'

    return [get_copy_call(a) for a in arr]



def map_and_combine(dataframe_list, data_type, metadata, entrez_map_file):
    """
    Maps and combines dataframes based on their data_type. It then merges 
    the final dataframe with provided metadata.
    """
    
    # Initialize the list to hold mapped dataframes
    df_list = []

    # Load mapping files
    genes = pd.read_csv(entrez_map_file)          # Map gene_name to entrez_id

    # Process each dataframe based on its data_type
    for df in dataframe_list:
        if data_type == "transcriptomics":
            mapped_df = df.merge(genes, left_on='gene_name', right_on='gene_symbol', how='left').reindex(
                            columns=['gene_id', 'gene_name', 'gene_type', 'tpm_unstranded', 'entrez_id', 'file_id'])
            mapped_df = mapped_df.rename(columns={'tpm_unstranded': 'transcriptomics'})
            mapped_df['source'] = 'GDC'
            mapped_df['study'] = 'HCMI'

        elif data_type == "copy_number":
            mapped_df = df.merge(genes, left_on='gene_name', right_on='gene_symbol', how='left').reindex(
                            columns=['entrez_id', 'copy_number', 'file_id'])
            mapped_df['copy_call'] = copy_num(mapped_df['copy_number'].tolist())  # Assuming copy_num is a defined function
            mapped_df['source'] = 'GDC'
            mapped_df['study'] = 'HCMI'

        elif data_type == "mutations":
            mapped_df = df.reindex(columns=['Entrez_Gene_Id', 'file_id', 'HGVSc', "Variant_Classification"])
            mapped_df.rename(columns={"HGVSc": "mutations"}, inplace=True)
            mapped_df['source'] = 'GDC'
            mapped_df['study'] = 'HCMI'
        
        df_list.append(mapped_df)

    # Concatenate the list of dataframes into one final dataframe
    final_dataframe = pd.concat(df_list)
    
    # Convert the metadata into a DataFrame
    metadata_list = [[item['id'], item['cases'][0]['case_id'], item['cases'][0]['samples'][0]['sample_id']]
                     for item in metadata['data']['hits']]
    df_metadata = pd.DataFrame(metadata_list, columns=['file_id', 'case_id', 'sample_id'])
    
    # Convert 'file_id' columns to string type for accurate merging
    final_dataframe['file_id'] = final_dataframe['file_id'].astype(str)
    df_metadata['file_id'] = df_metadata['file_id'].astype(str)

    # Merge the metadata DataFrame with the final dataframe based on 'file_id'
    final_dataframe = pd.merge(final_dataframe, df_metadata, on='file_id', how='left')
    
    return final_dataframe


def write_dataframe_to_csv(dataframe, outname):
    """
    Writes the given dataframe to a CSV file with the specified filename.
    Args:
    - dataframe (pd.DataFrame): The DataFrame to be saved.
    - outname (str): The name of the output CSV file.
    """
    dataframe.to_csv(outname, index=False)

    
    
    
# BASE_URL = "https://api.figshare.com/v2"

# def upload_to_figshare(file_path, title, token):
#     """
#     Upload a file to Figshare and make the article public.

#     Parameters:
#     - file_path (str): Path to the file to be uploaded.
#     - title (str): Title of the article on Figshare.
#     - token (str):  Personal access token from Figshare.

#     Returns:
#     - str: URL of the published article.
#     """

#     headers = {
#         "Authorization": f"token {token}",
#         "Content-Type": "application/json"
#     }

#     # 1. Create a new article (draft)
#     article_data = {
#         "title": title,
#         "defined_type": "dataset"
#     }

#     response = requests.post(f"{BASE_URL}/account/articles", headers=headers, json=article_data)
#     response.raise_for_status()  # Raise an exception if the request was unsuccessful
#     location = response.json()["location"]
#     print(response.json())
    
    
#     # Configure OAuth2 access token for authorization: OAuth2
#     swagger_client.configuration.access_token = 'YOUR_ACCESS_TOKEN'

#     # create an instance of the API class
#     api_instance = swagger_client.ArticlesApi()
#     articleId = 789 # Long | Article unique identifier
#     file =  # FileCreator | 

#     try: 
#         # Initiate Upload
#         api_response = api_instance.private_article_upload_initiate(articleId, file)
#         pprint(api_response)
#     except ApiException as e:
#         print("Exception when calling ArticlesApi->privateArticleUploadInitiate: %s\n" % e)

      
#     # 2. Upload the file
#     print("Upload the file")
#     with open(file_path, 'rb') as f:
#         upload_response = requests.put(f"{location}/files", headers=headers, data=f)
#         print(upload_response.text)
#         upload_response.raise_for_status()

#     print("Make the file public")
#     # 3. Make the article public
#     publish_response = requests.post(f"{location}/publish", headers=headers)
#     publish_response.raise_for_status()

    # Return the URL of the published article
#     return location
  


# def raw_issue_request(method, url, data=None, binary=False):
#     headers = {'Authorization': 'token ' + token}
#     if data is not None and not binary:
#         data = json.dumps(data)
#     response = requests.request(method, url, headers=headers, data=data)
#     try:
#         response.raise_for_status()
#         try:
#             data = json.loads(response.content)
#         except ValueError:
#             data = response.content
#     except HTTPError as error:
#         print('Caught an HTTPError: ', error.message)
#         print('Body:\n', response.content)
#         raise

#     return data


# def issue_request(method, endpoint, *args, **kwargs):
#     return raw_issue_request(method, BASE_URL.format(endpoint=endpoint), *args, **kwargs)


# def list_articles():
#     result = issue_request('GET', 'account/articles')
#     print('Listing current articles:')
#     if result:
#         for item in result:
#             print(**item)
#     else:
#         print('No articles.')


# def create_article(title):
#     data = {
#         'title': title  # You may add any other information about the article here as you wish.
#     }
#     result = issue_request('POST', 'account/articles', data=data)
#     print('Created article:', result['location'], '\n')

#     result = raw_issue_request('GET', result['location'])

#     return result['id']


# def list_files_of_article(article_id):
#     result = issue_request('GET', 'account/articles/{}/files'.format(article_id))
#     print('Listing files for article {}:',article_id)
#     if result:
#         for item in result:
#             print('  {id} - {name}',**item)
#     else:
#         print('  No files.')


# def get_file_check_data(file_name):
#     with open(file_name, 'rb') as fin:
#         md5 = hashlib.md5()
#         size = 0
#         data = fin.read(CHUNK_SIZE)
#         while data:
#             size += len(data)
#             md5.update(data)
#             data = fin.read(CHUNK_SIZE)
#         return md5.hexdigest(), size


# def initiate_new_upload(article_id, file_name):
#     endpoint = 'account/articles/{}/files'
#     endpoint = endpoint.format(article_id)

#     md5, size = get_file_check_data(file_name)
#     data = {'name': os.path.basename(file_name),
#             'md5': md5,
#             'size': size}

#     result = issue_request('POST', endpoint, data=data)
#     print('Initiated file upload:', result['location'], '\n')

#     result = raw_issue_request('GET', result['location'])

#     return result


# def complete_upload(article_id, file_id):
#     issue_request('POST', 'account/articles/{}/files/{}'.format(article_id, file_id))


# def upload_parts(file_info):
#     url = '{upload_url}'.format(**file_info)
#     result = raw_issue_request('GET', url)

#     print('Uploading parts:')
#     with open(FILE_PATH, 'rb') as fin:
#         for part in result['parts']:
#             upload_part(file_info, fin, part)


# def upload_part(file_info, stream, part):
#     udata = file_info.copy()
#     udata.update(part)
#     url = '{upload_url}/{partNo}'.format(**udata)

#     stream.seek(part['startOffset'])
#     data = stream.read(part['endOffset'] - part['startOffset'] + 1)

#     raw_issue_request('PUT', url, data=data, binary=True)
#     print('  Uploaded part {partNo} from {startOffset} to {endOffset}',**part)


# def push_to_figshare(filepath, title, token):
#     BASE_URL = 'https://api.figshare.com/v2/{endpoint}'
#     CHUNK_SIZE = 1048576
#     # We first create the article
#     print("a")
#     list_articles()
#     print("b")
#     article_id = create_article(title)
#     print("c")
#     list_articles()
#     print("d")
#     list_files_of_article(article_id)
#     print("e")
#     # Then we upload the file.
#     file_info = initiate_new_upload(article_id, filepath)
#     print("f")
#     # Until here we used the figshare API; following lines use the figshare upload service API.
#     upload_parts(file_info)
#     print("g")
#     # We return to the figshare API to complete the file upload process.
#     complete_upload(article_id, file_info['id'])
#     print("h")
#     list_files_of_article(article_id)
#     print("i")

    
   



def upload_to_figshare(token, title, filepath):
    BASE_URL = 'https://api.figshare.com/v2/{endpoint}'
    CHUNK_SIZE = 1048576
    def raw_issue_request(method, url, data=None, binary=False):
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
        return raw_issue_request(method, BASE_URL.format(endpoint=endpoint), *args, **kwargs)

    def create_article(title):
#         data = {'title': title}
#         data = {
#             'title': title,
#             'description': "Cancer Organoid Data",
#             'keywords': ["cancer","organoid","hcmi"],
#             'categories': [11]
#         }
        data = {
            'title': title,
            'description': "Cancer Organoid Data",
            'keywords': ["cancer","organoid","hcmi"],
#             'categories': [1],
            "categories_by_source_id": [
            "321101",
            "400207"
          ]
        }
        result = issue_request('POST', 'account/articles', data=data)
        result = raw_issue_request('GET', result['location'])
        return result['id']
    
    

    def get_file_check_data(file_name):
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
        endpoint = 'account/articles/{}/files'.format(article_id)
        md5, size = get_file_check_data(file_name)
        data = {'name': os.path.basename(file_name),
                'md5': md5,
                'size': size}
        result = issue_request('POST', endpoint, data=data)
        result = raw_issue_request('GET', result['location'])
        return result

    def complete_upload(article_id, file_id):
        issue_request('POST', 'account/articles/{}/files/{}'.format(article_id, file_id))

    def upload_parts(file_info):
        url = '{upload_url}'.format(**file_info)
        print("url: ", url)
        result = raw_issue_request('GET', url)
        print("result: ", result)
        with open(filepath, 'rb') as fin:
            for part in result['parts']:
                upload_part(file_info, fin, part)

    def upload_part(file_info, stream, part):
        udata = file_info.copy()
        udata.update(part)
        url = '{upload_url}/{partNo}'.format(**udata)
        stream.seek(part['startOffset'])
        data = stream.read(part['endOffset'] - part['startOffset'] + 1)
        raw_issue_request('PUT', url, data=data, binary=True)
        

    def change_article_status(article_id, status='public'):
        data = {'status': status}
        response = issue_request('PUT', f'account/articles/{article_id}', data=data)
        
    def publish_article(token, article_id):
        headers = {
            'Authorization': 'token ' + token,
            'Content-Type': 'application/json'
        }
    

        # This URL and method are hypothetical; refer to the API documentation for the actual values.
        url = BASE_URL.format(endpoint=f'account/articles/{article_id}/publish')

        response = requests.post(url, headers=headers)

        # Handle the response
        if response.status_code == 200:
            print("Article published successfully!")
        else:
            print("Error:", response.status_code)
            print(response.text)




    article_id = create_article(title)
    print("b")
    file_info = initiate_new_upload(article_id, filepath)
    print("c")
    upload_parts(file_info)
    print("d")
    complete_upload(article_id, file_info['id'])
    print("e")
    change_article_status(article_id, 'public')
    print("f")
    print("Make the file public")
    publish_article(token,article_id)
    print("g")
        


    
    

def main():
    """
    Main function to orchestrate the data processing.
    
    Example usage:
     - Without manifest files already downloaded:
            python getHCMIData.py -m path_to_manifest_file -t transcriptomics -o transcriptomics.csv
     - With manifest files already downloaded:
            python getHCMIData.py -M path_to_manifest_folder -t copy_number -o copy_number.csv
    """
    
    parser = argparse.ArgumentParser(description='Process data.')
    parser.add_argument('-m', '--manifest', help='Path to manifest file', required=True)
    parser.add_argument('-M', '--manifestfolder', help='Path to manifest folder', required=False)
    tc = ['transcriptomics', 'copy_number', 'mutations']
    parser.add_argument('-t', '--type', help='Type of data (e.g., transcriptomics, copy_number)',choices = tc, required=True)
# #     parser.add_argument('-f', '--figshare', help='Figshare data URL', required=True)
    parser.add_argument('-o', '--outname', help='Output CSV Name', required=True)
    args = parser.parse_args()
    
#     if bool(args.manifest) == bool(args.manifestfolder):  # Exactly one is required. Use Folder when data is already downloaded.
#         parser.error("Exactly one of --manifest or --manifestfolder is required.")

#     if bool(args.manifest):
#         download_option = True
#         manifest = args.manifest
#         # Ensure gdc-client is available and download if not.
#         ensure_gdc_client()
#     elif bool(args.manifestfolder):
#         download_option = False
#         manifest = args.manifestfolder
#         print("Using provided manifest folder without downloading data...")
        
        
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

    # Extract data files
    data_files = get_clean_files(args.type)

    # Retrieve figshare gene data for entrez map
    gene_url = "https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47"
    entrez_map_file = retrieve_figshare_data(gene_url)
    

    # Combine the data
    combined_data = map_and_combine(data_files, args.type, metadata, entrez_map_file)

    # Save to CSV
    write_dataframe_to_csv(combined_data, args.outname)
    
    
    print("Data processing complete!")
    token = "f356c7e4b66b50867b279bd16b78d2d81e0c420e9bda39074573e7b9426dae01535f51ba03190177a90b22dda514423135aaa13bdd347917a023cd2fe25a832c"
#     url = upload_to_figshare(args.outname, args.outname, token)
#     print(f"Article published at: {url}")

    print("attempting to push")
#     push_to_figshare(args.outname, args.outname, token)
    
   

    upload_to_figshare(token, args.outname, args.outname)



if __name__ == "__main__":
    main()
    