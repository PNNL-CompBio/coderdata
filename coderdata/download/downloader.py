# coderdata/download/downloader.py

import requests
import os
import json
import yaml


def download_data_by_prefix(dataset_prefix=None):
    """
    Download the most recent version of files from a Figshare dataset, filtered by a specific prefix or all files.

    This function queries the Figshare API to retrieve details of a dataset and then downloads files from it.
    Files can be filtered by a specified prefix such as hcmi, beataml, etc. If 'all', an empty string, or None is passed as the prefix,
    all files in the dataset are downloaded. The function identifies the most recent version of a file 
    by selecting the one with the highest ID among duplicates with the same name.

    Parameters
    ----------
    dataset_prefix : str, optional
        The prefix of the dataset to download (e.g., 'hcmi'). If 'all', an empty string, or None, 
        all files in the dataset are downloaded. Default is None.

    Returns
    -------
    None
        The function downloads files to the local repository and does not return any value.
    """
    
    # Get the dataset details
    url = "https://api.figshare.com/v2/articles/25030490"
    
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Failed to get dataset details from Figshare: {response.text}")

    data = response.json()

    # Filter files by the specified prefix
    if dataset_prefix and dataset_prefix.lower() != "all":
        filtered_files = [file for file in data['files'] if file['name'].startswith(dataset_prefix)]
    else:
        filtered_files = data['files']

    # Group files by name and select the one with the highest ID
    unique_files = {}
    for file in filtered_files:
        file_name = file['name']
        file_id = file['id']
        if file_name not in unique_files or file_id > unique_files[file_name]['id']:
            unique_files[file_name] = {'file_info': file, 'id': file_id}

    for file_name, file_data in unique_files.items():
        file_info = file_data['file_info']
        file_url = file_info['download_url']

        # Download the file
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            with open(file_name, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192): 
                    f.write(chunk)

        print(f"Downloaded {file_name} to local repository.")

    return

