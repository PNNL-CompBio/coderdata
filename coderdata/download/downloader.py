# coderdata/download/downloader.py

from importlib import resources
from pathlib import Path
from os import PathLike
import os
import requests
import warnings

import yaml

def download(
        name: str='all',
        local_path: PathLike=Path.cwd(),
        exist_ok: bool=False
        ):
    """
    Download the most recent version of files from a Figshare dataset,
    filtered by a specific prefix or all files.

    This function queries the Figshare API to retrieve details of a
    dataset and then downloads files from it. Files can be filtered by a
    specified prefix such as hcmi, beataml, etc. If 'all', an empty
    string, or None is passed as the prefix, all files in the dataset
    are downloaded. The function identifies the most recent version of a
    file by selecting the one with the highest ID among duplicates with
    the same name.

    Parameters
    ----------
    dataset_prefix : str, optional
        The prefix of the dataset to download (e.g., 'hcmi'). If 'all',
        an empty string, or None, all files in the dataset are
        downloaded. Default is None.

    Returns
    -------
    None
        The function downloads files to the local repository and does
        not return any value.
    """
    
    # Create Path object from `local_path`
    if type(local_path) != Path:
        local_path = Path(local_path)

    if not local_path.exists():
        Path.mkdir(local_path)
    # Get the dataset details
    with resources.open_text('coderdata', 'dataset.yml') as f:
        data_information = yaml.load(f, Loader=yaml.FullLoader)
    url = data_information['figshare']
    
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(
            f"Failed to get dataset details from Figshare: {response.text}"
            )

    data = response.json()

    # making sure that we are case insensitive
    name = name.casefold()

    # Filter files by the specified prefix
    if name != "all":
        filtered_files = [
            file 
            for file 
            in data['files'] 
            if file['name'].startswith(name) or 'genes' in file['name']
            ]
    else:
        filtered_files = data['files']

    # Group files by name and select the one with the highest ID
    unique_files = {}
    for file in filtered_files:
        file_name = local_path.joinpath(file['name'])
        file_id = file['id']
        if (
            file_name not in unique_files
            or file_id > unique_files[file_name]['id']
        ):
            unique_files[file_name] = {'file_info': file, 'id': file_id}

    for file_name, file_data in unique_files.items():
        file_info = file_data['file_info']
        file_url = file_info['download_url']

        # Download the file
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            if file_name.exists() and not exist_ok:
                warnings.warn(
                    f"{file_name} already exists. Use argument 'exist_ok=True'"
                    "to overwrite existing file."
                    )
            else:
                with open(file_name, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192): 
                        f.write(chunk)

        print(f"Downloaded '{file_url}' to '{file_name}'")

    return

