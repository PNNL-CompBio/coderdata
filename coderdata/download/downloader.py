# coderdata/download/downloader.py

from importlib import resources
from hashlib import md5
from pathlib import Path
from os import PathLike
import os
import requests
import warnings
import yaml
from typing import Iterable, List, Dict, Any, Optional
    
    
    
def _gather_files_from_response(resp: requests.Response) -> List[Dict[str, Any]]:
    """
    Normalize Figshare API responses into a list of file dicts.

    Supports:
      1) Article endpoint:   https://api.figshare.com/v2/articles/{id}
         -> JSON object with key 'files' (list)

      2) Files endpoint:     https://api.figshare.com/v2/articles/{id}/files[?...]
         -> JSON list of file objects (possibly paginated with Link headers)
    """
    data = resp.json()
    if isinstance(data, dict) and "files" in data and isinstance(data["files"], list):
        return data["files"]
    if isinstance(data, list):
        return data
    raise ValueError("Unexpected Figshare API response structure; expected dict with 'files' "
                     "or a list of file objects.")


def _iter_paginated_files(url: str, session: Optional[requests.Session] = None) -> Iterable[Dict[str, Any]]:
    """
    Iterate over all files, following 'Link: <...>; rel=\"next\"' pagination if present.
    Works for both the article endpoint (no pagination) and the files endpoint (may paginate).
    """
    sess = session or requests.Session()
    next_url = url

    while next_url:
        resp = sess.get(next_url)
        if resp.status_code != 200:
            raise Exception(f"Failed to get dataset details from Figshare: {resp.text}")

        for f in _gather_files_from_response(resp):
            yield f

        # RFC5988-style 'Link' header pagination
        link = resp.headers.get("Link") or resp.headers.get("link")
        next_url = None
        if link:
            parts = [p.strip() for p in link.split(",")]
            for part in parts:
                if 'rel="next"' in part:
                    start = part.find("<") + 1
                    end = part.find(">", start)
                    if start > 0 and end > start:
                        next_url = part[start:end]
                        break

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
        local_path.mkdir(parents=True, exist_ok=True)
    # Get the dataset details
    with resources.open_text('coderdata', 'dataset.yml') as f:
        data_information = yaml.load(f, Loader=yaml.FullLoader)
    url = data_information['figshare']

    name = (name or "all").casefold()
    session = requests.Session()
    all_files = list(_iter_paginated_files(url, session=session))

    if name != "all":
        filtered_files = [
            f for f in all_files
            if (f.get('name', '').casefold().startswith(name)) or ('genes' in f.get('name', '').casefold())
        ]
    else:
        filtered_files = all_files

    unique_files = {}
    for file in filtered_files:
        fname = file.get('name')
        fid = file.get('id')
        if fname is None or fid is None:
            continue
        file_name = local_path.joinpath(fname)
        if (file_name not in unique_files) or (fid > unique_files[file_name]['id']):
            unique_files[file_name] = {'file_info': file, 'id': fid}

    for file_name, file_data in unique_files.items():
        file_info = file_data['file_info']
        file_id = str(file_info['id'])
        file_url = f"https://api.figshare.com/v2/file/download/{file_id}"
        file_md5sum = file_info.get('supplied_md5')

        if file_name.exists() and not exist_ok:
            warnings.warn(
                f"{file_name} already exists. Use argument 'exist_ok=True' to overwrite the existing file."
            )

        retry_count = 10
        while retry_count > 0:
            with session.get(file_url, stream=True) as r:
                r.raise_for_status()
                with open(file_name, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            if file_md5sum:
                with open(file_name, 'rb') as f:
                    check_md5sum = md5(f.read()).hexdigest()
                if file_md5sum == check_md5sum:
                    break
                else:
                    retry_count -= 1
                    if retry_count > 0:
                        warnings.warn(
                            f"{file_name} failed MD5 verification "
                            f"(expected: {file_md5sum}, got: {check_md5sum}). Retrying..."
                        )
            else:
                break

        if retry_count == 0 and file_md5sum:
            warnings.warn(
                f"{file_name} could not be downloaded with a matching MD5 after retries."
            )
        else:
            print(f"Downloaded '{file_url}' to '{file_name}'")


