import os
import requests
import argparse
import hashlib
import json
import time
import yaml


def upload_to_figshare(token, title, directory, project_id, publish, article_id=None):
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
            return json.loads(response.content)
        except ValueError:
            return response.content

    def issue_request(method, endpoint, *args, **kwargs):
        """
        Sends a request to a specific Figshare API endpoint.
        """
        return raw_issue_request(method, BASE_URL.format(endpoint=endpoint), *args, **kwargs)

    def create_article(title, project_id):
        """
        Creates a new article on Figshare with a given title.
        """
        data = {
            'title': title,
            'description': "Cancer Organoid Data",
            'keywords': ["cancer", "organoid", "hcmi"],
            'defined_type': 'dataset',
            "categories_by_source_id": ["321101", "400207"]
        }
        result = issue_request('POST', f'account/projects/{project_id}/articles', data=data)
        result = raw_issue_request('GET', result['location'])
        return result['id']

    def get_file_check_data(file_path):
        """
        Check the MD5 checksum and size of a given file.
        """
        with open(file_path, 'rb') as fin:
            md5 = hashlib.md5()
            size = 0
            data = fin.read(CHUNK_SIZE)
            while data:
                size += len(data)
                md5.update(data)
                data = fin.read(CHUNK_SIZE)
            return md5.hexdigest(), size

    def initiate_new_upload(article_id, file_path):
        """
        Initiates the file upload process for a specific article.
        """
        endpoint = f'account/articles/{article_id}/files'
        md5, size = get_file_check_data(file_path)
        data = {'name': os.path.basename(file_path), 'md5': md5, 'size': size}
        result = issue_request('POST', endpoint, data=data)
        return raw_issue_request('GET', result['location'])

    def upload_parts(file_info, file_path):
        """
        Handles the multipart upload for large files.
        """
        url = file_info['upload_url']
        result = raw_issue_request('GET', url)
        with open(file_path, 'rb') as fin:
            for part in result['parts']:
                upload_part(file_info, fin, part)

    def upload_part(file_info, stream, part):
        """
        Uploads a specific part/chunk of the file.
        """
        udata = file_info.copy()
        udata.update(part)
        url = f"{udata['upload_url']}/{udata['partNo']}"
        stream.seek(udata['startOffset'])
        data = stream.read(udata['endOffset'] - udata['startOffset'] + 1)
        raw_issue_request('PUT', url, data=data, binary=True)

    def complete_upload(article_id, file_id):
        """
        Upload article to Figshare
        """
        issue_request('POST', f'account/articles/{article_id}/files/{file_id}')

    def change_article_status(article_id, status='public'):
        """
        Changes the visibility status of a specific article on Figshare.
        """
        issue_request('PUT', f'account/articles/{article_id}', data={'status': status})

    def publish_article(article_id):
        """
        Publish article on Figshare
        """
        issue_request('POST', f'account/articles/{article_id}/publish')

    def get_remote_file_info(article_id, file_name):
        """
        Extract info from files on Figshare
        """
        existing_files = issue_request('GET', f'account/articles/{article_id}/files')
        for file in existing_files:
            if file['name'] == file_name:
                return file
        return None

    def upload_file_with_retry(file_info, file_path, max_retries=15):
        """
        Upload files and retry 15 times if fails.  
        """
        for attempt in range(max_retries):
            try:
                upload_parts(file_info, file_path)
                return
            except requests.exceptions.ConnectionError as e:
                print(f"Upload failed, retrying {attempt + 1}/{max_retries}. Error: {e}")
                time.sleep(2 ** attempt) 
        raise Exception("Max retries reached. Upload failed.")
    

    def create_or_get_article(title, project_id, article_id):
        """
        Create new article or just use old article id. User selected
        """
        if article_id:
            return article_id
        else:
            return create_article(title, project_id)

    def delete_existing_file(article_id, file_id):
        """
        Delete file on figshare. Done when replacing a file.
        """
        issue_request('DELETE', f'account/articles/{article_id}/files/{file_id}')


    def write_figshare_details_to_yaml(article_id, project_id, title):
        """
        Write details of Figshare to yaml
        """
        #convert slashes and periods to underscores so the file links are generated correctly.
        title_updated = title.replace('/', '_')
        title_updated = title_updated.replace('.', '_')
        article_info = issue_request('GET', f'articles/{article_id}')
        # article_link = f"https://figshare.com/articles/dataset/{title}/{project_id}/file/{article_id}"
        article_link = f"https://figshare.com/articles/dataset/{title_updated}/{article_id}"

        # Retrieve the article details
        article_details_response = requests.get(article_info['url'])
        article_details_response.raise_for_status()
        article_details = article_details_response.json()

        # Construct the URLs
        file_url_links = {file['name']:f"https://figshare.com/articles/dataset/{title_updated}/{article_id}?file={file['id']}" for file in article_details['files']}
        file_download_link = {file['name']: file['download_url'] for file in article_details['files']}
        yaml_data = {
            'article_link': article_link,
            'file_url': file_url_links,
            'file_download': file_download_link
        }

        with open('/tmp/figshare_latest.yml', 'w') as file:
            yaml.dump(yaml_data, file, default_flow_style=False)


    article_id = create_or_get_article(title, project_id, article_id)
    all_files_uploaded = True

    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)

        remote_file_info = get_remote_file_info(article_id, file_name)
        if remote_file_info:
            local_md5, local_size = get_file_check_data(file_path)
            if remote_file_info['size'] != local_size or remote_file_info['computed_md5'] != local_md5:
                print(f"Updating file {file_name} in Figshare...")
                delete_existing_file(article_id, remote_file_info['id'])
                file_info = initiate_new_upload(article_id, file_path)
                upload_file_with_retry(file_info, file_path)
                complete_upload(article_id, file_info['id'])
            else:
                print(f"File {file_name} already exists and is up to date, skipping...")
                continue
        else:
            file_info = initiate_new_upload(article_id, file_path)
            upload_file_with_retry(file_info, file_path)
            complete_upload(article_id, file_info['id'])

    if all_files_uploaded and publish:
        change_article_status(article_id, 'public')
        publish_article(article_id)
        print("Article published successfully!")
    elif not all_files_uploaded:
        print("Not all files were uploaded successfully. Article was not published.")
    elif all_files_uploaded and not publish:
        print("Files uploaded successfully but not published.")
        
    if all_files_uploaded:
        write_figshare_details_to_yaml(article_id, project_id,title)

def main():
    parser = argparse.ArgumentParser(description='Upload files to Figshare.')
    parser.add_argument('-z', '--token', help='Figshare token ID', required=True)
    parser.add_argument('-t', '--title', help='Title for the Figshare article', required=True)
    parser.add_argument('-d', '--directory', help='Directory containing files to upload', required=True)
    parser.add_argument('-p', '--publish', help='Publish the article', action='store_true')
    parser.add_argument('-j', '--project_id', help='Existing Figshare project ID', required=True)
    parser.add_argument('-a', '--article_id', help='Existing Figshare article ID', required=False, default=None)
    args = parser.parse_args()

    upload_to_figshare(args.token, args.title, args.directory, args.project_id, args.publish, args.article_id)

if __name__ == "__main__":
    main()