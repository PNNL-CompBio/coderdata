import os
import requests
import argparse
import hashlib
import json

def create_project(token, title):
    BASE_URL = 'https://api.figshare.com/v2/{endpoint}'
    headers = {'Authorization': 'token ' + token, 'Content-Type': 'application/json'}

    data = {
        "title": title,
        "description": "Cancer Organoid Data",
        "keywords": ["cancer", "organoid", "hcmi"],
        "categories_by_source_id": [
            "321101",
            "400207"
        ]
    }

    response = requests.post(BASE_URL.format(endpoint='account/projects'), headers=headers, data=json.dumps(data))
    response.raise_for_status()
    return response.json()['entity_id']
        
def upload_to_figshare(token, title, filepath, project_id,publish):
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

    def create_article(title,project_id):
        """
        Creates a new article on Figshare with a given title.
        """
        data = {
            'title': title,
            'description': "Cancer Organoid Data",
            'keywords': ["cancer","organoid","hcmi"],
            'defined_type': 'dataset',
            "categories_by_source_id": [
            "321101",
            "400207"
          ]
        }
        result = issue_request('POST', f'account/projects/{project_id}/articles', data=data)
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
    # 'account/projects/187725/articles
        url = BASE_URL.format(endpoint=f'account/articles/{article_id}/publish')
        response = requests.post(url, headers=headers)

        # Handle the response
        if response.status_code == 201:
            print("Article published successfully!")
        else:
            print("Error:", response.status_code)
            print(response.text)
            
    article_id = create_article(title,project_id)
    file_info = initiate_new_upload(article_id, filepath)
    upload_parts(file_info)
    complete_upload(article_id, file_info['id'])
    change_article_status(article_id, 'public')
    if publish == True:
        print("Make the file public")
        publish_article(token,article_id)
        
    
    
def main():
    parser = argparse.ArgumentParser(description='Upload files to Figshare.')
    parser.add_argument('-z', '--token', help='Figshare token ID', required=True)
    parser.add_argument('-t', '--title', help='Title for the Figshare article', required=True)
    parser.add_argument('-d', '--directory', help='Directory containing files to upload', required=True)
    parser.add_argument('-p', '--publish', help='Should this be published?', required=True)
    args = parser.parse_args()

    project_id = create_project(args.token,args.title)
            
    for file in os.listdir(args.directory):
        file = os.path.join(args.directory, file)
        print(f"Running 'upload_to_figshare' for {file}")
        upload_to_figshare(args.token, file, file, project_id,args.publish)

if __name__ == "__main__":
    main()

