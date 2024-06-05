import argparse
import yaml
import sys
import subprocess
import os

def update_version(file_path, new_version):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    with open(file_path, 'w') as file:
        for line in lines:
            if "version=" in line:
                line = f"    version='{new_version}',"
            file.write(line)

def update_article_url(yaml_file_path, downloader_file_path):
    # Read the article link from the YAML file
    with open(yaml_file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)
    
    article_link = yaml_data.get('article_link', '')
    # Extract the article ID from the URL
    article_id = article_link.rsplit('/', 1)[-1]

    # Construct the new API URL
    new_url = f"https://api.figshare.com/v2/articles/{article_id}"

    # Update the URL in the downloader.py file
    with open(downloader_file_path, 'r') as file:
        lines = file.readlines()

    with open(downloader_file_path, 'w') as file:
        for line in lines:
            if 'url =' in line and 'figshare.com' in line:
                file.write(f'    url = "{new_url}"\n')
            else:
                file.write(line)

def package_and_upload():
    # Command to create distribution packages
    package_cmd = ['python3', 'setup.py', 'sdist', 'bdist_wheel']
    # Command to upload packages using Twine
    upload_cmd = ['twine', 'upload', 'dist/*', '--verbose', '-u', '__token__', '-p', os.getenv('PYPI_TOKEN')]

    # Run package command
    result = subprocess.run(package_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("Package creation failed:", result.stderr)
        return

    # Run upload command
    result = subprocess.run(upload_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("Upload failed:", result.stderr)
    else:
        print("Upload successful:", result.stdout)

def main():
    parser = argparse.ArgumentParser(description="Update article URL in downloader script.")
    parser.add_argument('-y',"--yaml", help="Path to the figshare_latest.yml file")
    parser.add_argument('-d',"--downloader", help="Path to the downloader.py file")
    parser.add_argument('-v',"--version", help="Version to update setup.py with")

    args = parser.parse_args()

    update_version('setup.py', args.version)
    update_article_url(args.yaml, args.downloader)
    package_and_upload()

if __name__ == "__main__":
    main()