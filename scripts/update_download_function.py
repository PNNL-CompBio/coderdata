import argparse
import yaml

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

def main():
    parser = argparse.ArgumentParser(description="Update article URL in downloader script.")
    parser.add_argument('-y',"--yaml", help="Path to the figshare_latest.yml file")
    parser.add_argument('-d',"--downloader", help="Path to the downloader.py file")
    args = parser.parse_args()

    update_article_url(args.yaml, args.downloader)

if __name__ == "__main__":
    main()
