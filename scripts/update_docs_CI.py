import os
import pandas as pd
import yaml

def read_yaml_file(file_path):
    """
    Read and parse the YAML file.
    """
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def list_files(directory, yaml_data):
    """
    List all files in the given directory and categorize them.
    Include links from the YAML data.
    """
    files = os.listdir(directory)
    categorized_files = {}
    for file in files:
        category = file.split('_')[0]
        file_info = {
            'name': file,
            'download_link': yaml_data['file_download'].get(file, ''),
            'file_url': yaml_data['file_url'].get(file, ''),
            'type': os.path.splitext(file)[1]
        }
        categorized_files.setdefault(category, []).append(file_info)
    return categorized_files

def create_docs_section(category, files, directory):
    """
    Create a Markdown section for a given category of files.
    """
    section = f"## {category.capitalize()} Data\n\n"
    for file_info in files:
        file_path = os.path.join(directory, file_info['name'])
        if file_info['type'] in ['.png']:
            # Embed PNG images
            section += f"![{file_info['name']}]({file_path})\n"
        elif file_info['type'] in ['.csv']:
            # Display CSV content as a table
            try:
                df = pd.read_csv(file_path)
                section += df.to_markdown(index=False) + "\n"
            except Exception as e:
                section += f"Error reading {file_info['name']}: {e}\n"
        else:
            # For other file types, provide a link
            section += f"- [{file_info['name']}](./docs/{file_info['name']})\n"
            if file_info['download_link']:
                section += f"  - [Download]({file_info['download_link']})\n"
            if file_info['file_url']:
                section += f"  - [View on Figshare]({file_info['file_url']})\n"
    return section

def generate_documentation(directory, yaml_file):
    """
    Generate the documentation from files in the specified directory.
    """
    yaml_data = read_yaml_file(yaml_file)
    categorized_files = list_files(directory, yaml_data)
    documentation = "# Cancer Omics and Drug Experiment Response Data (`coderdata`) Python Package\n\n"

    # Add predefined sections
    documentation += "### Installation\n\n...\n\n### Use\n\nHow to use the package\n\n"

    # Add sections for each category
    for category, files in categorized_files.items():
        documentation += create_docs_section(category, files, directory)

    # Add combined summary section
    documentation += "## Combined Data Summary\n\nTODO: create script that generates figures\n\n"

    # Save to a Markdown file
    with open("docs/index.md", "w") as doc_file:
        doc_file.write(documentation)

    print("Documentation generated successfully.")

# Replace 'docs' with the path to your documentation directory and 'docs/figshare_latest.yml' with the path to your YAML file
generate_documentation('docs', 'docs/figshare_latest.yml')
