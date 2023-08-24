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
        "mutation": "ensemble_masked.maf.gz"
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
        a = 2**float(a)
        if a < 0.5210507:
            return 'deep del'
        elif a < 0.7311832:
            return 'het loss'
        elif a < 1.214125:
            return 'normal'
        elif a < 1.731183:
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

        elif data_type == "mutation":
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
#     parser.add_argument('-f', '--figshare', help='Figshare data URL', required=True)
    parser.add_argument('-o', '--outname', help='Output CSV Name', required=True)
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

    # Extract data files
    print("Extracting data files...")
    data_files = get_clean_files(args.type)

    # Retrieve figshare gene data for entrez map
    print("Retrieving figshare gene data for entrez map...")
    gene_url = "https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47"
    entrez_map_file = retrieve_figshare_data(gene_url)
    

    # Combine the data
    print("Joining the dataframes...")
    combined_data = map_and_combine(data_files, args.type, metadata, entrez_map_file)

    # Save to CSV
    print("Saving data to CSV...")
    write_dataframe_to_csv(combined_data, args.outname)
    
    print("Data processing complete!")


if __name__ == "__main__":
    main()
    
