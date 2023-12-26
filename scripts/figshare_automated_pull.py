import requests
import wget
import os
import shutil
import gzip
import pandas as pd
import re
import argparse

def retrieve_figshare_data(datatype, modeltype):
    
    """
    *Need to be in a directory void of data files*
    
    Downloads data from FigShare urls into your directory.
    
    Args:
        datatype: list of data types to download. (options: "samples", "genes", "experiments", "drugs",
                                                "rppa". "proteomics", "miRNA", "mutation",
                                                "transcriptomics", "copy_number", "methylation",
                                                "drugs_by_structure", "newid_experiments",
                                                "reduce_experiments", "CNV", "proteomics")
        modeltype: list of model types to download. (options: "cellline", "patient", "HCMI", "beatAML")
        
        
    Example Usage: 
    retrieve_figshare_data(['samples'], ['cellline', 'patient'])
    
    Returns: 
    a list of downloaded files.
    """
    
    
    figshare_urls = {'cellline_samples': 'https://figshare.com/ndownloader/files/40576103?private_link=525f7777039f4610ef47',
                 'cellline_genes': 'https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47',
                 'cellline_experiments': 'https://figshare.com/ndownloader/files/41259270?private_link=525f7777039f4610ef47',
                 'cellline_drugs': 'https://figshare.com/ndownloader/files/41259273?private_link=525f7777039f4610ef47',
                 'cellline_rppa': 'https://figshare.com/ndownloader/files/41466699?private_link=525f7777039f4610ef47',
                 'cellline_proteomics': 'https://figshare.com/ndownloader/files/41466702?private_link=525f7777039f4610ef47',
                 'cellline_miRNA': 'https://figshare.com/ndownloader/files/42120534?private_link=525f7777039f4610ef47',
                 'cellline_mutations': 'https://figshare.com/ndownloader/files/42131268?private_link=525f7777039f4610ef47',
                 'cellline_transcriptomics': 'https://figshare.com/ndownloader/files/42131304?private_link=525f7777039f4610ef47',
                 'cellline_copy_number': 'https://figshare.com/ndownloader/files/42131325?private_link=525f7777039f4610ef47',
                 'cellline_methylation': 'https://figshare.com/ndownloader/files/42131337?private_link=525f7777039f4610ef47',
                 'cellline_drugs_by_structure': 'https://figshare.com/ndownloader/files/42357210?private_link=525f7777039f4610ef47',
                 'cellline_newid_experiments': 'https://figshare.com/ndownloader/files/42357213?private_link=525f7777039f4610ef47',
                 'cellline_reduce_experiments': 'https://figshare.com/ndownloader/files/42357216?private_link=525f7777039f4610ef47',
                 'patient_data_samples': 'https://figshare.com/ndownloader/files/42147513?private_link=7ffe48478ec907b36dfb',
                 'patient_data_somatic_mutation': 'https://figshare.com/ndownloader/files/42147516?private_link=7ffe48478ec907b36dfb',
                 'patient_data_CNV': 'https://figshare.com/ndownloader/files/42147519?private_link=7ffe48478ec907b36dfb',
                 'patient_data_transcriptomics': 'https://figshare.com/ndownloader/files/42147522?private_link=7ffe48478ec907b36dfb',
                 'patient_data_proteomics': 'https://figshare.com/ndownloader/files/42147525?private_link=7ffe48478ec907b36dfb',
                 'HCMI_copy_number': 'https://figshare.com/ndownloader/files/42211392?private_link=46a4aefc42e47fe1fb6d',
                 'HCMI_mutations': 'https://figshare.com/ndownloader/files/42211395?private_link=46a4aefc42e47fe1fb6d',
                 'HCMI_transcriptomics': 'https://figshare.com/ndownloader/files/42211398?private_link=46a4aefc42e47fe1fb6d',
                 'beatAML_samples': 'https://figshare.com/ndownloader/files/42289053',
                 'beatAML_drugs': 'https://figshare.com/ndownloader/files/42357918',
                 'beatAML_experiments': 'https://figshare.com/ndownloader/files/42357921',
                 'beatAML_mutations': 'https://figshare.com/ndownloader/files/42357924v',
                 'beatAML_proteomics': 'https://figshare.com/ndownloader/files/42357927',
                 'beatAML_transciptomics': 'https://figshare.com/ndownloader/files/42357930'
                }
    
    #collecting url's for datasets
    url_list = []
    
    for key in figshare_urls:
        for x in datatype:
            if x in key:
                for y in modeltype:
                    if y in key:
                        url_list.append(figshare_urls[key])

    #downloading the datasets and adding them to a list
    files = []
    
    for url in url_list:
        files_0 = os.listdir()
        wget.download(url)
        files_1 = os.listdir()
        figdir = str(next(iter((set(files_1) - set(files_0)))))
        files.append(figdir)
    return files

def merger(*data_types, directory=".", outname="Merged_Data.csv", no_duplicate=True, drop_na=False):
    """
    combines datasets of chosen data types.
    
    Args:
        data_types: Type of data sets to merge. (options: "samples", "genes", "experiments", "drugs",
                                                "rppa". "proteomics", "miRNA", "mutation",
                                                "transcriptomics", "copy_number", "methylation",
                                                "drugs_by_structure", "newid_experiments",
                                                "reduce_experiments", "CNV", "proteomics")
        directory: Directory where the data resides.
        outname: Name of the output file.
        no_duplicate: Drop duplicate rows if set to True. Default is True.
        drop_na: Drop rows with NA values if set to True. Default is False.
        
    Example Usage:
    python figshare_pull.py samples mutations copy_number methylation proteomics -d . -o Merged_Data.csv --datatype ["samples", "mutations", "copy_number", "methylation"] --modeltype ["cellline", "HCMI"]  --no_duplicate True --drop_na False

    Returns:
        DataFrame: Merged dataset of desired data types.
    """

    dfs = {}    
    merged_datasets = []

    for data_type in data_types:
        files = [f for f in os.listdir(directory) if data_type in f and f.endswith(('.csv', '.tsv', '.csv.gz', '.tsv.gz'))]
        
        if not files:
            print(f"No files found for data type: {data_type}. This data type will not be included.")
            continue

        selected_files = files
        print(f"Selected file(s) for {data_type}: {selected_files}. Proceeding with merge.")

        
        datatype_merged_list = []
        for selected_file in selected_files:
            path = os.path.join(directory, selected_file)
            compression = 'gzip' if selected_file.endswith('.gz') else None
            delimiter = "\t" if selected_file.endswith((".tsv", ".tsv.gz")) else ","
            chunk_iter = pd.read_csv(path, sep=delimiter, compression=compression, chunksize=10**5, low_memory=False)
            df_parts = [chunk for chunk in chunk_iter]
            dfs[data_type] = pd.concat(df_parts, ignore_index=True)
            
            single_df = dfs[data_type]
            datatype_merged_list.append(single_df) 
    
        datatype_datasets_merged = pd.concat(datatype_merged_list, ignore_index = True)
        merged_datasets.append(datatype_datasets_merged)
    
    merged_df = pd.concat(merged_datasets, axis = 0, ignore_index = True)
    if no_duplicate:
        merged_df.drop_duplicates(inplace=True)
    if drop_na:
        merged_df.dropna(inplace=True)

    merged_df.to_csv(outname, index=False)
    return merged_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge datasets by data types from specified directory.")
    parser.add_argument("data_types", nargs="+", help="type of datasets to merge. (Example: 'transcriptomics', 'mutations')")
    parser.add_argument("datatype", nargs="+", help="list of data types to download. (Example: ['transcriptomics', 'mutations'])")
    parser.add_argument("modeltype", nargs="+", help="list of model types to download. (Example: '[cellline', 'HCMI'])")
    parser.add_argument("-d", "--directory", default=".", help="Directory where the data resides.")
    parser.add_argument("-o", "--outname", help="Name of the output file.")
    parser.add_argument("-u", "--url", help="URL to figshare.")
    parser.add_argument("--no_duplicate", action="store_true", help="Drop duplicate rows.")
    parser.add_argument("--drop_na", action="store_true", help="Drop rows with NA values.")

    args = parser.parse_args()
#     type_list = ['samples', 'mutation', 'copy_number']
#     model_list = ['celline', 'HCMI']
#     cell_line_files = retrieve_figshare_data(type_list, model_list)
    retrieve_figshare_data(datatype = args.datatype, modeltype = args.modeltype)
    print("\n")
    merger(*args.data_types, directory=args.directory, outname=args.outname, no_duplicate=args.no_duplicate, drop_na=args.drop_na)

