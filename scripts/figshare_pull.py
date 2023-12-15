import requests
import wget
import os
import shutil
import gzip
import pandas as pd
import re
import argparse

def retrieve_figshare_data(url):
    """
    Download data from FigShare url.
    
    Returns: a list of newly unpacked files.
    """
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    figdir = str(next(iter((set(files_1) - set(files_0)))))
    shutil.unpack_archive(figdir)
    files_2 = os.listdir()
    files_3 = list(set(files_2) - set(files_1))
    files_4 = os.listdir()
    files_5 = list(set(files_4) - set(files_1))
    return files_5

def merger(*data_types, directory=".", outname="Merged_Data.csv", no_duplicate=True, drop_na=False):
    """
    Merge datasets by data types.
    
    Args:
        data_types: Type of data sets to merge. (Example: "transcriptomics", "mutations")
        directory: Directory where the data resides.
        outname: Name of the output file.
        no_duplicate: Drop duplicate rows if set to True. Default is True.
        drop_na: Drop rows with NA values if set to True. Default is False.
        
    Example Usage:
    python figshare_pull.py transcriptomics mutations copy_number methylation proteomics -d . -o Merged_Data.csv --url "https://figshare.com/ndownloader/articles/22822286?private_link=525f7777039f4610ef47"  --no_duplicate --drop_na

    Returns:
        DataFrame: Merged dataset.
    """

    def get_prefix_number(filename):
        """Extract prefix number or return 0 if none exists."""
        match = re.match(r"(\d+)", filename)
        return int(match.group(1)) if match else 0

    dfs = {}
    primary_data_type = None
    
    for data_type in data_types:
        files = [f for f in os.listdir(directory) if data_type in f and f.endswith(('.csv', '.tsv', '.csv.gz', '.tsv.gz'))]
        
        if not files:
            print(f"No files found for data type: {data_type}. This data type will not be included.")
            continue

        selected_file = max(files, key=get_prefix_number)
        print(f"Selected file for {data_type}: {selected_file}. Proceeding with merge.")
        
        if not primary_data_type:
            primary_data_type = data_type

        path = os.path.join(directory, selected_file)
        compression = 'gzip' if selected_file.endswith('.gz') else None
        delimiter = "\t" if selected_file.endswith((".tsv", ".tsv.gz")) else ","
        chunk_iter = pd.read_csv(path, sep=delimiter, compression=compression, chunksize=10**5, low_memory=False)
        df_parts = [chunk for chunk in chunk_iter]
        dfs[data_type] = pd.concat(df_parts, ignore_index=True)

    if not primary_data_type:
        print("No suitable data found for any specified data type.")
        return

    merged_df = dfs[primary_data_type]

    for data_type in data_types[1:]:
        if data_type in dfs:
            merge_cols = ["improve_sample_id", "entrez_id"]
            if "source" in merged_df.columns and "source" in dfs[data_type].columns:
                merge_cols.append("source")
            if "study" in merged_df.columns and "study" in dfs[data_type].columns:
                merge_cols.append("study")
            merged_df = merged_df.merge(dfs[data_type], on=merge_cols, how="outer")
    if no_duplicate:
        merged_df.drop_duplicates(inplace=True)
    if drop_na:
        merged_df.dropna(inplace=True)

    merged_df.to_csv(outname, index=False)
    return merged_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge datasets by data types from specified directory.")
    parser.add_argument("data_types", nargs="+", help="Types of data sets to merge. (Example: 'transcriptomics', 'mutations')")
    parser.add_argument("-d", "--directory", default=".", help="Directory where the data resides.")
    parser.add_argument("-o", "--outname", help="Name of the output file.")
    parser.add_argument("-u", "--url", help="URL to figshare.")
    parser.add_argument("--no_duplicate", action="store_true", help="Drop duplicate rows.")
    parser.add_argument("--drop_na", action="store_true", help="Drop rows with NA values.")

    args = parser.parse_args()
#     cell_line_data = "https://figshare.com/ndownloader/articles/22822286?private_link=525f7777039f4610ef47"
#     cell_line_files = retrieve_figshare_data(cell_line_data)
    retrieve_figshare_data(args.url)
    print("\n")
    merger(*args.data_types, directory=args.directory, outname=args.outname, no_duplicate=args.no_duplicate, drop_na=args.drop_na)

