import synapseclient
import pandas as pd
import os
import wget
from synapseclient import Project, Folder, File, Link
import requests
import numpy as np
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep, time
from datetime import datetime
import pubchem_retrieval


def retrieve_figshare_data(url):
    """
    Retrieve data from a given figshare URL and identify the downloaded file by its name.

    This function uses the wget tool to download a file from the provided figshare URL.
    By comparing the directory contents before and after the download, 
    it identifies the newly downloaded file's name.

    Parameters
    ----------
    url : str
        The figshare URL pointing to the data to be downloaded.

    Returns
    -------
    str
        Name of the downloaded file, e.g., "genes.csv".
    """
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    new_file = str(next(iter(set(files_1) - set(files_0))))
    return new_file

def download_from_github(raw_url, save_path):
    """
    Download a file from a raw GitHub URL and save to the specified path.

    Parameters
    ----------
    raw_url : str
        The raw GitHub URL pointing to the file to be downloaded.
    save_path : str
        The local path where the downloaded file will be saved.

    Returns
    -------
    None
    """
    response = requests.get(raw_url)
    with open(save_path, 'wb') as f:
        f.write(response.content)
    return



def generate_raw_drug_file(original_drug_file, sample_mapping_file, updated_raw_drug_file, wave_file):
    """
    Generate the formatted raw drug file to input into the curve fitting algorithm.
    Outputs {updated_raw_drug_file}.0 into the local directory

    Parameters
    ----------
    original_drug_file : str
        File path to the original unaltered drug file.
    sample_mapping_file : str
        File path to the sample mapping. Used to get labID which later maps to improve.
    updated_raw_drug_file : str
        File path to the newly generated raw drug file.
    wave_file: Supplimentary file with wave information
    Returns
    -------
    None
    """
    raw_drug = pd.read_csv(original_drug_file, sep="\t")
    drug_mod = raw_drug.rename(columns={"well_concentration": "DOSE",
                                        "normalized_viability": "GROWTH",
                                        "inhibitor": "DRUG"})
    
    sample_map = pd.read_excel(sample_mapping_file)
    primary_mapping = drug_mod.merge(sample_map[['dbgap_rnaseq_sample', 'labId']], 
                                     left_on='dbgap_rnaseq_sample', 
                                     right_on='dbgap_rnaseq_sample', 
                                     how='left')
    secondary_mapping = drug_mod.merge(sample_map[['dbgap_dnaseq_sample', 'labId']], 
                                       left_on='dbgap_dnaseq_sample', 
                                       right_on='dbgap_dnaseq_sample', 
                                       how='left')
    drug_mod['CELL'] = primary_mapping['labId'].combine_first(secondary_mapping['labId'])
    
    cohort_data = pd.read_excel(wave_file, skiprows=1)
    cohort_merge = drug_mod.merge(cohort_data[['dbgap_rnaseq_sample', 'cohort']],
                                  left_on='dbgap_rnaseq_sample',
                                  right_on='dbgap_rnaseq_sample',
                                  how='left')
    cohort_merge_secondary = drug_mod.merge(cohort_data[['dbgap_dnaseq_sample', 'cohort']],
                                            left_on='dbgap_dnaseq_sample',
                                            right_on='dbgap_dnaseq_sample',
                                            how='left')
    
    drug_mod['SOURCE'] = cohort_merge['cohort'].combine_first(cohort_merge_secondary['cohort'])
    drug_mod['STUDY'] = "BeatAML"
    drug_mod = drug_mod[['CELL', 'DRUG', 'DOSE', 'GROWTH', 'SOURCE', 'STUDY']]
    drug_mod = drug_mod.sort_values(by=['DRUG', 'CELL'])
    drug_mod = drug_mod.groupby(['DRUG', 'CELL']).filter(lambda x: len(x) >= 5)
    drug_mod.to_csv(updated_raw_drug_file, index=False, sep="\t")
    return




def format_drug_df(drug_path):
    """
    Format and process the drug dataframe from a given CSV file path.

    Reads a CSV file, processes its content, extracts chem_name and other_name 
    from the 'inhibitor' column, and returns the modified dataframe.

    Parameters
    ----------
    drug_path : str
        Path to the drug CSV file.

    Returns
    -------
    pd.DataFrame
        Formatted drug dataframe.
    """
    d_df = pd.read_csv(drug_path, index_col=None,sep="\t")
    d_df[['chem_name', 'other_name']] = d_df['DRUG'].str.extract(r'^(.*?)\s*(?:\((.+)\))?$')
    d_df["chem_name"] = d_df["chem_name"].str.replace('\s-\s', ':')
    return d_df

original_drug_file = "beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
original_drug_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
download_from_github(original_drug_url, original_drug_file)

additional_mapping_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx"
sample_mapping_file = "beataml_waves1to4_sample_mapping.xlsx"
download_from_github(additional_mapping_url, sample_mapping_file)

supplementary_url = 'https://ars.els-cdn.com/content/image/1-s2.0-S1535610822003129-mmc2.xlsx'
supplimentary_file = '1-s2.0-S1535610822003129-mmc2.xlsx'
download_from_github(supplementary_url, supplimentary_file)
updated_raw_drug_file = "beatAML_drug_raw.tsv"


generate_raw_drug_file(original_drug_file,sample_mapping_file, updated_raw_drug_file,supplimentary_file)

df = format_drug_df(updated_raw_drug_file)

chem_list = pd.concat([df.chem_name,df.other_name]).unique().tolist()

pubchem_retrieval.update_dataframe_and_write_tsv(unique_names=chem_list, output_filename="drugs.tsv")


