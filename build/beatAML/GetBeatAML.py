import synapseclient
import pandas as pd
import os
import wget
from synapseclient import Project, Folder, File, Link
import requests
import numpy as np
import subprocess
import argparse
import time

# def download_from_github(raw_url, save_path):
#     """
#     Download a file from a raw GitHub URL and save to the specified path.

#     Parameters
#     ----------
#     raw_url : str
#         The raw GitHub URL pointing to the file to be downloaded.
#     save_path : str
#         The local path where the downloaded file will be saved.

#     Returns
#     -------
#     None
#     """
#     response = requests.get(raw_url)
#     with open(save_path, 'wb') as f:
#         f.write(response.content)
#     return

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

def remove_sixth_from_last_quote(s):
    """
    Modify an incorrectly formatted file. 
    """
    reversed_s = s[::-1]  # Reverse the string
    count = 0
    index_to_remove = None
    for i, char in enumerate(reversed_s):
        if char == '"':
            count += 1
        if count == 6:
            index_to_remove = i
            break
    # Remove the 6th occurrence of the quote character
    if index_to_remove is not None:
        reversed_s = reversed_s[:index_to_remove] + reversed_s[index_to_remove+1:]
    return reversed_s[::-1]  # Reverse again to obtain the original order

def modify_patient_file(file_path):
    """
    Modify an incorrectly formatted file. 
    """
    with open(file_path, 'r', encoding='ISO-8859-1') as f:
        lines = f.readlines()
    # Modify the 110th line (0-based index is 109)
    lines[109] = remove_sixth_from_last_quote(lines[109])
    with open(file_path, 'w', encoding='utf-8') as f:
        f.writelines(lines)
        
def generate_samples_file(prev_samples_path):
    """
    Generate samples file by reading, processing and merging multiple source files.

    This function reads two Excel files, processes the contained data,
    and generates a combined sample file in CSV format.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the combined samples data.
    """

    # Read the beataml_waves1to4_sample_mapping.xlsx file
    beataml_samples = pd.read_excel('beataml_waves1to4_sample_mapping.xlsx')
    additional_samples = beataml_samples[~beataml_samples['rna_control'].isin(['No', '']) & pd.notna(beataml_samples['rna_control'])]
    
    # Map rna_control to common_name and retain columns LabId, dbgap_rnaseq_sample
    additional_samples = additional_samples[['labId', 'dbgap_rnaseq_sample', 'rna_control']]

    # Read the other file
    mapped_info = pd.read_excel('1-s2.0-S1535610822003129-mmc2.xlsx', skiprows=1)

    # Merge with the original samples on dbgap_rnaseq_sample
    original_samples = mapped_info.merge(beataml_samples, on="dbgap_rnaseq_sample", how="inner")

    # Concatenate the original samples and the additional samples
    full_samples = pd.concat([original_samples, additional_samples], ignore_index=True)

    full_samples = full_samples[["labId", "specificDxAtInclusion", "specimenType","rna_control"]]
    full_samples['common_name'] = np.where(full_samples['specimenType'].isna(), full_samples['rna_control'], full_samples['specimenType'])
    full_samples.drop(columns=["specimenType","rna_control"], inplace=True)
    full_samples.rename(columns={"labId": "other_id"}, inplace=True)
    full_samples.rename(columns={"specificDxAtInclusion": "other_names"}, inplace=True)
    full_samples['other_names'] = full_samples['other_names'].fillna('Control')
    full_samples["cancer_type"] = "Acute Myeloid Leukaemia"
    full_samples["model_type"] = "ex vivo"
    full_samples["other_id_source"] = "beatAML"
    full_samples.drop_duplicates(subset='other_id', keep='first', inplace=True)
    
    prot_sample_file = 'PNNL_clinical_summary_12_08_2021_updated_OS_4patients_02_28_2022.txt'
    modify_patient_file(prot_sample_file)
    prot_samples = pd.read_csv(prot_sample_file, delimiter=' ', skipinitialspace=True, encoding='ISO-8859-1', quotechar='"')
    prot_samples = prot_samples[["labId","specificDxAtInclusion","specimenType"]]
    prot_samples.rename(columns={"labId": "other_id"}, inplace=True)
    prot_samples.rename(columns={"specificDxAtInclusion": "other_names"}, inplace=True)
    prot_samples.rename(columns={"specimenType": "common_name"}, inplace=True)
    prot_samples["cancer_type"] = "Acute Myeloid Leukaemia"
    prot_samples["model_type"] = "ex vivo"
    prot_samples["other_id_source"] = "beatAML"
    
    all_samples = pd.concat([prot_samples, full_samples])
    all_samples['species'] = 'Homo sapiens (Human)'
    if prev_samples_path == "":
        maxval = 0
    else:
        maxval = max(pd.read_csv(prev_samples_path).improve_sample_id)
    mapping = {labId: i for i, labId in enumerate(all_samples['other_id'].unique(), start=(int(maxval)+1))}
    all_samples['improve_sample_id'] = all_samples['other_id'].map(mapping)
    all_samples.insert(1, 'improve_sample_id', all_samples.pop('improve_sample_id'))
    all_samples.to_csv("/tmp/beataml_samples.csv", index=False)
    return all_samples

    
def retrieve_drug_info(compound_name):
    """
    Retrieve detailed information for a given compound name using the PubChem API.

    Parameters
    ----------
    compound_name : str
        The name of the compound for which detailed information is needed.

    Returns
    -------
    tuple
        A tuple containing pubchem_id, CanonicalSMILES, IsomericSMILES, InChIKey,
        MolecularFormula, and MolecularWeight of the compound.
    """
    if pd.isna(compound_name):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    ##limit is 1 call per 5 seconds. add in wait call.
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight/JSON"
    response = requests.get(url)

    if response.status_code != 200:
        print(response.text)
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    data = response.json()
    if "PropertyTable" in data:
        properties = data["PropertyTable"]["Properties"][0]
        pubchem_id = properties.get('CID',np.nan)
        canSMILES = properties.get("CanonicalSMILES", np.nan)
        # isoSMILES = properties.get("IsomericSMILES", np.nan)
        InChIKey = properties.get("InChIKey", np.nan)
        formula = properties.get("MolecularFormula", np.nan)
        weight = properties.get("MolecularWeight", np.nan)

        return pubchem_id, canSMILES, InChIKey, formula, weight
    else:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    

def update_dataframe_with_pubchem(d_df):
    """
    Update the provided dataframe with drug information from PubChem.

    For each chem_name in the dataframe, retrieve drug details using the PubChem API 
    and update the dataframe columns with the fetched information.

    Parameters
    ----------
    d_df : pd.DataFrame
        The dataframe containing a 'chem_name' column which will be used to fetch information.

    Returns
    -------
    pd.DataFrame
        Updated dataframe with drug details from PubChem.
    """
    # Get unique chem_name values and retrieve their data
    chem_names = d_df['chem_name'].dropna().unique()
    chem_data_dict = {}
    for name in chem_names:
        print("Attempting to call pubchem API for chem_name: ", name)
        chem_data_dict[name] = retrieve_drug_info(name)
        time.sleep(0.2)
    failed_chem_names = {k for k, v in chem_data_dict.items() if all(pd.isna(val) for val in v)}
    other_names = d_df[d_df['chem_name'].isin(failed_chem_names)]['other_name'].dropna().unique()
    other_data_dict = {}
    for name in other_names:
        print("Attempting to call pubchem API for other_name: ", name)
        other_data_dict[name] = retrieve_drug_info(name)
        time.sleep(0.2)

    # Combine both dictionaries for easy lookup
    data_dict = {**chem_data_dict, **other_data_dict}

    #print(data_dict)
#    print(data_dict['isoSMILES'])
    # Update the DataFrame using the data dictionary
    for idx, row in d_df.iterrows():
        if row['chem_name'] in data_dict and not all(pd.isna(val) for val in data_dict[row['chem_name']]):
            values = data_dict[row['chem_name']]
        else:
            values = data_dict.get(row['other_name'], (np.nan, np.nan, np.nan, np.nan, np.nan))

        d_df.at[idx, 'pubchem_id'] = values[0]
        d_df.at[idx, "canSMILES"] = values[1]
        # d_df.at[idx, "isoSMILES"] = values[2]
        d_df.at[idx, "InChIKey"] = values[2]
        d_df.at[idx, "formula"] = values[3]
        d_df.at[idx, "weight"] = values[4]
    
    return d_df

def merge_drug_info(d_df,drug_map):
    """
    Merge drug information from a given drug mapping dataframe to the main drug dataframe.

    Parameters
    ----------
    d_df : pd.DataFrame
        Main drug dataframe containing drug-related columns.
    drug_map : pd.DataFrame
        Mapping dataframe containing drug information and the column 'canSMILES'.

    Returns
    -------
    pd.DataFrame
        The merged dataframe containing combined drug information.
    """
    # print(d_df['isoSMILES'].dtype, drug_map['isoSMILES'].dtype)
    d_df['canSMILES'] = d_df['canSMILES'].astype(str)
    drug_map['canSMILES'] = drug_map['canSMILES'].astype(str)
    result_df = d_df.merge(drug_map[['canSMILES', 'improve_drug_id']], on='canSMILES', how='left')
    return result_df

def format_drug_map(drug_map_path):
    """
    Format and clean up the drug mapping file.

    Reads a drug map file, removes duplicates based on the 'canSMILES' column,
    and returns the cleaned dataframe.

    Parameters
    ----------
    drug_map_path : str
        Path to the drug mapping file.

    Returns
    -------
    pd.DataFrame
        Formatted and cleaned drug mapping dataframe.
    """
    if drug_map_path:
        drug_map = pd.read_csv(drug_map_path, sep = "\t")
        drug_map = drug_map.drop_duplicates(subset='canSMILES', keep='first')
    else:
        drug_map = pd.DataFrame(columns=[
            'improve_drug_id', 'chem_name', 'pubchem_id',
            'canSMILES', 'InChIKey', 'formula', 'weight'
        ])
    return drug_map

#Drug Response
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
    d_df[['chem_name', 'other_name']] = d_df['inhibitor'].str.extract(r'^(.*?)\s*(?:\((.+)\))?$')
    d_df["chem_name"] = d_df["chem_name"].str.replace('\s-\s', ':',regex=True)
    d_df['chem_name'] = [a.lower() for a in d_df['chem_name']]
    return d_df

def add_improve_id(previous_df, new_df):
    """
    Add 'improve_drug_id' to the new dataframe based on unique 'canSMILES' not present in the previous dataframe.

    Parameters
    ----------
    previous_df : pd.DataFrame
        Dataframe containing previously mapped drug ids.
    new_df : pd.DataFrame
        Dataframe where new ids need to be added.

    Returns
    -------
    pd.DataFrame
        New dataframe with 'improve_drug_id' added.
    """
    if not previous_df.empty and 'improve_drug_id' in previous_df.columns:
        id_list = [int(val.replace('SMI_', '')) for val in previous_df['improve_drug_id'].tolist() if pd.notnull(val) and val.startswith('SMI_')]
        max_id = max(id_list) if id_list else 0
    else:
        max_id = 0
    # Identify canSMILES in the new dataframe that don't exist in the old dataframe
    unique_new_smiles = set(new_df['canSMILES']) - set(previous_df['canSMILES'])
    # Identify rows in the new dataframe with canSMILES that are unique and where improve_drug_id is NaN
    mask = (new_df['canSMILES'].isin(unique_new_smiles)) & (new_df['improve_drug_id'].isna())
    id_map = {}
    for smiles in unique_new_smiles:
        max_id += 1
        id_map[smiles] = f"SMI_{max_id}"
    # Apply the mapping to the new dataframe for rows with unique canSMILES and NaN improve_drug_id
    new_df.loc[mask, 'improve_drug_id'] = new_df['canSMILES'].map(id_map)
    return new_df


def map_exp_to_improve(exp_path):#df,improve_map_file):
    """
    Map 'sample_id' in the given dataframe to 'improve_sample_id' using a provided mapping file.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing 'sample_id' which needs mapping.
    improve_map_file : str
        Path to the CSV file containing the mapping between 'sample_id' and 'improve_sample_id'.

    Returns
    -------
    pd.DataFrame
        Mapped dataframe with 'improve_sample_id' added and 'sample_id' removed.
    """
    mapped_df = pd.read_csv(exp_path,sep='\t')
    mapped_df['source'] = 'synapse'
    mapped_df['study'] = 'BeatAML'
    return mapped_df



def map_and_combine(df, data_type, entrez_map_file, improve_map_file, map_file=None):
    """
    Map the sample ids, combine dataframes and output the merged dataframe based on the given data type.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing data which needs to be mapped and combined.
    data_type : str
        String specifying the type of data ('expression' or 'mutation').
    entrez_map_file : str
        Path to the file containing mapping from gene name to entrez gene id.
    improve_map_file : str
        Path to the file containing mapping between 'sample_id' and 'improve_sample_id'.
    map_file : str, optional
        Additional mapping file required based on the data type. Default is None.

    Returns
    -------
    pd.DataFrame
        Mapped and combined dataframe.
    """
    # Load mapping files
    genes = pd.read_csv(entrez_map_file)          # Map gene_name to entrez_id
    improve = pd.read_csv(improve_map_file)        # Map sample_id to improve_sample_id
    if map_file:
        mapped_ids = pd.read_excel(map_file)
    
    # Process each dataframe based on its data_type
    if data_type == "transcriptomics":
        df['Gene'] = df['Gene'].str.replace(r'\.\d+$', '', regex=True)
        mapped_df = df.merge(genes, left_on='Gene', right_on='other_id', how='left').reindex(
                        columns=['transcriptomics', 'entrez_id', "sample_id","Gene"])
        mapped_df = mapped_df.merge(mapped_ids[['dbgap_rnaseq_sample', 'labId']], 
                         left_on='sample_id', 
                         right_on='dbgap_rnaseq_sample', 
                         how='left')
        mapped_df.drop(columns=['sample_id', 'dbgap_rnaseq_sample'], inplace=True)

        mapped_df.rename(columns={'labId': 'sample_id'}, inplace=True)
    
    elif data_type == "mutations":
        df = df[['dbgap_sample_id','hgvsc', 'hgvsp', 'gene', 'variant_classification','t_vaf', 'refseq', 'symbol']]
        mapped_df = df.merge(genes, left_on='symbol', right_on='gene_symbol', how='left').reindex(
                        columns=['hgvsc', 'entrez_id', "dbgap_sample_id","variant_classification"])
        mapped_df = mapped_df.merge(mapped_ids[['dbgap_dnaseq_sample', 'labId']], 
                         left_on='dbgap_sample_id', 
                         right_on='dbgap_dnaseq_sample', 
                         how='left')

        mapped_df.rename(columns={"hgvsc": "mutation"}, inplace=True)
        mapped_df.rename(columns={"labId": "sample_id"}, inplace=True)
        mapped_df.rename(columns={"Entrez_Gene_Id": "entrez_id"}, inplace=True)

        variant_mapping = {
            'frameshift_variant': 'Frameshift_Variant',
            'missense_variant': 'Missense_Mutation',
            'stop_gained': 'Nonsense_Mutation',
            'inframe_deletion': 'In_Frame_Del',
            'protein_altering_variant': 'Protein_Altering_Variant',
            'splice_acceptor_variant': 'Splice_Site',
            'splice_donor_variant': 'Splice_Site',
            'start_lost': 'Start_Codon_Del',
            'inframe_insertion': 'In_Frame_Ins',
            'stop_lost': 'Nonstop_Mutation'
        }

        mapped_df['variant_classification'] = mapped_df['variant_classification'].map(variant_mapping)

    elif data_type == "proteomics":
        mapped_ids['sampleID'] = mapped_ids['sampleID'].str.split('_').apply(lambda x: x[2])
        mapped_ids.rename(columns={"sampleID": "id"}, inplace=True)
        mapped_ids = mapped_ids.reindex(columns=['labId', "id"])        
        df = df.merge(mapped_ids,
                      left_on='id',
                      right_on='id',
                      how='left')
        df = df.reindex(columns=["Protein", "proteomics", "labId"])
        mapped_df = df.merge(genes,
                             left_on='Protein',
                             right_on='gene_symbol',
                             how='left').reindex(
            columns=['proteomics', 'entrez_id', 'labId'])
        mapped_df.rename(
            columns={
                "labId": "sample_id",
            },
            inplace=True
        )
        
    mapped_df = pd.merge(mapped_df, improve[['other_id', 'improve_sample_id']], 
                         left_on='sample_id', 
                         right_on='other_id',
                         how='left')
    mapped_df.insert(0, 'improve_sample_id', mapped_df.pop('improve_sample_id'))
    
    # Replace NaNs, round values, and convert to integers for specified columns
    columns_to_convert = ['improve_sample_id', 'entrez_id']
    mapped_df[columns_to_convert] = mapped_df[columns_to_convert].fillna(0).round().astype('int32')
    
    mapped_df['source'] = 'synapse'
    mapped_df['study'] = 'BeatAML'
    mapped_df =mapped_df.drop_duplicates()

    final_dataframe = mapped_df.dropna()
    return final_dataframe


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
    drug_mod['DRUG'] = [a.lower() for a in drug_mod['DRUG']]
    
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
    
    drug_mod['source'] = cohort_merge['cohort'].combine_first(cohort_merge_secondary['cohort'])
    drug_mod['study'] = "BeatAML"
    drug_mod = drug_mod[['CELL', 'DRUG', 'DOSE', 'GROWTH', 'source', 'study']]
    drug_mod = drug_mod.sort_values(by=['DRUG', 'CELL'])
    drug_mod = drug_mod.groupby(['DRUG', 'CELL']).filter(lambda x: len(x) >= 5)
    drug_mod['time'] = 72
    drug_mod['time_unit'] ='hrs'
    drug_mod.to_csv(updated_raw_drug_file, index=False, sep="\t")
    return

def generate_drug_list(drug_map_path,drug_path):
    '''
    generates mapping of AML files to drugs
    '''
    # Drug and Experiment Data
    print("Starting Drug Data")
    drug_map = format_drug_map(drug_map_path) ##read in original/prior drugs in db
    d_df = format_drug_df(drug_path) ##format new drug data from beataml
    d_df = update_dataframe_with_pubchem(d_df)
    d_res = merge_drug_info(d_df, drug_map)
    d_res = add_improve_id(drug_map, d_res)
    #Drug Data
    #print(d_res)
    drug_res = d_res[["improve_drug_id","chem_name","pubchem_id","formula","weight","InChIKey","canSMILES"]]
    drug_res = drug_res.drop_duplicates()
    drug_res.to_csv("/tmp/beataml_drugs.tsv",sep="\t", index=False)

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script handles all aspects of the beat AML data but was designed to be fun in four modes (samples, drugs, omics, exp). If no argument is provided it will try to run all four at once.')
    parser.add_argument('-t', '--token', type=str, help='Synapse Token')
    ##the next three arguments determine what we'll do

    parser.add_argument('-s', '--samples', action = 'store_true', help='Only generate samples, requires previous samples',default=False)
    parser.add_argument('-p', '--prevSamples', nargs='?',type=str, default='', const='', help='Use this to provide previous sample file, will run sample file generation')
    
    parser.add_argument('-d', '--drugs',action='store_true', default=False,help='Query drugs only, requires drug file')
    parser.add_argument('-r', '--drugFile',nargs='?',type=str, default='', const='',help='Path to existing drugs.tsv file to query')
    
    parser.add_argument('-o', '--omics',action='store_true',default=False,help='Set this flag to query omics, requires current samples')
    parser.add_argument('-c', '--curSamples', type=str, help='Add path if you want to generate data')
    parser.add_argument('-g', '--genes',type=str, help='Path to gene file, required for omics processing')
    
    parser.add_argument('-e', '--exp', action='store_true', default=False,help='Set this to generate dose response curves. requires drug file and sample file')



    args = parser.parse_args()

    #######
    

    print("Logging into Synapse")
    syn = synapseclient.Synapse()
    PAT = args.token
    syn.login(authToken=PAT)
    # drug response: syn51674470
    # gene: syn25714248
    # mutatons: syn32533104
    # rna: syn32529921
    # clinical summary: syn26642974
    # Proteomics: syn26427390
    entity_ids = [
#         'syn51674470', 
        'syn25714248', 
#         'syn32533104', 
#         'syn32529921', 
        'syn26642974',
        'syn26427390',
        'syn64126458',
        'syn64126462',
        'syn64126463',
        'syn64126464',
        'syn64126468'
    ]
    print("Downloading Files from Synapse")
    for entity_id in entity_ids:
        entity = syn.get(entity_id, downloadLocation='.')
        
    # Download required files. Files in github repo also required.
    #gene_url = "https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47"
    #entrez_map_file = retrieve_figshare_data(gene_url)

    # additional_mapping_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx"
    sample_mapping_file = "beataml_waves1to4_sample_mapping.xlsx"
    # download_from_github(additional_mapping_url, sample_mapping_file)

    # supplementary_url = 'https://ars.els-cdn.com/content/image/1-s2.0-S1535610822003129-mmc2.xlsx'
    supplimentary_file = '1-s2.0-S1535610822003129-mmc2.xlsx'
    # download_from_github(supplementary_url, supplimentary_file)
    
    
    if args.samples:
        if args.prevSamples is None or args.prevSamples=='':
            print("No Previous Samples file was found. Data will not align with other datasets. Use ONLY for testing purposes.")
        else:
            print("Previous Samples File Provided. Running BeatAML Sample File Generation")
        #Generate Samples File
        generate_samples_file(args.prevSamples)
    if args.drugs:
        if args.drugFile is None or args.drugFile=='':
            print("Prior Drug File not provided. Data will not align with other datasets. Use ONLY for testing purposes.")
        else:
            print("Drug File Provided. Proceeding with build.")
        original_drug_file = "beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
        # original_drug_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
        # download_from_github(original_drug_url, original_drug_file)
        generate_drug_list(args.drugFile, original_drug_file) 
    if args.omics:
        if args.genes is None or args.curSamples is None:
            print('Cannot process omics without sample mapping and gene mapping files')
            exit()
        else:
            improve_map_file = args.curSamples
            transcriptomics_file = "beataml_waves1to4_counts_dbgap.txt" #"beataml_waves1to4_norm_exp_dbgap.txt" ##this is the wrong file, these are the normalize values
            # transcriptomics_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_counts_dbgap.txt" #"https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt"
            # download_from_github(transcriptomics_url, transcriptomics_file)
            
            mutations_file = "beataml_wes_wv1to4_mutations_dbgap.txt"
            # mutations_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_wes_wv1to4_mutations_dbgap.txt"
            # download_from_github(mutations_url, mutations_file)
            
            mutation_map_file = "beataml_waves1to4_sample_mapping.xlsx"
            # mutation_map_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx"
            # download_from_github(mutation_map_url, mutation_map_file)
            # New Transcriptomics Data
            print("Starting Transcriptomics Data")
            ##first run conversion tool
            os.system("python tpmFromCounts.py --counts "+transcriptomics_file)
            
            
            t_df = pd.read_csv('tpm_'+transcriptomics_file, sep = '\t')
            t_df = t_df.reset_index().rename(columns={'stable_id': 'Gene'})
            t_df = pd.melt(t_df, id_vars=['Gene'], var_name='sample_id', value_name='transcriptomics')
            print(improve_map_file)
            t_df = map_and_combine(t_df, "transcriptomics", args.genes, improve_map_file, sample_mapping_file)
            t_df = t_df[t_df.entrez_id.notna()]
            t_df = t_df[["improve_sample_id","transcriptomics","entrez_id","source","study"]].drop_duplicates()
            t_df.to_csv("/tmp/beataml_transcriptomics.csv.gz",index=False,compression='gzip')

            # New Proteomics Data
            print("Starting Proteomics Data")
            proteomics_map = "Data Available for Proteomic Samples.xlsx"
            p_df = pd.read_csv("ptrc_ex10_crosstab_global_gene_corrected.txt", sep = '\t')
            p_df = p_df.reset_index().rename(columns={'index': 'Protein'})
            p_df = pd.melt(p_df, id_vars=['Protein'], var_name='id', value_name='proteomics')
            p_df = map_and_combine(p_df, "proteomics", args.genes, improve_map_file, proteomics_map)
            p_df = p_df[["improve_sample_id","proteomics","entrez_id","source","study"]]
            p_df.to_csv("/tmp/beataml_proteomics.csv.gz",index=False,compression='gzip')
        
            # New Mutation Data
            print("Starting Mutation Data")
            m_df = pd.read_csv(mutations_file, sep = '\t')
            
            m_df = map_and_combine(m_df, "mutations", args.genes,improve_map_file, mutation_map_file)
            m_df = m_df[["improve_sample_id","mutation", "entrez_id","variant_classification","source","study"]]
            m_df.to_csv("/tmp/beataml_mutations.csv.gz",index=False,compression='gzip')
        
    if args.exp:
        if args.curSamples is None or args.drugFile is None:
            print("Cannot run curve fitting without drug mapping and sample mapping")
            exit()
        else:
            imp_samp_map = pd.read_csv(args.curSamples)
            imp_drug_map = pd.read_csv(args.drugFile,sep='\t')
            original_drug_file = "beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
            # original_drug_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"    
            # Generate Raw Drugs File to use in Curve fitting algorithm
            # download_from_github(original_drug_url, original_drug_file)
             # Experiment Data
            updated_raw_drug_file = "beatAML_drug_raw.tsv"
            generate_raw_drug_file(original_drug_file,sample_mapping_file, updated_raw_drug_file,supplimentary_file)
            d_df = pd.read_csv(updated_raw_drug_file,sep='\t')
            d_res = d_df.rename(columns={"CELL":"other_id","AUC":"fit_auc",'DRUG':'chem_name'})
            d_res = d_res.merge(imp_samp_map, on='other_id')
            d_res = d_res.merge(imp_drug_map,on='chem_name')
            d_res = d_res.rename(columns = {'improve_drug_id':'Drug'}) ## stupid but we have to change aks later
            d_res.to_csv(updated_raw_drug_file,sep='\t')

            print("Starting Curve Fitting Algorithm")
            # Run Curve fitting algorithm from scripts directory.
            # Note the file path to fit_curve.py may need to be changed.
            command = ['python3', 'fit_curve.py' ,'--input', 'beatAML_drug_raw.tsv', '--output', 'beatAML_drug_processed.tsv', '--beataml']
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                print("Curve Fitting executed successfully!")
            else:
                print("Curve fitting failed.")
                print("Out:", result.stdout)
                print("Error:", result.stderr)
            print("Starting Experiment Data")
            drug_path = "beatAML_drug_processed.tsv.0"
            exp_res = map_exp_to_improve(drug_path)
            exp_res.to_csv("/tmp/beataml_experiments.tsv", index=False, sep='\t')
          
#        print("Finished Pipeline")
    
