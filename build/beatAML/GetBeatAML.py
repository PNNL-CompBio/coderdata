import synapseclient
import pandas as pd
import os
import wget
from synapseclient import Project, Folder, File, Link
import requests
import numpy as np
import subprocess
import argparse


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
    full_samples["cancer_type"] = "ACUTE MYELOID LEUKAEMIA"
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
    prot_samples["cancer_type"] = "ACUTE MYELOID LEUKAEMIA"
    prot_samples["model_type"] = "ex vivo"
    prot_samples["other_id_source"] = "beatAML"    
    
    all_samples = pd.concat([prot_samples, full_samples])
    maxval = max(pd.read_csv(prev_samples_path).improve_sample_id)
    mapping = {labId: i for i, labId in enumerate(all_samples['other_id'].unique(), start=(int(maxval)+1))}
    all_samples['improve_sample_id'] = all_samples['other_id'].map(mapping)
    all_samples.insert(1, 'improve_sample_id', all_samples.pop('improve_sample_id'))
    all_samples.to_csv("beataml_samples.csv", index=False)
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
        A tuple containing CanonicalSMILES, IsomericSMILES, InChIKey,
        MolecularFormula, and MolecularWeight of the compound.
    """
    if pd.isna(compound_name):
        return np.nan, np.nan, np.nan, np.nan, np.nan
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight/JSON"
    response = requests.get(url)

    if response.status_code != 200:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    
    data = response.json()
    if "PropertyTable" in data:
        properties = data["PropertyTable"]["Properties"][0]
        canSMILES = properties.get("CanonicalSMILES", np.nan)
        isoSMILES = properties.get("IsomericSMILES", np.nan)
        inchikey = properties.get("InChIKey", np.nan)
        formula = properties.get("MolecularFormula", np.nan)
        weight = properties.get("MolecularWeight", np.nan)

        return canSMILES, isoSMILES, inchikey, formula, weight
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
    failed_chem_names = {k for k, v in chem_data_dict.items() if all(pd.isna(val) for val in v)}
    other_names = d_df[d_df['chem_name'].isin(failed_chem_names)]['other_name'].dropna().unique()
    other_data_dict = {}
    for name in other_names:
        print("Attempting to call pubchem API for other_name: ", name)
        other_data_dict[name] = retrieve_drug_info(name)

    # Combine both dictionaries for easy lookup
    data_dict = {**chem_data_dict, **other_data_dict}

    # Update the DataFrame using the data dictionary
    for idx, row in d_df.iterrows():
        if row['chem_name'] in data_dict and not all(pd.isna(val) for val in data_dict[row['chem_name']]):
            values = data_dict[row['chem_name']]
        else:
            values = data_dict.get(row['other_name'], (np.nan, np.nan, np.nan, np.nan, np.nan))
        
        d_df.at[idx, "canSMILES"] = values[0]
        d_df.at[idx, "isoSMILES"] = values[1]
        d_df.at[idx, "inchikey"] = values[2]
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
        Mapping dataframe containing drug information and the column 'isoSMILES'.

    Returns
    -------
    pd.DataFrame
        The merged dataframe containing combined drug information.
    """
    result_df = d_df.merge(drug_map[['isoSMILES', 'improve_drug_id']], on='isoSMILES', how='left')
    return result_df

def format_drug_map(drug_map_path):
    """
    Format and clean up the drug mapping file.

    Reads a drug map file, removes duplicates based on the 'isoSMILES' column,
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
    drug_map = pd.read_csv(drug_map_path, sep = "\t")
    drug_map = drug_map.drop_duplicates(subset='isoSMILES', keep='first')
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
    d_df[['chem_name', 'other_name']] = d_df['DRUG'].str.extract(r'^(.*?)\s*(?:\((.+)\))?$')
    d_df["chem_name"] = d_df["chem_name"].str.replace('\s-\s', ':')
    return d_df

def add_improve_id(previous_df, new_df):
    """
    Add 'improve_drug_id' to the new dataframe based on unique 'isoSMILES' not present in the previous dataframe.

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
    max_id = max([int(val.replace('SMI_', '')) for val in previous_df['improve_drug_id'].tolist() if pd.notnull(val) and val.startswith('SMI_')])
    # Identify isoSMILES in the new dataframe that don't exist in the old dataframe
    unique_new_smiles = set(new_df['isoSMILES']) - set(previous_df['isoSMILES'])
    # Identify rows in the new dataframe with isoSMILES that are unique and where improve_drug_id is NaN
    mask = (new_df['isoSMILES'].isin(unique_new_smiles)) & (new_df['improve_drug_id'].isna())
    id_map = {}
    for smiles in unique_new_smiles:
        max_id += 1
        id_map[smiles] = f"SMI_{max_id}"
    # Apply the mapping to the new dataframe for rows with unique isoSMILES and NaN improve_drug_id
    new_df.loc[mask, 'improve_drug_id'] = new_df['isoSMILES'].map(id_map)
    return new_df


def map_exp_to_improve(df,improve_map_file):
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
    improve = pd.read_csv(improve_map_file)        # Map sample_id to improve_sample_id
    mapped_df = pd.merge(df, improve[['other_id', 'improve_sample_id']], left_on='sample_id', right_on='other_id', how='left')
    mapped_df.drop(columns=['sample_id', 'other_id'], inplace=True)
    mapped_df.insert(0, 'improve_sample_id', mapped_df.pop('improve_sample_id'))
    mapped_df['source'] = 'synapse'
    mapped_df['study'] = 'BeatAML'
    mapped_df= mapped_df.rename(columns={'IC50':'ic50',
                              'EC50':'ec50',
                             'EC50se':'ec50se',
                             'Einf':'einf',
                             'HS':'hs',
                             'AAC1':'aac1',
                             'AUC1':'auc1',
                             'DSS1':'dss1',
                             'R2fit':'r2fit'
                             }
                    )
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
        mapped_df = df.merge(genes, left_on='Gene', right_on='gene_symbol', how='left').reindex(
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

        mapped_df.rename(columns={"hgvsc": "mutations"}, inplace=True)
        mapped_df.rename(columns={"labId": "sample_id"}, inplace=True)
        mapped_df.rename(columns={"Entrez_Gene_Id": "entrez_id"}, inplace=True)
        
    elif data_type == "mutations":
        df = df[['dbgap_sample_id','hgvsc', 'hgvsp', 'gene', 'variant_classification','t_vaf', 'refseq', 'symbol']]
        mapped_df = df.merge(genes, left_on='symbol', right_on='gene_symbol', how='left').reindex(
                        columns=['hgvsc', 'entrez_id', "dbgap_sample_id","variant_classification"])


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
    mapped_df['source'] = 'synapse'
    mapped_df['study'] = 'BeatAML'

    final_dataframe = mapped_df
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers and a string.')
    parser.add_argument('-t', '--token', type=str, help='Synapse Token')
    args = parser.parse_args()
    
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
        'syn26427390'
    ]
    print("Downloading Files")
    for entity_id in entity_ids:
        entity = syn.get(entity_id, downloadLocation='.')
        
    proteomics_map = "Data Available for Proteomic Samples.xlsx"
    
    # Download required files. Files in github repo also required.
    gene_url = "https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47"
    entrez_map_file = retrieve_figshare_data(gene_url)

    additional_mapping_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx"
    sample_mapping_file = "beataml_waves1to4_sample_mapping.xlsx"
    download_from_github(additional_mapping_url, sample_mapping_file)
    
    original_drug_file = "beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
    original_drug_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt"
    download_from_github(original_drug_url, original_drug_file)
    
    updated_raw_drug_file = "beatAML_drug_raw.tsv"
    
    drug_path = "beatAML_drug_processed.tsv.0"
    drug_map_path = retrieve_figshare_data("https://figshare.com/ndownloader/files/43112314?private_link=0ea222d9bd461c756fb0")
    
    transcriptomics_file = "beataml_waves1to4_norm_exp_dbgap.txt"
    transcriptomics_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt"
    download_from_github(transcriptomics_url, transcriptomics_file)
    
    mutations_file = "beataml_wes_wv1to4_mutations_dbgap.txt"
    mutations_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_wes_wv1to4_mutations_dbgap.txt"
    download_from_github(mutations_url, mutations_file)
   
    mutation_map_file = "beataml_waves1to4_sample_mapping.xlsx"
    mutation_map_url = "https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx"
    download_from_github(mutation_map_url, mutation_map_file)
    
    supplementary_url = 'https://ars.els-cdn.com/content/image/1-s2.0-S1535610822003129-mmc2.xlsx'
    supplimentary_file = '1-s2.0-S1535610822003129-mmc2.xlsx'
    download_from_github(supplementary_url, supplimentary_file)
    
    prev_samples_path = "../hcmi/hcmi_samples.csv"
    
    #Generate Samples File
    generate_samples_file(prev_samples_path)
    improve_map_file = "beataml_samples.csv"
    
    print("Starting Raw Drug File Generation ")
    # Generate Raw Drugs File to use in Curve fitting algorithm
    generate_raw_drug_file(original_drug_file,sample_mapping_file, updated_raw_drug_file,supplimentary_file)

    print("Starting Curve Fitting Algorithm")
    # Run Curve fitting algorithm from scripts directory.
    # Note the file path to fit_curve.py may need to be changed.
    command = ['python', '../utils/fit_curve.py' ,'--input', 'beatAML_drug_raw.tsv', '--output', 'beatAML_drug_processed.tsv']
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        print("Curve Fitting executed successfully!")
    else:
        print("Curve fitting failed.")

    # New Transcriptomics Data
    print("Starting Transcriptomics Data")
    t_df = pd.read_csv(transcriptomics_file, sep = '\t')
    t_df.index = t_df.display_label
    t_df = t_df.iloc[:, 4:]
    t_df = t_df.reset_index().rename(columns={'display_label': 'Gene'})
    t_df = pd.melt(t_df, id_vars=['Gene'], var_name='sample_id', value_name='transcriptomics')
    t_df = map_and_combine(t_df, "transcriptomics", entrez_map_file, "beataml_samples.csv", sample_mapping_file)
    t_df = t_df[t_df.entrez_id.notna()]
    t_df = t_df[["improve_sample_id","transcriptomics","entrez_id","source","study"]]
    t_df.to_csv("beataml_transcriptomics.csv",index=False)

    # New Proteomics Data
    print("Starting Proteomics Data")
    p_df = pd.read_csv("ptrc_ex10_crosstab_global_gene_corrected.txt", sep = '\t')
    p_df = p_df.reset_index().rename(columns={'index': 'Protein'})
    p_df = pd.melt(p_df, id_vars=['Protein'], var_name='id', value_name='proteomics')
    p_df = map_and_combine(p_df, "proteomics", entrez_map_file, improve_map_file, proteomics_map)
    p_df = p_df[["improve_sample_id","proteomics","entrez_id","source","study"]]
    p_df.to_csv("beataml_proteomics.csv",index=False)
    
    # New Mutation Data
    print("Starting Mutation Data")
    m_df = pd.read_csv(mutations_file, sep = '\t')
    m_df = map_and_combine(m_df, "mutations", entrez_map_file, "beataml_samples.csv", mutation_map_file)
    m_df = m_df[["improve_sample_id","mutations", "entrez_id","variant_classification","source","study"]]
    m_df.to_csv("beataml_mutations.csv",index=False)
    
    # Drug and Experiment Data
    print("Starting Drug Data")
    drug_map = format_drug_map(drug_map_path)
    d_df = format_drug_df(drug_path)
    d_df = update_dataframe_with_pubchem(d_df)
    d_res = merge_drug_info(d_df, drug_map)
    d_res = add_improve_id(drug_map, d_res)
    #Drug Data
    drug_res = d_res[["improve_drug_id","chem_name","formula","weight","inchikey","canSMILES","isoSMILES"]]
    drug_res.rename(columns={"inchikey": "inCHIKey"}, inplace=True)
    drug_res.to_csv("beataml_drugs.tsv",sep="\t", index=False)
    
    print("Starting Experiment Data")
    # Experiment Data
    d_res = d_res.rename(columns={"CELL":"sample_id","AUC":"auc"})
    exp_res = map_exp_to_improve(d_res,"beataml_samples.csv")
    exp_res = exp_res[["source","improve_sample_id","improve_drug_id","study","auc","ic50","ec50","ec50se","r2fit","einf","hs","aac1","auc1","dss1"]]
    exp_res.to_csv("beataml_experiments.csv", index=False)
    print("Finished Pipeline")