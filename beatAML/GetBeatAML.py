import synapseclient
import pandas as pd
import os
import wget
import requests
import numpy as np
from synapseclient import Project, Folder, File, Link

syn = synapseclient.Synapse()
PAT = "ENTER YOUR PERSONAL ACCESS TOKEN HERE"
syn.login(authToken=PAT)

def download_from_github(raw_url, save_path):
    response = requests.get(raw_url)
    with open(save_path, 'wb') as f:
        f.write(response.content)
    return

def retrieve_figshare_data(url):
    files_0 = os.listdir()
    wget.download(url)
    files_1 = os.listdir()
    new_file = str(next(iter(set(files_1) - set(files_0))))
    return new_file


def map_and_combine(df, data_type, entrez_map_file, improve_map_file, sample_map=None):
    """
    Maps and combines dataframes based on their data_type. It then merges 
    the final dataframe with provided metadata.
    """
    # Initialize the list to hold mapped dataframes

    # Load mapping files
    samples = sample_map
    genes = pd.read_csv(entrez_map_file)          # Map gene_name to entrez_id
    
    improve = pd.read_csv(improve_map_file)        # Map sample_id to improve_sample_id

    # Process each dataframe based on its data_type
#     for df in dataframe_list:
    if data_type == "transcriptomics":
        mapped_df = df.merge(genes, left_on='Gene', right_on='gene_symbol', how='left').reindex(
                        columns=['transcriptomics', 'entrez_id', "sample_id"])
#         mapped_df.drop(columns=['gene_id', 'gene_name', 'gene_type'], inplace=True)

    elif data_type == "mutations":
#         mapped_df = df.reindex(columns=['Entrez_Gene_Id', 'symbol', 'hgvsc'])
        mapped_df = df.merge(genes, left_on='symbol', right_on='gene_symbol', how='left').reindex(
                        columns=['hgvsc', 'entrez_id', "alt_ID"])
        mapped_df.rename(columns={"hgvsc": "mutations"}, inplace=True)
        mapped_df.rename(columns={"alt_ID": "sample_id"}, inplace=True)
        mapped_df.rename(columns={"Entrez_Gene_Id": "entrez_id"}, inplace=True)

    elif data_type == "proteomics":

        # Update 'sampleID' in mapped_ids dataframe
        sample_map['sampleID'] = sample_map['sampleID'].str.split('_').apply(lambda x: x[2])

        # Rename 'sampleID' to 'id'
        sample_map.rename(columns={"sampleID": "id"}, inplace=True)
#         print("sample_map", sample_map)
        # Reindex the columns of mapped_ids dataframe
        sample_map = mapped_ids.reindex(columns=['labId', "id"])        
        
        # Merge p_df with mapped_ids
        df = df.merge(mapped_ids, left_on='id', right_on='id', how='left')

        # Reindex the columns of p_df
        df = df.reindex(columns=["Protein", "proteomics", "labId"])

        # Merge df with genes and reindex + rename columns
        mapped_df = df.merge(genes, left_on='Protein', right_on='gene_symbol', how='left').reindex(
            columns=['proteomics', 'entrez_id', 'labId'])

        mapped_df.rename(
            columns={
                "labId": "sample_id",
            },
            inplace=True
        )
    # Add improve ID then Drop the sample_id and other_id columns
    mapped_df = pd.merge(mapped_df, improve[['other_id', 'improve_sample_id']], left_on='sample_id', right_on='other_id', how='left')
    mapped_df.drop(columns=['sample_id', 'other_id'], inplace=True)
    mapped_df.insert(0, 'improve_sample_id', mapped_df.pop('improve_sample_id'))
    mapped_df['source'] = 'synapse'
    mapped_df['study'] = 'BeatAML'

    # Concatenate the list of dataframes into one final dataframe
    final_dataframe = mapped_df
    return final_dataframe


def retrieve_drug_info(compound_name):
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
    # Get unique chem_name values and retrieve their data
    chem_names = d_df['chem_name'].dropna().unique()
    chem_data_dict = {}
    for name in chem_names:
        print("Attempting to call pubchem API for chem_name: ", name)
        chem_data_dict[name] = retrieve_drug_info(name)

    # Identify rows whose chem_name failed and retrieve the other_name for those rows
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
    result_df = d_df.merge(drug_map[['isoSMILES', 'improve_drug_id']], on='isoSMILES', how='left')
    result_df = result_df[["sample_id","improve_drug_id","chem_name","auc","formula","weight","canSMILES","inchikey","isoSMILES"]]
    return result_df

def format_drug_map(drug_map_path):
    drug_map = pd.read_csv(drug_map_path, sep = "\t")
    drug_map = drug_map.drop_duplicates(subset='isoSMILES', keep='first')
    return drug_map

#Drug Response
def format_drug_df(drug_path):
    d_df = pd.read_csv(drug_path, index_col=None)
    d_df.drop("Unnamed: 0",axis='columns',inplace=True)
    d_df[['chem_name', 'other_name']] = d_df['inhibitor'].str.extract(r'^(.*?)\s*(?:\((.+)\))?$')
    d_df["chem_name"] = d_df["chem_name"].str.replace('\s-\s', ':')
    return d_df


def add_improve_id(previous_df, new_df):
    
    # Extract the maximum value of improve_drug_id from the previous dataframe
    # Assuming the format is SMI_number. We get rid of 'SMI_' and then convert to int.
    max_id = max([int(val.replace('SMI_', '')) for val in previous_df['improve_drug_id'].tolist() if pd.notnull(val) and val.startswith('SMI_')])

    # Identify isoSMILES in the new dataframe that don't exist in the old dataframe
    unique_new_smiles = set(new_df['isoSMILES']) - set(previous_df['isoSMILES'])
    
    # Identify rows in the new dataframe with isoSMILES that are unique and where improve_drug_id is NaN
    mask = (new_df['isoSMILES'].isin(unique_new_smiles)) & (new_df['improve_drug_id'].isna())
    
    # Create a mapping of unique new isoSMILES to improve_drug_id
    id_map = {}
    for smiles in unique_new_smiles:
        max_id += 1
        id_map[smiles] = f"SMI_{max_id}"
    
    # Apply the mapping to the new dataframe for rows with unique isoSMILES and NaN improve_drug_id
    new_df.loc[mask, 'improve_drug_id'] = new_df['isoSMILES'].map(id_map)
    
    return new_df

def map_exp_to_improve(df,improve_map_file):
    improve = pd.read_csv(improve_map_file)        # Map sample_id to improve_sample_id
    mapped_df = pd.merge(df, improve[['other_id', 'improve_sample_id']], left_on='sample_id', right_on='other_id', how='left')
    mapped_df.drop(columns=['sample_id', 'other_id'], inplace=True)
    mapped_df.insert(0, 'improve_sample_id', mapped_df.pop('improve_sample_id'))
    mapped_df['source'] = 'synapse'
    mapped_df['study'] = 'BeatAML'
    mapped_df['ic50'] = np.nan
    mapped_df['ec50'] = np.nan
    mapped_df['ec50se'] = np.nan
    mapped_df['einf'] = np.nan
    mapped_df['hs'] = np.nan
    mapped_df['aac1'] = np.nan
    mapped_df['auc1'] = np.nan
    mapped_df['dss1'] = np.nan
    return mapped_df



if __name__ == "__main__":
    
    # Download required files. Files in github repo also required.
    gene_url = "https://figshare.com/ndownloader/files/40576109?private_link=525f7777039f4610ef47"
    entrez_map_file = retrieve_figshare_data(gene_url)

    samples_url = "https://figshare.com/ndownloader/files/42289053?private_link=f05259d19b8c34a9a53d"
    improve_map_file = retrieve_figshare_data(samples_url)

    # Transcriptomics Data
    t_df = pd.read_csv("BeatAML_Waves1to4_RNA_data_normalized.txt", sep = '\t')
    t_df = t_df.reset_index().rename(columns={'index': 'Gene'})
    t_df = pd.melt(t_df, id_vars=['Gene'], var_name='sample_id', value_name='transcriptomics')
    t_df = map_and_combine(t_df, "transcriptomics", entrez_map_file, improve_map_file)
    print("Writing transcriptomics.csv")
    t_df.to_csv("transcriptomics.csv", index=False)

    # Mutation Data
    p_df = pd.read_csv("ptrc_ex10_crosstab_global_gene_corrected.txt", sep = '\t')
    p_df = p_df.reset_index().rename(columns={'index': 'Protein'})
    p_df = pd.melt(p_df, id_vars=['Protein'], var_name='id', value_name='proteomics')
    mapped_ids = pd.read_csv("Proteomic_sample_map_beatAML.csv")
    p_df = map_and_combine(p_df, "proteomics", entrez_map_file, improve_map_file, mapped_ids)
    print("Writing proteomics.csv")
    p_df.to_csv("proteomics.csv", index=False)

    # Mutation Data
    m_df = pd.read_csv("beataml_wes_wv1to4_mutations_177_samples.txt", sep = '\t')
    m_df = map_and_combine(m_df, "mutations", entrez_map_file, improve_map_file)
    print("Writing mutations.csv")
    m_df.to_csv("mutations.csv", index=False)


    # Drug and Experiment Data
    drug_path = "drug_response.csv"
    drug_map_path = "drugs_map.tsv.gz"
    drug_url = "https://raw.githubusercontent.com/PNNL-CompBio/candleDataProcessing/main/cell_line/drugs_by_structure.tsv.gz"

    download_from_github(drug_url, drug_map_path)
    drug_map = format_drug_map(drug_map_path)
    d_df = format_drug_df(drug_path)
    d_df = update_dataframe_with_pubchem(d_df)
    d_res = merge_drug_info(d_df, drug_map)
    d_res = add_improve_id(drug_map, d_res)

    #Drug Data
    drug_res = d_res[["improve_drug_id","chem_name","formula","weight","inchikey","canSMILES","isoSMILES"]]
    drug_res.rename(columns={"inchikey": "inCHIKey"}, inplace=True)

    drug_res
    print("Writing drugs.tsv")
    drug_res.to_csv("drugs.tsv",sep="\t", index=False)

    # Experiment Data
    exp_res = d_res[["sample_id","improve_drug_id","auc"]]
    exp_res = map_exp_to_improve(exp_res,improve_map_file)
    print("Writing experiments.csv")
    exp_res.to_csv("experiments.csv", index=False)