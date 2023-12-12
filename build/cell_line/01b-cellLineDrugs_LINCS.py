import pandas as pd
import requests
import numpy as np

"""
Overview:
1. Get LINCS drug information from GEO
2. Get current drug data from FigShare
3. Get PubChem information for new drugs
4. Create new IMPROVE IDs for new drugs
5. Generate new file with drug information
"""

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

def update_LINCS_dataframe_with_pubchem(d_df):
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
    data_dict = {}
    for name in chem_names:
        print("Attempting to call pubchem API for chem_name: ", name)
        data_dict[name] = retrieve_drug_info(name)

    # Update the DataFrame using the data dictionary
    for idx, row in d_df.iterrows():
        if row['chem_name'] in data_dict and not all(pd.isna(val) for val in data_dict[row['chem_name']]):
            values = data_dict[row['chem_name']]

        d_df.at[idx, "canSMILES"] = values[0]
        d_df.at[idx, "isoSMILES"] = values[1]
        d_df.at[idx, "InChIKey"] = values[2]
        d_df.at[idx, "formula"] = values[3]
        d_df.at[idx, "weight"] = values[4]
    
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

if __name__ == "__main__":
    # 1. Get LINCS drug information from GEO
    LINCS_drugs_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5Fpert%5Finfo.txt.gz"
    LINCS_drugs = pd.read_table(LINCS_drugs_url, delimiter = "\t") # cols: pert_id, canonical_smiles, inchi_key, pert_iname, pert_type
    LINCS_drugs = LINCS_drugs[LINCS_drugs['pert_type'] == "trt_cp"] # removes DMSO controls, leaving just drug treatments
    LINCS_drugs = LINCS_drugs[["pert_iname","inchi_key"]].drop_duplicates()
    
    # 2. Get current drug data from FigShare
    drugs_url = "https://figshare.com/ndownloader/files/42357210"
    drugs = pd.read_table(drugs_url, compression="gzip")
    
    # 3. Get PubChem information for new drugs
    LINCS_drugs.rename(columns={"pert_iname": "chem_name"}, inplace=True)
    LINCS_drugs.rename(columns={"inchi_key": "InChIKey"}, inplace=True)
    new_drugs = LINCS_drugs[~LINCS_drugs['InChIKey'].isin(drugs['InChIKey'])]
    old_drugs = drugs[drugs['InChIKey'].isin(LINCS_drugs['inchi_key'])]
    new_drugs['improve_drug_id'] = np.nan
    new_drugs = update_LINCS_dataframe_with_pubchem(new_drugs)
    
    # 4. Create new IMPROVE IDs for new drugs
    new_drugs = add_improve_id(drugs, new_drugs)
    more_old_drugs = drugs[drugs['isoSMILES'].isin(new_drugs['isoSMILES'])]
    new_drugs = new_drugs.dropna() # remove drugs not given new IDs because isoSMILES matched existing drug
    
    # 5. Generate new file with drug information
    new_drugs = new_drugs[["formula","weight","canSMILES","isoSMILES","InChIKey","chem_name","pubchem_id","improve_drug_id"]]
    LINCS_drugs = pd.concat([old_drugs, new_drugs]).drop_duplicates()
    LINCS_drugs.to_csv("lincs_drugs.tsv", sep="\t", index=False)
