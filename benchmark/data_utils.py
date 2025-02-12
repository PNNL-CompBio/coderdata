import os
from typing import Union, List
import urllib
import pandas as pd
from rdkit import Chem
import os
import pandas as pd
import numpy as np
#------------------
# 1. download data
#------------------

# version = 'benchmark-data-pilot1'
# ftp_dir = f'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/{version}/csa_data/raw_data/'
# version = 'benchmark-data-imp-2023'
# ftp_dir = f'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/{version}/csa_data/'

# candle_data_dict = {
#     'ccle_candle': "CCLE",
#     'ctrpv2_candle':"CTRPv2",
#     'gdscv1_candle':"GDSCv1",
#     'gdscv2_candle':"GDSCv2",
#     'gcsi_candle': "gCSI"}

# class Downloader:
#     def __init__(self, version):
#         self.version = version
#         if version == 'benchmark-data-pilot1':
#             self.ftp_dir = f'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/{version}/csa_data/raw_data/'
#         elif version == 'benchmark-data-imp-2023':
#             self.ftp_dir = f'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/{version}/csa_data/'

#     def download_candle_data(self, data_type="CCLE", split_id=100, data_dest='Data/'):
#         self.download_candle_split_data(data_type=data_type, split_id=split_id, data_dest=data_dest)
#         self.download_candle_resp_data(data_dest=data_dest)
#         self.download_candle_gexp_data(data_dest=data_dest)
#         self.download_candle_mut_data(data_dest=data_dest)
#         self.download_candle_smiles_data(data_dest=data_dest)
#         self.download_candle_drug_ecfp4_data(data_dest=data_dest)
#         self.download_landmark_genes(data_dest=data_dest)


#     def download_deepttc_vocabs(self, data_dest='Data/'):

#         src_dir = 'https://raw.githubusercontent.com/jianglikun/DeepTTC/main/ESPF/' 
#         fnames = ['drug_codes_chembl_freq_1500.txt', 'subword_units_map_chembl_freq_1500.csv']

#         for fname in fnames:
#             src = os.path.join(src_dir, fname)
#             dest = os.path.join(data_dest, fname)
#             if not os.path.exists(dest):
#                 urllib.request.urlretrieve(src, dest)



#     def download_candle_split_data(self, data_type="CCLE", split_id=0, data_dest='Data/'):
#         print(f'downloading {data_type} split {split_id} data')


#         split_src = os.path.join(self.ftp_dir, 'splits')
#         train_split_name = f'{data_type}_split_{split_id}_train.txt'
#         val_split_name = f'{data_type}_split_{split_id}_val.txt'
#         test_split_name = f'{data_type}_split_{split_id}_test.txt'

        
#         # download split data
#         for file in [train_split_name, val_split_name, test_split_name]:
#             src = os.path.join(split_src, file)
#             dest = os.path.join(data_dest, file)

#             if not os.path.exists(dest):
#                 urllib.request.urlretrieve(src, dest)

#     def download_candle_resp_data(self, data_dest='Data/'):
#         # ftp_dir = 'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/benchmark-data-pilot1/csa_data/raw_data/'

#         print('downloading response data')
#         # download response data
#         resp_name = 'response.tsv'
#         if self.version=='benchmark-data-imp-2023':
#             resp_name = 'response.txt'

#         src = os.path.join(self.ftp_dir, 'y_data', resp_name)
#         dest = os.path.join(data_dest, resp_name)

#         if not os.path.exists(dest):
#             urllib.request.urlretrieve(src, dest)

#     def download_candle_gexp_data(self, data_dest='Data/'):
#         print('downloading expression data')
#         gexp_name = 'cancer_gene_expression.tsv'
#         if self.version=='benchmark-data-imp-2023':
#             gexp_name = 'cancer_gene_expression.txt'

#         src = os.path.join(self.ftp_dir, 'x_data', gexp_name)
#         dest = os.path.join(data_dest, gexp_name)

#         if not os.path.exists(dest):
#             urllib.request.urlretrieve(src, dest)

#     def download_candle_mut_data(self, data_dest='Data/'):
#         # ftp_dir = 'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/benchmark-data-pilot1/csa_data/raw_data/'
#         # gene mutation data
#         print('downloading mutation data')
#         gmut_name = 'cancer_mutation_count.tsv'
#         if self.version=='benchmark-data-imp-2023':
#             gmut_name = 'cancer_mutation_count.txt'

#         src = os.path.join(self.ftp_dir, 'x_data', gmut_name)
#         dest = os.path.join(data_dest, gmut_name)
#         if not os.path.exists(dest):
#             urllib.request.urlretrieve(src, dest)

#     def download_candle_smiles_data(self, data_dest='Data/'):
#         # ftp_dir = 'https://ftp.mcs.anl.gov/pub/candle/public/improve/benchmarks/single_drug_drp/benchmark-data-pilot1/csa_data/raw_data/'
#         # gene mutation data
#         print('downloading smiles data')

#         smiles_name = 'drug_SMILES.tsv'
#         if self.version=='benchmark-data-imp-2023':
#             smiles_name = 'drug_SMILES.txt'

#         src = os.path.join(self.ftp_dir, 'x_data', smiles_name)
#         dest = os.path.join(data_dest, smiles_name)
#         if not os.path.exists(dest):
#             urllib.request.urlretrieve(src, dest)

#     def download_candle_drug_ecfp4_data(self, data_dest='Data/'):
#         # gene mutation data
#         print('downloading drug_ecfp4 data')
#         name = 'drug_ecfp4_nbits512.tsv'
#         if self.version=='benchmark-data-imp-2023':
#             name = 'drug_ecfp4_512bit.txt'
            
#         src = os.path.join(self.ftp_dir, 'x_data', name)
#         dest = os.path.join(data_dest, name)
#         if not os.path.exists(dest):
#             urllib.request.urlretrieve(src, dest)

#     def download_landmark_genes(self, data_dest='Data/'):
#         urllib.request.urlretrieve('https://raw.githubusercontent.com/gihanpanapitiya/GraphDRP/to_candle/landmark_genes', data_dest+'/landmark_genes')


def add_smiles(smiles_df, df, metric):
    """
    Add smiles information to the given dataframe based on the improve_drug_id column.

    Args:
        smiles_df (pandas.DataFrame): DataFrame containing the smiles data.
        df (pandas.DataFrame): DataFrame to which smiles information needs to be added.
        metric (str): Name of the column containing the metric data.

    Returns:
        pandas.DataFrame: DataFrame with added smiles information.
    """
    data_smiles_df = pd.merge(df, smiles_df, on="improve_drug_id", how='left')
    data_smiles_df = data_smiles_df.dropna(subset=[metric])
    data_smiles_df = data_smiles_df[['improve_sample_id', 'canSMILES', 'improve_drug_id', metric]]
    data_smiles_df = data_smiles_df.drop_duplicates()
    data_smiles_df.dropna(inplace=True)
    data_smiles_df = data_smiles_df.reset_index(drop=True)

    return data_smiles_df


class DataProcessor:
    def __init__(self, version):
        self.version = version

    def load_drug_response_data(self, data_path, data_type="CCLE",
    split_id=100, split_type='train', response_type='ic50', sep="\t",
    dropna=True):
        """
        Returns datarame with cancer ids, drug ids, and drug response values. Samples
        from the original drug response file are filtered based on the specified
        sources. 

        Args:
            source (str or list of str): DRP source name (str) or multiple sources (list of strings)
            split(int or None): split id (int), None (load all samples)
            split_type (str or None): one of the following: 'train', 'val', 'test'
            y_col_name (str): name of drug response measure/score (e.g., AUC, IC50)

        Returns:
            pd.Dataframe: dataframe that contains drug response values
            columns: [source, improve_drug_id, improve_sample_id, response_type]
        """
        # TODO: at this point, this func implements the loading a single source
        y_file_path = os.path.join(data_path, 'response.tsv')
        # if self.version=='benchmark-data-imp-2023':
        #     y_file_path = os.path.join(data_path, 'response.txt')

        df = pd.read_csv(y_file_path, sep=sep)

        # import pdb; pdb.set_trace()
        if isinstance(split_id, int):
            # Get a subset of samples
            ids = self.load_split_file(data_path, data_type, split_id, split_type)
            df = df.loc[ids]
        else:
            # Get the full dataset for a given source
            df = df[df["source"].isin([data_type])]

        cols = ["source",
                "improve_drug_id",
                "improve_sample_id",
                response_type]
        df = df[cols]  # [source, drug id, cancer id, response]
        if dropna:
            df.dropna(axis=0, inplace=True)
        df = df.reset_index(drop=True)
        return df


    # def load_split_file(self,
    #     data_path: str,
    #     data_type: str,
    #     split_id: Union[int, None]=None,
    #     split_type: Union[str, List[str], None]=None) -> list:
    #     """
    #     Args:
    #         source (str): DRP source name (str)

    #     Returns:
    #         ids (list): list of id integers
    #     """
    #     if isinstance(split_type, str):
    #         split_type = [split_type]

    #     # Check if the split file exists and load
    #     ids = []
    #     for st in split_type:
    #         fpath = os.path.join(data_path, f"{data_type}_split_{split_id}_{st}.txt")
    #         # assert fpath.exists(), f"Splits file not found: {fpath}"
    #         ids_ = pd.read_csv(fpath, header=None)[0].tolist()
    #         ids.extend(ids_)
    #     return ids

#-----------------------------------
# 2. preprocess data to swnet format
#-----------------------------------
# def process_response_data(df_resp, response_type='ic50'):    
#     # df = pd.read_csv('response.tsv', sep='\t')
#     drd = df_resp[['improve_drug_id', 'improve_sample_id', response_type]]
#     drd.columns =['drug','cell_line','IC50']
#     # drd = drd.dropna(axis=0)
#     drd.reset_index(drop=True, inplace=True)
#     # drd.to_csv('tmp/Paccmann_MCA/Data/response.csv')
#     return drd

    def load_smiles_data(self,
        data_dir,
        sep: str="\t",
        verbose: bool=True) -> pd.DataFrame:
        """
        IMPROVE-specific func.
        Read smiles data.
        src_raw_data_dir : data dir where the raw DRP data is stored
        """
        
        smiles_path = os.path.join(data_dir, 'drug_SMILES.tsv')
        if self.version=='benchmark-data-imp-2023':
            smiles_path = os.path.join(data_dir, 'drug_SMILES.txt')

        df = pd.read_csv(smiles_path, sep=sep)

        # TODO: updated this after we update the data
        df.columns = ["improve_drug_id", "smiles"]

        if verbose:
            print(f"SMILES data: {df.shape}")
            # print(df.dtypes)
            # print(df.dtypes.value_counts())
        return df



    # def load_morgan_fingerprint_data(self,
    #     data_dir,
    #     sep: str="\t",
    #     verbose: bool=True) -> pd.DataFrame:
    #     """
    #     Return Morgan fingerprints data.
    #     """

    #     path = os.path.join(data_dir, 'drug_ecfp4_nbits512.tsv')
    #     if self.version=='benchmark-data-imp-2023':
    #         path = os.path.join(data_dir, 'drug_ecfp4_512bit.txt')

    #     df = pd.read_csv(path, sep=sep)
    #     df = df.set_index('improve_drug_id')
    #     return df

    # preprocess coderdata rna-seq data.


    # def load_gene_expression_data(self,
    #                         data_dir,
    #                         gene_system_identifier: Union[str, List[str]]="Entrez",
    #                         sep: str="\t",
    #                         verbose: bool=True) -> pd.DataFrame:
    #     """
    #     Returns gene expression data.

    #     Args:
    #         gene_system_identifier (str or list of str): gene identifier system to use
    #             options: "Entrez", "Gene_Symbol", "Ensembl", "all", or any list
    #                     combination of ["Entrez", "Gene_Symbol", "Ensembl"]

    #     Returns:
    #         pd.DataFrame: dataframe with the omic data
    #     """
    #     gene_expression_file_path = "cancer_gene_expression.tsv"
    #     canc_col_name= "improve_sample_id"
    #     # level_map encodes the relationship btw the column and gene identifier system
    #     level_map = {"Ensembl": 0, "Entrez": 1, "Gene_Symbol": 2}
    #     header = [i for i in range(len(level_map))]

    #     df = pd.read_csv(gene_expression_file_path, sep=sep, index_col=0, header=header)

    #     df.index.name = canc_col_name  # assign index name
    #     df = set_col_names_in_multilevel_dataframe(df, level_map, gene_system_identifier)
    #     if verbose:
    #         print(f"Gene expression data: {df.shape}")
    #         # print(df.dtypes)
    #         # print(df.dtypes.value_counts())
    #     return df
    
    def convert_long_to_wide_format(input_file_path: str, sep: str = '\t') -> None:
        """
        Converts a long-format gene expression data file to wide format.

        Args:
            input_file_path (str): Path to the input CSV file in long format.
            sep (str, optional): Separator used in the output file. Defaults to '\t'.
        """

        # Load the long-format gene expression data
        if input_file_path.endswith('.gz'):
            df_long = pd.read_csv(input_file_path, compression='gzip')
        else:
            df_long = pd.read_csv(input_file_path)
        
        # Check if 'Sanger' exists in the 'source' column
        if 'Sanger' in df_long['source'].values:
            # Apply log transformation to 'transcriptomics' values where the 'source' is 'Sanger'
            df_long.loc[df_long['source'] == 'Sanger', 'transcriptomics'] = np.log1p(df_long.loc[df_long['source'] == 'Sanger', 'transcriptomics'])

        # Pivot the DataFrame to wide format
        df_wide = df_long.pivot_table(index='improve_sample_id', columns='entrez_id', values='transcriptomics', aggfunc='first')

        # Reset the index to make `improve_sample_id` a column again, and fill NaN values with 0
        df_wide.reset_index(inplace=True)
        df_wide.fillna(0, inplace=True)

        # Save the wide-format DataFrame to the specified output file in the "shared_input" directory
        input_file_name = os.path.basename(input_file_path)
        # output_file_name = input_file_name.replace(".csv.gz", "_wide.tsv")
        output_file_path = os.path.join("./shared_input", input_file_name+"_wide.tsv")
        df_wide.to_csv(output_file_path, sep=sep, index=False)

    def CNV_convert_long_to_wide_format(input_file_path: str, sep: str = '\t') -> None:
        # Load the long-format copy number data
        if input_file_path.endswith('.gz'):
            df_long = pd.read_csv(input_file_path, compression='gzip')
        else:
            df_long = pd.read_csv(input_file_path)

        # if the df_long has a column named 'copy_number', change column name to 'CNV'
        if 'copy_number' in df_long.columns:
            df_long.rename(columns={'copy_number': 'CNV'}, inplace=True)

        # Pivot the DataFrame to wide format
        df_wide = df_long.pivot_table(index='improve_sample_id', columns='entrez_id', values='CNV', aggfunc='first')

        # Reset the index to make `improve_sample_id` a column again, and fill NaN values with 0
        df_wide.reset_index(inplace=True)
        df_wide.fillna(0, inplace=True)

        # Process to remove decimal points from column names (if any)
        new_columns = {}
        for col in df_wide.columns:
            try:
                # Attempt to convert column names to float then to int and back to string
                # This will be done for columns that are not 'improve_sample_id'
                if col != 'improve_sample_id':
                    new_col = str(int(float(col)))
                    new_columns[col] = new_col
            except ValueError:
                # In case the column name is not a float, keep it as is
                pass

        # Update the DataFrame columns if there are any to update
        if new_columns:
            df_wide.rename(columns=new_columns, inplace=True)

        # Save the wide-format DataFrame to the specified output file in the "shared_input" directory
        input_file_name = os.path.basename(input_file_path)
        output_file_name = input_file_name.replace(".csv.gz", "_wide.tsv").replace(".csv", "_wide.tsv")
        output_file_path = os.path.join("./shared_input", output_file_name)
        df_wide.to_csv(output_file_path, sep=sep, index=False)

    # def load_cell_mutation_data(self,
    #             data_dir,
    #                             gene_system_identifier: Union[str, List[str]]="Gene_Symbol",
    #                             sep: str="\t", verbose: bool=True) -> pd.DataFrame:
    #     """
    #     Returns gene expression data.

    #     Args:
    #         gene_system_identifier (str or list of str): gene identifier system to use
    #             options: "Entrez", "Gene_Symbol", "Ensembl", "all", or any list
    #                     combination of ["Entrez", "Gene_Symbol", "Ensembl"]

    #     Returns:
    #         pd.DataFrame: dataframe with the omic data
    #     """
    #     cell_mutation_file_path = os.path.join(data_dir, 'cancer_mutation_count.tsv')
    #     if self.version=='benchmark-data-imp-2023':
    #         cell_mutation_file_path = os.path.join(data_dir, 'cancer_mutation_count.txt')
    #     canc_col_name= "improve_sample_id"
    #     # level_map encodes the relationship btw the column and gene identifier system
    #     level_map = {"Ensembl": 0, "Entrez": 1, "Gene_Symbol": 2}
    #     header = [i for i in range(len(level_map))]

    #     df = pd.read_csv(cell_mutation_file_path, sep=sep, index_col=0, header=header)

    #     df.index.name = canc_col_name  # assign index name
    #     df = set_col_names_in_multilevel_dataframe(df, level_map, gene_system_identifier)
    #     if verbose:
    #         print(f"cell mutation data: {df.shape}")
    #         # print(df.dtypes)
    #         # print(df.dtypes.value_counts())
    #     return df


    # def load_landmark_genes(self, data_path):
    #     genes = pd.read_csv(os.path.join(data_path, 'landmark_genes'), header=None)
    #     genes = genes.values.ravel().tolist()
    #     return genes

def set_col_names_in_multilevel_dataframe(
    df: pd.DataFrame,
    level_map: dict,
    gene_system_identifier: Union[str, List[str]]="Gene_Symbol") -> pd.DataFrame:
    """ Util function that supports loading of the omic data files.
    Returns the input dataframe with the multi-level column names renamed as
    specified by the gene_system_identifier arg.

    Args:
        df (pd.DataFrame): omics dataframe
        level_map (dict): encodes the column level and the corresponding identifier systems
        gene_system_identifier (str or list of str): gene identifier system to use
            options: "Entrez", "Gene_Symbol", "Ensembl", "all", or any list
                     combination of ["Entrez", "Gene_Symbol", "Ensembl"]

    Returns:
        pd.DataFrame: the input dataframe with the specified multi-level column names
    """
    df = df.copy()

    level_names = list(level_map.keys())
    level_values = list(level_map.values())
    n_levels = len(level_names)
    
    if isinstance(gene_system_identifier, list) and len(gene_system_identifier) == 1:
        gene_system_identifier = gene_system_identifier[0]

    # print(gene_system_identifier)
    # import pdb; pdb.set_trace()
    if isinstance(gene_system_identifier, str):
        if gene_system_identifier == "all":
            df.columns = df.columns.rename(level_names, level=level_values)  # assign multi-level col names
        else:
            df.columns = df.columns.get_level_values(level_map[gene_system_identifier])  # retian specific column level
    else:
        assert len(gene_system_identifier) <= n_levels, f"'gene_system_identifier' can't contain more than {n_levels} items."
        set_diff = list(set(gene_system_identifier).difference(set(level_names)))
        assert len(set_diff) == 0, f"Passed unknown gene identifiers: {set_diff}"
        kk = {i: level_map[i] for i in level_map if i in gene_system_identifier}
        # print(list(kk.keys()))
        # print(list(kk.values()))
        df.columns = df.columns.rename(list(kk.keys()), level=kk.values())  # assign multi-level col names
        drop_levels = list(set(level_map.values()).difference(set(kk.values())))
        df = df.droplevel(level=drop_levels, axis=1)
    return df


def remove_smiles_with_noneighbor_frags(smiles_df):

    remove_smiles=[]
    for i in smiles_df.index:
        smiles = smiles_df.loc[i, 'smiles']
        has_atoms_wothout_neighbors = check_for_atoms_without_neighbors(smiles)
        if has_atoms_wothout_neighbors:
            remove_smiles.append(smiles)

    smiles_df = smiles_df[~smiles_df.smiles.isin(remove_smiles)]
    smiles_df.dropna(inplace=True)
    smiles_df.reset_index(drop=True, inplace=True)

    return smiles_df    

def check_for_atoms_without_neighbors(smiles):

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    frags = Chem.GetMolFrags(mol, asMols=True)
    frag_atoms = [i.GetNumAtoms() for i in frags]
    has_atoms_wothout_neighbors = any([i==1 for i in frag_atoms])

    
    return has_atoms_wothout_neighbors



# def load_generic_expression_data(file_path,
#                         gene_system_identifier = "Gene_Symbol",
#                         sep: str="\t",
#                         verbose: bool=True) -> pd.DataFrame:
#     """
#     Returns gene expression data.

#     Args:
#         gene_system_identifier (str or list of str): gene identifier system to use
#             options: "Entrez", "Gene_Symbol", "Ensembl", "all", or any list
#                     combination of ["Entrez", "Gene_Symbol", "Ensembl"]

#     Returns:
#         pd.DataFrame: dataframe with the omic data
#     """

    
#     canc_col_name= "improve_sample_id"
#     # level_map encodes the relationship btw the column and gene identifier system
#     level_map = {"Ensembl": 0, "Entrez": 1, "Gene_Symbol": 2}
#     header = [i for i in range(len(level_map))]

#     df = pd.read_csv(file_path, sep=sep, index_col=0, header=header)

#     df.index.name = canc_col_name  # assign index name
#     df = set_col_names_in_multilevel_dataframe(df, level_map, gene_system_identifier)
#     return df

# function that will average auc values for the same improve_sample_id and improve_drug_id
def average_auc(df):
    # This will convert non-numeric values to NaN
    df['auc'] = pd.to_numeric(df['auc'], errors='coerce')
    df['auc'] = df.groupby(['improve_sample_id', 'improve_drug_id'])['auc'].transform('mean')
    df = df.drop_duplicates(subset=['improve_sample_id', 'improve_drug_id'])
    return df

def average_ic50(df):
    # This will convert non-numeric values to NaN
    df['ic50'] = pd.to_numeric(df['ic50'], errors='coerce')
    df['ic50'] = df.groupby(['improve_sample_id', 'improve_drug_id'])['ic50'].transform('mean')
    df = df.drop_duplicates(subset=['improve_sample_id', 'improve_drug_id'])
    return df


# updated version of coderdata
def average_dose_response_value(df):
    # This will convert non-numeric values to NaN
    df['dose_response_value'] = pd.to_numeric(df['dose_response_value'], errors='coerce')
    df['dose_response_value'] = df.groupby(['improve_sample_id', 'improve_drug_id'])['dose_response_value'].transform('mean')
    df = df.drop_duplicates(subset=['improve_sample_id', 'improve_drug_id'])
    return df

def filter_exp_data(exp, study_description, dose_response_metric):
    # List available studies from the 'study' column
    available_studies = exp['study'].unique()
    
    # Check if the user-provided study description is in the available studies
    if study_description not in available_studies:
        # If not, return an error with a list of available studies
        raise ValueError(f"Study description '{study_description}' not found. Available studies: {', '.join(available_studies)}")
    
    # Subset based on the study description
    filtered_exp = exp[exp['study'].str.contains(study_description)]
    
    # Remove any NaN in 'dose_response_metric' [update]
    filtered_exp = filtered_exp.dropna(subset=['dose_response_metric'])
    
    # Now, filter based on dose response metric
    # This assumes 'metric' is the column where dose_response_metric values are stored
    if dose_response_metric not in filtered_exp['dose_response_metric'].unique():
        # If the dose_response_metric is not found in the filtered dataset, return an error
        available_metrics = ', '.join(filtered_exp['dose_response_metric'].unique())
        raise ValueError(f"Dose response metric '{dose_response_metric}' not found. Available metrics in '{study_description}': {available_metrics}")
    
    # Further filtering by dose response metric
    final_filtered_exp = filtered_exp[filtered_exp['dose_response_metric'] == dose_response_metric]
    
    return final_filtered_exp
