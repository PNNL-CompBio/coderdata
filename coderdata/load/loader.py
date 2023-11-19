# coderdata/load/loader.py

import pandas as pd
import gzip
import os



class DatasetLoader:
    def __init__(self, dataset_type, data_directory="."):
        """
        Load datasets of a specific type into predefined attributes of this class instance.

        Each attribute will be a Pandas DataFrame corresponding to a file with the dataset prefix.
        Attributes include transcriptomics, proteomics, mutations, etc.

        Parameters
        ----------
        dataset_type : str
            The type of dataset to load (e.g., 'hcmi', 'beataml').
        data_directory : str
            The directory where dataset files are stored.
        """
        # Initialize attributes as empty DataFrames
        self.transcriptomics = pd.DataFrame()
        self.proteomics = pd.DataFrame()
        self.mutations = pd.DataFrame()
        self.copy_number = pd.DataFrame()
        self.miRNA = pd.DataFrame()
        self.samples = pd.DataFrame()
        self.drugs = pd.DataFrame()
        self.experiments = pd.DataFrame()
        self.methylation = pd.DataFrame()
        self.metabolomics = pd.DataFrame()
        self.full = pd.DataFrame()
        self.load_datasets(dataset_type, data_directory)

    def load_datasets(self, dataset_type, data_directory):
        for file_name in os.listdir(data_directory):
            if file_name.startswith(dataset_type) and (file_name.endswith(('.csv', '.tsv', '.csv.gz', '.tsv.gz'))):
                file_path = os.path.join(data_directory, file_name)
#                 dataset_name = os.path.splitext(os.path.splitext(file_name)[0])[0].split('_')[1]
                dataset_name = file_name.split('_', 1)[1].split('.')[0]
                print("Loading:",dataset_name)
                if hasattr(self, dataset_name):
                    dataframe = self.load_file(file_path)
                    setattr(self, dataset_name, dataframe)

    @staticmethod
    def load_file(file_path):
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as file:
                return pd.read_csv(file, delimiter=DatasetLoader.determine_delimiter(file_path))
        else:
            return pd.read_csv(file_path, delimiter=DatasetLoader.determine_delimiter(file_path))

    @staticmethod
    def determine_delimiter(file_path):
        return '\t' if file_path.endswith('.tsv') or file_path.endswith('.tsv.gz') else ','

    def merge_datasets(self, drop_na=False):
        """
        Merge all non-empty available datasets into a single DataFrame.
        The 'samples' attribute is merged last, after removing duplicates based on 'improve_sample_id'.
        Merges on common keys ('improve_sample_id', 'entrez_id', 'improve_drug_id') where available.
        """
        # Gather all non-empty DataFrame attributes
        dataframes = {attr: getattr(self, attr) for attr in dir(self)
                      if isinstance(getattr(self, attr), pd.DataFrame)
                      and not getattr(self, attr).empty
                      and not attr.startswith('__')
                      and attr != 'full'}

        # Save the 'samples' DataFrame and remove it from the initial merging process
        samples_df = dataframes.pop('samples', pd.DataFrame())

        # List of potential merge keys
        merge_keys = ['improve_sample_id', 'entrez_id', 'improve_drug_id']
        print("Merging Datasets into full.")
        # Merge non-samples dataframes
        for attr, df in dataframes.items():
            if self.full.empty:
                self.full = df
            else:
                common_keys = [key for key in merge_keys if key in df.columns and key in self.full.columns]
                if common_keys:
                    self.full = pd.merge(self.full, df, on=common_keys, how='outer')

        # Merge 'samples' DataFrame last
        if not samples_df.empty:
            samples_df = samples_df.drop_duplicates(subset=['improve_sample_id'])
            common_keys = [key for key in merge_keys if key in samples_df.columns and key in self.full.columns]
            if common_keys:
                self.full = pd.merge(self.full, samples_df, on=common_keys, how='outer')


        merge_suffixes = ['_x', '_y', '_z']  # Add more if needed
        for base_col in set(col.rsplit('_', 1)[0] for col in self.full.columns if any(col.endswith(suffix) for suffix in merge_suffixes)):
            # Identify all suffix columns for the base column
            suffix_cols = [col for col in self.full.columns if col.startswith(base_col) and any(col.endswith(suffix) for suffix in merge_suffixes)]

            # Start with the first suffix column
            merged_col = self.full[suffix_cols[0]].copy()

            # Merge remaining suffix columns
            for col in suffix_cols[1:]:
                merged_col = merged_col.fillna(self.full[col])

            # Replace the original column (if it exists) with the merged column, else create a new column
            if base_col in self.full.columns:
                self.full[base_col] = self.full[base_col].fillna(merged_col)
            else:
                self.full[base_col] = merged_col

            # Drop the suffix columns
            self.full.drop(suffix_cols, axis=1, inplace=True)
