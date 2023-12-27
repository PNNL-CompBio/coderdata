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
            The type of dataset to load (e.g., 'hcmi', 'beataml', 'cptac', 'cell_line', 'lincs').
        data_directory : str
            The directory where dataset files are stored.
        """
        # Initialize attributes as empty DataFrames
        self.transcriptomics = pd.DataFrame()
        self.proteomics = pd.DataFrame()
        self.mutations = pd.DataFrame()
        self.copy_number = pd.DataFrame()
        self.samples = pd.DataFrame()
        self.drugs = pd.DataFrame()
        self.mirna = pd.DataFrame()
        self.experiments = pd.DataFrame()
        self.methylation = pd.DataFrame()
        self.metabolomics = pd.DataFrame()
        self.full = pd.DataFrame()
        self.dataset_type = dataset_type
        self.data_directory = data_directory
        self.load_datasets(dataset_type, data_directory)
        self.data_format_params = {
        'transcriptomics': ('improve_sample_id', 'entrez_id', 'transcriptomics'),
        'proteomics': ('improve_sample_id', 'entrez_id', 'proteomics'),
        'mutations': ('improve_sample_id', 'entrez_id', 'mutations'),
        'copy_number': ('improve_sample_id', 'entrez_id', 'copy_number'),
        'methylation': ('improve_sample_id', 'entrez_id', 'methylation')
    }

    def load_datasets(self, dataset_type, data_directory):
        for file_name in os.listdir(data_directory):
            if file_name.startswith(dataset_type) and (file_name.endswith(('.csv', '.tsv', '.csv.gz', '.tsv.gz'))):
                file_path = os.path.join(data_directory, file_name)
                dataset_name = file_name[len(dataset_type):].split('.')[0].strip('_')
                print("Loading:", dataset_name)
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
            
    def reformat_dataset(self, *args):
        """
        Reformat a given dataset or all datasets to either 'long' or 'wide' format.
        
        Parameters
        ----------
        *args : variable length argument list
            Can be either one argument (format) to reformat all datasets, or two arguments
            (dataset_name, format) to reformat a specific dataset.
        """
        # Validate and parse arguments
        if len(args) == 1:
            format = args[0]
            datasets = self.data_format_params.keys()
        elif len(args) == 2:
            dataset_name, format = args
            datasets = [dataset_name]
        else:
            raise ValueError("Invalid number of arguments. Expecting one or two arguments.")

        if format not in ['long', 'wide']:
            raise ValueError("Format must be 'long' or 'wide'.")

        # Reformat datasets
        for dataset_name in datasets:
            if hasattr(self, dataset_name):
                dataset = getattr(self, dataset_name)
                if not dataset.empty and dataset_name in self.data_format_params:
                    id_vars, var_name, value_name = self.data_format_params[dataset_name]

                    # Check if the dataset is already in the desired format
                    if format == 'wide' and self.is_wide_format(dataset, id_vars, var_name):
                        print(f"Dataset '{dataset_name}' is already in wide format.")
                        continue
                    elif format == 'long' and self.is_long_format(dataset, id_vars, var_name, value_name):
                        print(f"Dataset '{dataset_name}' is already in long format.")
                        continue

                    # Convert to the desired format
                    if format == 'wide':
                        setattr(self, dataset_name, self.to_wide_format(dataset, id_vars, var_name, value_name))
                        print(f"{dataset_name} successfully converted to wide format")
                    elif format == 'long':
                        setattr(self, dataset_name, self.to_long_format(dataset, id_vars, var_name, value_name))
                        print(f"{dataset_name} successfully converted to long format")
            else:
                print(f"Dataset '{dataset_name}' is empty or already in desired format.")


    def to_wide_format(self, dataset, index, columns, values):
        # Reset index if needed
        if dataset.index.name == columns:
            dataset = dataset.reset_index()

        # Check if the values column is numeric and apply the appropriate aggregation method. 
        if pd.api.types.is_numeric_dtype(dataset[values]):
            dataset_agg = dataset.groupby([index, columns])[values].mean().reset_index()
        else:
            dataset_agg = dataset.groupby([index, columns])[values].first().reset_index()

        # Pivot the aggregated dataset
        wide_format = dataset_agg.pivot(index=index, columns=columns, values=values)
        return wide_format.reset_index()

    def to_long_format(self, dataset, id_vars, var_name, value_name):
        long_format = pd.melt(dataset, id_vars=id_vars, var_name=var_name, value_name=value_name)
        return long_format
    
    def is_wide_format(self, dataset, id_var, column_var):
        """
        Check if the dataset is in wide format.
        """
        return (dataset.columns.name == column_var and dataset.columns[0] == id_var)

    def is_long_format(self, dataset, id_vars, var_name, value_name):
        """
        Check if the dataset is in long format.
        """
        required_columns = set([id_vars, var_name, value_name])
        return required_columns.issubset(set(dataset.columns))
    
    def reload_datasets(self, dataset_name=None):
        """
        Reload a specific dataset or all datasets.

        Parameters
        ----------
        dataset_name : str, optional
            The name of the dataset to reload. If None, all datasets are reloaded.
        """
        if dataset_name:
            # Reload a specific dataset
            if hasattr(self, dataset_name):
                file_path = self._get_file_path(dataset_name)
                if file_path:
                    dataframe = self.load_file(file_path)
                    setattr(self, dataset_name, dataframe)
                    print(f"Dataset '{dataset_name}' reloaded.")
                else:
                    print(f"File for dataset '{dataset_name}' not found.")
            else:
                print(f"Dataset '{dataset_name}' does not exist.")
        else:
            # Reload all datasets
            for attr in self.data_format_params.keys():
                if hasattr(self, attr):
                    file_path = self._get_file_path(attr)
                    if file_path:
                        dataframe = self.load_file(file_path)
                        setattr(self, attr, dataframe)
                        print(f"Dataset '{attr}' reloaded.")
                    else:
                        print(f"File for dataset '{attr}' not available.")
                        
    def _get_file_path(self, dataset_name):
        """
        Utility method to get the file path for a given dataset name.
        This method assumes a specific naming convention and directory structure.
        """
        for file_name in os.listdir(self.data_directory):
            if file_name.startswith(self.dataset_type) and file_name.endswith(('.csv', '.tsv', '.csv.gz', '.tsv.gz')):
                if dataset_name in file_name:
                    return os.path.join(self.data_directory, file_name)
        return None