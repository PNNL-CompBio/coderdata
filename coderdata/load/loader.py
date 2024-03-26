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
            The type of dataset to load (e.g., 'hcmi', 'beataml', 'cptac', 'depmap', "mpnst).
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
        self.genes = pd.DataFrame()
        # self.perturbations = pd.DataFrame()
        self.full = pd.DataFrame()
        self.dataset_type = dataset_type
        self.data_directory = data_directory
        self.load_datasets(dataset_type, data_directory)
        self.data_format_params = {
        'samples': ('improve_sample_id', 'cancer_type', 'model_type',"common_name","other_id", "other_names", "id_source", "species"),
        'transcriptomics': ('improve_sample_id', 'entrez_id', 'transcriptomics'),
        'proteomics': ('improve_sample_id', 'entrez_id', 'proteomics'),
        'mutations': ('improve_sample_id', 'entrez_id', 'mutation'),
        'copy_number': ('improve_sample_id', 'entrez_id', 'copy_number'),
        'methylation': ('improve_sample_id', 'entrez_id', 'methylation'),
        'experiments': ('improve_sample_id', 'improve_drug_id', 'dose_response_value'),
        'drugs': ('improve_drug_id', 'chem_name', 'isoSMILES'),
        'genes': ('entrez_id', 'gene_symbol', 'other_id')
    }

    def load_datasets(self, dataset_type, data_directory):
        print("Processing Data...")
        # Load genes.csv or genes.csv.gz if present
        genes_file_path = None
        for file_name in os.listdir(data_directory):
            if file_name == 'genes.csv' or file_name == 'genes.csv.gz':
                genes_file_path = os.path.join(data_directory, file_name)
                self.genes = self.load_file(genes_file_path)
                print("Loaded genes dataset.")

            if file_name.startswith(dataset_type) and (file_name.endswith(('.csv', '.tsv', '.csv.gz', '.tsv.gz'))):
                file_path = os.path.join(data_directory, file_name)
                dataset_name = file_name[len(dataset_type):].split('.')[0].strip('_')
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
        
    def info(self):
        """
        Display information about the datasets and datatypes.
        """
        # Dataset descriptions
        descriptions = {
            'depmap': 'The cell line datasets were collected from numerous resources such as the LINCS project, DepMap, and the Sanger Institute.',
            'cptac': 'The Clinical Proteomic Tumor Analysis Consortium (CPTAC) project is a collaborative network funded by the National Cancer Institute (NCI).',
            'hcmi': 'Human Cancer Models Initiative (HCMI) data was collected though the National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Portal.',
            'beataml': 'Beat acute myeloid leukemia (BeatAML) data was collected though GitHub and Synapse.',
            'mpnst': 'A collection of NF1-MPNST patient-derived xenografts, organoids, and tumors. Data hosted on synapse.'
        }

        # Check if this is a joined dataset
        if hasattr(self, 'source_of_datatype'):
            print("This is a joined dataset comprising of:")
            for source in set(sum(self.source_of_datatype.values(), [])):
                if source in descriptions:
                    print(f"- {source}: {descriptions[source]}")
        else:
            print(f"Dataset Type: {self.dataset_type}")
            if self.dataset_type in descriptions:
                print(descriptions[self.dataset_type])
        
        # Print information about available datatypes and their formats
        print("\nAvailable Datatypes and Their Formats:")
        for attr in dir(self):
            if isinstance(getattr(self, attr), pd.DataFrame) and not getattr(self, attr).empty:
                format_params = self.data_format_params.get(attr, [])
                if len(format_params) == 3:  # Assuming format_params has three elements: id_vars, var_name, value_name
                    id_vars, var_name, value_name = format_params
                    format_type = 'long' if self.is_long_format(getattr(self, attr), id_vars, var_name, value_name) else 'wide'
                    print(f"- {attr}: {format_type} format")
                else:
                    print(f"- {attr}: Format not specified")

def join_datasets(*args):
    """
    Joins datasets from multiple DatasetLoader instances or dataset type strings. 
    Includes datasets even if they exist in only one instance.

    Parameters
    ----------
    *args : variable number of DatasetLoader instances or dataset type strings
        
    Returns
    -------
    DatasetLoader object
        A new DatasetLoader instance with the joined datasets.
    """
    if not args:
        raise ValueError("At least one dataset or loader must be provided.")

    loaders = []

    for arg in args:
        if isinstance(arg, str):
            loaders.append(DatasetLoader(arg))
        elif isinstance(arg, DatasetLoader):
            loaders.append(arg)
        else:
            raise ValueError("Arguments must be either dataset type strings or DatasetLoader instances.")

    joined_loader = DatasetLoader(loaders[0].dataset_type)
    joined_loader.source_of_datatype = {}

    all_data_types = set()
    for loader in loaders:
        all_data_types.update(set(loader.data_format_params.keys()))

    for attr in all_data_types:
        merged_df = pd.DataFrame()
        sources = set()

        for loader in loaders:
            df = getattr(loader, attr, pd.DataFrame())
            if not df.empty:
                format_params = loader.data_format_params.get(attr, [])
                if len(format_params) >= 3:  # Assuming format_params has three elements: id_vars, var_name, value_name
                    id_vars, var_name, value_name = format_params[:3]
                    if not loader.is_long_format(df, id_vars, var_name, value_name):
                        print(f"Cannot merge {attr} as it is not in long format in {loader.dataset_type}.")
                        return None  

                merge_keys = loader.data_format_params.get(attr, [])

                for key in merge_keys:
                    if key in df.columns:
                        if merged_df.empty or key not in merged_df.columns:
                            continue

                        # Check if data types are different
                        if df[key].dtype != merged_df[key].dtype:
                            try:
                                # Try converting both columns to float
                                df[key] = df[key].astype(float)
                                merged_df[key] = merged_df[key].astype(float)
                            except ValueError:
                                # If conversion to float fails, convert to string
                                df[key] = df[key].astype(str)
                                merged_df[key] = merged_df[key].astype(str)

                if merged_df.empty:
                    merged_df = df
                else:
                    merged_df = pd.concat((merged_df, df))
                    
                if hasattr(loader, 'source_of_datatype') and attr in loader.source_of_datatype:
                    sources.update(loader.source_of_datatype[attr])
                else:
                    sources.add(loader.dataset_type)

        if not merged_df.empty:
            setattr(joined_loader, attr, merged_df)
            joined_loader.source_of_datatype[attr] = list(sources)

    return joined_loader
