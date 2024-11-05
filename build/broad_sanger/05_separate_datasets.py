import pandas as pd
import numpy as np
import os
import coderdata as cd
from coderdata import DatasetLoader
import psutil

# Function to monitor memory usage
def print_memory_usage(message):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2)  # Convert bytes to MB
    print(f"{message}: {mem:.2f} MB")

# Function to split the experiments by study
def split_experiments_by_study(data: DatasetLoader) -> dict:
    """
    Splits the DatasetLoader object into multiple smaller DatasetLoader objects
    according to the `study` recorded in the `.experiments` table in the DatasetLoader object.

    Parameters
    ----------
    data : DatasetLoader
        The DatasetLoader object containing the data set loaded into memory
        via `coderdata.DatasetLoader()`.

    Returns
    -------
    dict
        A dictionary dict[study, data] where keys `study` are the names 
        of the study in the `.experiments` part of the imported 
        DatasetLoader object and values `data` are the filtered smaller
        DatasetLoader objects containing only data corresponding to the 
        study.
    """

    df_ret = {}
    experiments = data.experiments
    # creating the groups based on 'study' to iterate over
    groups = experiments.groupby('study')
    for name, group in groups:

        # extracting improve sample and drug ids from the provided split
        sample_ids = list(np.unique(group['improve_sample_id'].values))
        drug_ids = list(np.unique(group['improve_drug_id'].values))

        # creating new DatasetLoader objects that contain only data
        # pertaining to the study defined by the previous grouping
        df_ret[name] = _filter(
            data=data, sample_ids=sample_ids, drug_ids=drug_ids, study=name
        )

    return df_ret

# Updated _filter function without deepcopy
def _filter(
        data: DatasetLoader,
        sample_ids: list,
        drug_ids: list,
        study: str = None,
        ) -> DatasetLoader:
    """
    Helper function to filter down the DatasetLoader object(s) to create
    independent, more concise DatasetLoader objects for further processing.
    This can be either splitting a dataset according to the different 
    drug response studies (e.g., the broad_sanger dataset) or if small 
    subsets need to be extracted (e.g., training/testing splits for 
    machine learning).

    Parameters
    ----------
    data : DatasetLoader
        Contains a full DatasetLoader object imported/loaded via 
        `cd.DatasetLoader`
    sample_ids : list
        A list of improve_sample_id[s] that the DatasetLoader object should
        be filtered to
    drug_ids : list
        A list of improve_drug_id[s] that the DatasetLoader object should 
        be filtered to
    study : str, default = None
        The drug response study that the DatasetLoader object should be 
        filtered to. This argument is only important for filtering the
        broad_sanger dataset if the splitting/filtering of the data 
        set is based on the drug response study

    Returns
    -------
    DatasetLoader
        The filtered DatasetLoader object
    """

    # Create a new DatasetLoader instance
    data_ret = DatasetLoader()

    # Filtering each individual data type down by only the improve 
    # sample/drug ids that are present in the study
    if not data.copy_number.empty:
        data_ret.copy_number = data.copy_number[
            data.copy_number['improve_sample_id'].isin(sample_ids)
        ]
    if not data.drugs.empty:
        data_ret.drugs = data.drugs[
            data.drugs['improve_drug_id'].isin(drug_ids)
        ]
    if not data.mutations.empty:
        data_ret.mutations = data.mutations[
            data.mutations['improve_sample_id'].isin(sample_ids)
        ]
    if not data.proteomics.empty:
        data_ret.proteomics = data.proteomics[
            data.proteomics['improve_sample_id'].isin(sample_ids)
        ]
    if not data.samples.empty:
        data_ret.samples = data.samples[
            data.samples['improve_sample_id'].isin(sample_ids)
        ]
    if not data.transcriptomics.empty:
        data_ret.transcriptomics = data.transcriptomics[
            data.transcriptomics['improve_sample_id'].isin(sample_ids)
        ]
    if not data.experiments.empty:
        data_ret.experiments = data.experiments[
            data.experiments['study'] == study
        ]

    return data_ret

# Load the main dataset
bs = cd.DatasetLoader("broad_sanger")
print_memory_usage("After loading 'broad_sanger'")

# Split the experiments by study
split = split_experiments_by_study(bs)
print_memory_usage("After splitting experiments by study")

# Delete 'bs' to free memory
del bs

datasets_to_process = ["CCLE", "CTRPv2", "PRISM", "GDSCv1", "GDSCv2", "FIMM", "gCSI", "NCI60"]

dataset_sources = {
    "CCLE": ["Broad"],
    "CTRPv2": ["Broad"],
    "PRISM": ["Broad"],
    "GDSCv1": ["Sanger"],
    "GDSCv2": ["Sanger"],
    "FIMM": ["Broad"],
    "gCSI": ["Broad"],  # gCSI generates its own omics data but it is comparable to CCLE. In future, retrieve gCSI omics.
    "NCI60": ["Broad"]
}

omics_datatypes = ["copy_number", "mutations", "proteomics", "samples", "transcriptomics"]
tsv_data_types = ["drugs", "experiments"]
output_datatypes = ["copy_number", "mutations", "proteomics", "samples", "transcriptomics", "drugs", "experiments"]
output_dir = "/tmp"

for dataset_name in datasets_to_process:
    print(f"Processing dataset: {dataset_name}")
    dataset = split[dataset_name]
    print_memory_usage(f"After loading dataset {dataset_name}")
    sources = dataset_sources[dataset_name]

    # Filter each datatype in the dataset
    for datatype in omics_datatypes:
        data = getattr(dataset, datatype)
        if data.empty:
            continue

        if 'source' in data.columns:
            filtered_data = data[data['source'].isin(sources)].copy()
            setattr(dataset, datatype, filtered_data)
        elif 'study' in data.columns:
            filtered_data = data[data['study'].isin(sources)].copy()
            setattr(dataset, datatype, filtered_data)
        else:
            continue

        print_memory_usage(f"After filtering {datatype} in {dataset_name}")

    # Save and delete data to free memory
    for data_type in output_datatypes:
        data = getattr(dataset, data_type)
        if data.empty:
            continue

        if data_type in tsv_data_types:
            file_extension = ".tsv"
            sep = '\t'
        else:
            file_extension = ".csv"
            sep = ','

        filename = f"{dataset_name}_{data_type}{file_extension}".lower()
        filepath = os.path.join(output_dir, filename)

        data.to_csv(filepath, sep=sep, index=False)
        print(f"Saved {data_type} from {dataset_name} to {filepath}")
        print_memory_usage(f"After saving {data_type} in {dataset_name}")

        # Delete data to free memory
        del data
        setattr(dataset, data_type, None)
        print_memory_usage(f"After deleting {data_type} in {dataset_name}")

    # Delete the dataset to free memory
    del dataset
    del split[dataset_name]
    print_memory_usage(f"After deleting dataset {dataset_name}")

print("Processing complete.")
