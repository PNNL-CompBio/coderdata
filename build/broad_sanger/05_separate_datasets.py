import pandas as pd
import numpy as np
import os
from copy import deepcopy
import coderdata as cd
from coderdata import DatasetLoader

#from Yannick stats.py file in utils.
def split_experiments_by_study(data: DatasetLoader) -> dict:
    """
    Splits the CoderData object into multiple smaller CoderData objects
    according to the `study` recorded in the ``.experiments`` table in 
    the CoderData object.

    Parameters
    ----------
    data : DatasetLoader
        The CoderData object containing the data set loaded into memory
        via ``coderdata.DatasetLoader()``.

    Returns
    -------
    dict
        A dictionary dict[study, data] where keys `study` are the names 
        of the study in the ``.experiments`` part of the imported 
        CoderData object and values `data` are the filtered smaller
        CoderData objects containing only data corresponding to the 
        study. 
    """

    df_ret = {}
    experiments = data.experiments
    print(experiments)
    print(experiments.source)
    print(experiments.study)
    # creating the groups based on 'study' to itterate over 
    groups = experiments.groupby('study')
    for name, group in groups:

        # extracting improve sample and drug ids from the provided split
        sample_ids = list(np.unique(group['improve_sample_id'].values))
        drug_ids = list(np.unique(group['improve_drug_id'].values))
        
        # creating new CoderData objects that contain only data
        # pertaining to the study defined by the previous grouping
        df_ret[name] = _filter(
            data=data, sample_ids=sample_ids, drug_ids=drug_ids, study=name
            )
    
    return df_ret


def _filter(
        data: DatasetLoader,
        sample_ids: list,
        drug_ids: list,
        study: str=None,
        ) -> DatasetLoader:
    """
    Helper function to filter down the CoderData object(s) to create
    independent more concise CoderData objects for further processing.
    This can be either splitting a dataset according to the different 
    drug response studies (e.g. the broad_sanger dataset) or if small 
    subsets need to be extracted (e.g. training / testing splits for 
    machine learning)

    Parameters
    ----------
    data : DatasetLoader
        Contains a full CoderData object imported/loaded via 
        ``cd.DataLoader``
    sample_ids : list
        A list of improve_sample_id[s] that the CoderData object should
        be filtered to
    drug_ids : list
        A list of improve_drug_id[s] that the CoderData object should 
        be filtered to
    study : str, default = None
        The drug response study that the CoderData object should be 
        filtered to. This argument is only important for filtering the
        broad_sanger dataset if the splitting / filtering of the data 
        set is based on the drug response study

    Returns
    -------
    DatasetLoader
        The filtered CoderData object
    
    Notes
    -----

    Different data types of the CoderData object are going to be 
    filtered using either the improve_sample_id or the improve_drug_id.
    
    - cd.copynumber -> reduce based on ``improve_sample_id``
    - cd.drugs -> reduce based on ``improve_drug_id``
    - cd.experiments -> reduce based on ``study`` (only applicable if 
      the dataset is broad_sanger)
    - cd.mutations -> reduce based on ``improve_sample_id``
    - cd.proteomics -> reduce based on ``improve_sample_id``
    - cd.samples -> reduce based on ``improve_sample_id``
    - cd.transcriptomics -> reduce based on ``improve_sample_id``
    
    """

    # creating a deep copy of the CoderData object such that any 
    # further operations on the object are not changing the original
    # object / data
    data_ret = deepcopy(data)

    # filtering each individual data type down by only the improve 
    # sample / drug ids that are present in the study
    if not data_ret.copy_number.empty:
        data_ret.copy_number = data_ret.copy_number[
            data_ret.copy_number['improve_sample_id'].isin(sample_ids)
        ]
    if not data_ret.drugs.empty:
        data_ret.drugs = data_ret.drugs[
            data_ret.drugs['improve_drug_id'].isin(drug_ids)
            ]
    if not data_ret.mutations.empty:
        data_ret.mutations = data_ret.mutations[
            data_ret.mutations['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.proteomics.empty:
        data_ret.proteomics = data_ret.proteomics[
            data_ret.proteomics['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.samples.empty:
        data_ret.samples = data_ret.samples[
            data_ret.samples['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.transcriptomics.empty:
        data_ret.transcriptomics = data_ret.transcriptomics[
            data_ret.transcriptomics['improve_sample_id'].isin(sample_ids)
            ]
    if not data_ret.experiments.empty:
        data_ret.experiments = data_ret.experiments[
            data_ret.experiments['study'] == study
        ]
    
    return data_ret


bs = cd.DatasetLoader("broad_sanger")
split = split_experiments_by_study(bs)


datasets = {
    "CCLE": split["CCLE"],
    "CTRPv2": split["CTRPv2"],
    "FIMM": split["FIMM"],
    "GDSCv1": split["GDSCv1"],
    "GDSCv2": split["GDSCv2"],
    "NCI60": split["NCI60"],
    "PRISM": split["PRISM"],
    "gCSI": split["gCSI"]
}

dataset_sources = {
    "CCLE": ["Broad"],
    "CTRPv2": ["Broad"],
    "PRISM": ["Broad"],
    "GDSCv1": ["Sanger"],
    "GDSCv2": ["Sanger"],
    "FIMM": ["Broad"],
    "gCSI": ["Broad"],  # gCSI generates its own omics data but it is comparable to CCLE. In future, retrive gCSI omics.
    "NCI60": ["Broad"]
}

# Define the datasets and omics_types to filter
datasets_to_process = ["CCLE", "CTRPv2", "PRISM", "GDSCv1", "GDSCv2", "FIMM", "gCSI", "NCI60"]
omics_datatypes = ["copy_number","mutations","proteomics","samples","transcriptomics"]


for dataset_name in datasets_to_process:
    dataset = datasets[dataset_name]
    sources = dataset_sources[dataset_name]
    

    # Filter each datatype in the dataset
    for datatype in omics_datatypes:
        data = getattr(dataset, datatype)
        
        if 'source' in data.columns:
            filtered_data = data[data['source'].isin(sources)].copy()
            setattr(dataset, datatype, filtered_data)
        elif 'study' in data.columns:
            # Some data may have 'study' instead of 'source'
            filtered_data = data[data['study'].isin(sources)].copy()
            setattr(dataset, datatype, filtered_data)
        else:
            # If neither 'source' nor 'study' is present, keep the data as is.
            continue
        
        
output_dir = "/tmp"
tsv_data_types = ["drugs", "experiments"]

output_datatypes = ["copy_number","mutations","proteomics","samples","transcriptomics","drugs","experiments"]
for dataset_name in datasets_to_process:
    dataset = datasets[dataset_name]
    
    for data_type in output_datatypes:
        data = getattr(dataset, data_type)
        
        # Determine the file extension and separator
        if data_type in tsv_data_types:
            file_extension = ".tsv"
            sep = '\t'
        else:
            file_extension = ".csv"
            sep = ','
        
        # Construct the filename
        # Convert to lowercase
        filename = f"{dataset_name}_{data_type}{file_extension}".lower()
        filepath = os.path.join(output_dir, filename)
        # Save the DataFrame to file
        data.to_csv(filepath, sep=sep, index=False)
        print(f"Saved {data_type} from {dataset_name} to {filepath}")
