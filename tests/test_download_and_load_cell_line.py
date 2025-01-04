# tests/test_download_depmap.py

from coderdata.download.downloader import download
from coderdata.load.loader import DatasetLoader
import os
import glob
import pandas as pd

def test_download_data_depmap():

    #depmap
    download('depmap')
    
    depmap_samples = glob.glob('depmap_samples*')
    assert len(depmap_samples) > 0, "File depmap_samples does not exist."
    
    depmap_experiments = glob.glob('depmap_experiments*')
    assert len(depmap_experiments) > 0, "File depmap_experiments does not exist."
    
    depmap_mutations = glob.glob('depmap_mutations*')
    assert len(depmap_mutations) > 0, "File depmap_mutations does not exist."
    
    depmap_proteomics = glob.glob('depmap_proteomics*')
    assert len(depmap_proteomics) > 0, "File depmap_proteomics does not exist."
    
    depmap_transcriptomics = glob.glob('depmap_transcriptomics*')
    assert len(depmap_transcriptomics) > 0, "File depmap_transcriptomics does not exist."
    
    depmap_drugs = glob.glob('depmap_drugs*')
    assert len(depmap_drugs) > 0, "File depmap_drugs does not exist."
    
    depmap_copy_number = glob.glob('depmap_copy_number*')
    assert len(depmap_copy_number) > 0, "File depmap_copy_number does not exist."
    
    
    dataset_type = 'depmap'
    expected_data_types = ['mutations', 'samples', 'experiments', 'transcriptomics','proteomics','drugs','copy_number']
    
    # Initialize DatasetLoader with the temporary directory
    loader = DatasetLoader(dataset_type)

    # Check if the correct datasets are loaded
    for data_type in expected_data_types:
        assert hasattr(loader, data_type), f"{data_type} dataset not loaded for {dataset_type}"
        loaded_data = getattr(loader, data_type)
        assert isinstance(loaded_data, pd.DataFrame), f"{data_type} is not a DataFrame"
        assert not loaded_data.empty, f"{data_type} DataFrame is empty"