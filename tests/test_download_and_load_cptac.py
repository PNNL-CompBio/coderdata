
# tests/test_download_cptac.py

from coderdata.download.downloader import download_data_by_prefix
from coderdata.load.loader import DatasetLoader
import os
import glob
import pandas as pd

def test_download_data_cptac():

    #CPTAC
    download_data_by_prefix('cptac')
    
    cptac_copy_number = glob.glob('cptac_copy_number*')
    assert len(cptac_copy_number) > 0, "File cptac_copy_number does not exist."
    
    cptac_proteomics = glob.glob('cptac_proteomics*')
    assert len(cptac_proteomics) > 0, "File cptac_proteomics does not exist."
    
    cptac_samples = glob.glob('cptac_samples*')
    assert len(cptac_samples) > 0, "File cptac_samples does not exist."
    
    cptac_mutations = glob.glob('cptac_mutations*')
    assert len(cptac_mutations) > 0, "File cptac_mutations does not exist."
    
    cptac_transcriptomics = glob.glob('cptac_transcriptomics*')
    assert len(cptac_transcriptomics) > 0, "File cptac_transcriptomics does not exist."
    
    dataset_type = 'cptac'
    expected_data_types = ['mutations', 'samples', 'transcriptomics','proteomics','copy_number']
    
    # Initialize DatasetLoader with the temporary directory
    loader = DatasetLoader(dataset_type)

    # Check if the correct datasets are loaded
    for data_type in expected_data_types:
        assert hasattr(loader, data_type), f"{data_type} dataset not loaded for {dataset_type}"
        loaded_data = getattr(loader, data_type)
        assert isinstance(loaded_data, pd.DataFrame), f"{data_type} is not a DataFrame"
        assert not loaded_data.empty, f"{data_type} DataFrame is empty"