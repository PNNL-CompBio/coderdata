# tests/test_download_hcmi.py

from coderdata.download.downloader import download_data_by_prefix
from coderdata.load.loader import DatasetLoader
import os
import glob
import pandas as pd

def test_download_data_hcmi():

    #HCMI
    download_data_by_prefix('hcmi')
    
    hcmi_mutations = glob.glob('hcmi_mutations*')
    assert len(hcmi_mutations) > 0, "File hcmi_mutations does not exist."
    
    hcmi_samples = glob.glob('hcmi_samples*')
    assert len(hcmi_samples) > 0, "File hcmi_samples does not exist."
    
    hcmi_transcriptomics = glob.glob('hcmi_transcriptomics*')
    assert len(hcmi_transcriptomics) > 0, "File hcmi_transcriptomics does not exist."
    
    hcmi_copynum = glob.glob('hcmi_copy_number*')
    assert len(hcmi_copynum) > 0, "File hcmi_copy_number  does not exist."
    
    dataset_type = 'hcmi'
    expected_data_types = ['mutations', 'samples', 'copy_number', 'transcriptomics']
    
    # Initialize DatasetLoader with the temporary directory
    loader = DatasetLoader(dataset_type)

    # Check if the correct datasets are loaded
    for data_type in expected_data_types:
        assert hasattr(loader, data_type), f"{data_type} dataset not loaded for {dataset_type}"
        loaded_data = getattr(loader, data_type)
        assert isinstance(loaded_data, pd.DataFrame), f"{data_type} is not a DataFrame"
        assert not loaded_data.empty, f"{data_type} DataFrame is empty"