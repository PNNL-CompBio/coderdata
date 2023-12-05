# tests/test_download_beataml.py

from coderdata.download.downloader import download_data_by_prefix
from coderdata.load.loader import DatasetLoader
import os
import glob
import pandas as pd

def test_download_data_beataml():

    #BeatAML
    download_data_by_prefix('beataml')
    
    beataml_drugs = glob.glob('beataml_drugs*')
    assert len(beataml_drugs) > 0, "File beataml_drugs  does not exist."
    
    beataml_experiments = glob.glob('beataml_experiments*')
    assert len(beataml_experiments) > 0, "File beataml_experiments does not exist."
    
    beataml_mutations = glob.glob('beataml_mutations*')
    assert len(beataml_mutations) > 0, "File beataml_mutations does not exist."
    
    beataml_proteomics = glob.glob('beataml_proteomics*')
    assert len(beataml_proteomics) > 0, "File beataml_proteomics  does not exist."
    
    beataml_samples = glob.glob('beataml_samples*')
    assert len(beataml_samples) > 0, "File beataml_samples does not exist."
    
    beataml_transcriptomics = glob.glob('beataml_transcriptomics*')
    assert len(beataml_transcriptomics) > 0, "File beataml_transcriptomics does not exist."


    dataset_type = 'beataml'
    expected_data_types = ['mutations', 'samples', 'experiments', 'transcriptomics','proteomics','drugs']
    
    loader = DatasetLoader(dataset_type)

    # Check if the correct datasets are loaded
    for data_type in expected_data_types:
        assert hasattr(loader, data_type), f"{data_type} dataset not loaded for {dataset_type}"
        loaded_data = getattr(loader, data_type)
        assert isinstance(loaded_data, pd.DataFrame), f"{data_type} is not a DataFrame"
        assert not loaded_data.empty, f"{data_type} DataFrame is empty"