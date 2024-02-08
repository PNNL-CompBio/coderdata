# tests/test_reload_all.py


import pandas as pd
import os
from coderdata.load.loader import DatasetLoader
from coderdata.download import download_data_by_prefix


def test_reload_all_datasets():
    download_data_by_prefix("hcmi")
    loader = DatasetLoader("hcmi")
    
    # Ensure the datasets are initially empty
    loader.transcriptomics = pd.DataFrame()
    loader.copy_number = pd.DataFrame()
    loader.samples = pd.DataFrame()
    loader.mutations = pd.DataFrame()

    # Reload all datasets
    loader.reload_datasets()
    
    # Verify that each dataset is reloaded correctly
    reloaded_transcriptomics = pd.read_csv('hcmi_transcriptomics.csv.gz')
    reloaded_copy_number = pd.read_csv('hcmi_copy_number.csv.gz')
    reloaded_samples = pd.read_csv('hcmi_samples.csv')
    reloaded_mutations = pd.read_csv('hcmi_mutations.csv.gz')

    pd.testing.assert_frame_equal(loader.transcriptomics, reloaded_transcriptomics)
    pd.testing.assert_frame_equal(loader.copy_number, reloaded_copy_number)
    pd.testing.assert_frame_equal(loader.samples, reloaded_samples)
    pd.testing.assert_frame_equal(loader.mutations, reloaded_mutations)