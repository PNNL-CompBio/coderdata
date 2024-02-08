# tests/test_reload_one.py


import pandas as pd
import os
from coderdata.load.loader import DatasetLoader
from coderdata.download import download_data_by_prefix


def test_reload_specific_dataset():
    download_data_by_prefix("hcmi")
    loader = DatasetLoader("hcmi")
    
    # Ensure the dataset is initially empty
    loader.transcriptomics = pd.DataFrame()

    # Reload a specific dataset
    loader.reload_datasets('transcriptomics')
    # Load in the actual data again to compare
    reloaded_data = pd.read_csv('hcmi_transcriptomics.csv.gz')
    pd.testing.assert_frame_equal(loader.transcriptomics, reloaded_data)
