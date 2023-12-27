# tests/test_download_cell_line.py

from coderdata.download.downloader import download_data_by_prefix
from coderdata.load.loader import DatasetLoader
import os
import glob
import pandas as pd

def test_download_data_cell_line():

    #cell_line
    download_data_by_prefix('cell_line')
    
    cell_line_samples = glob.glob('cell_line_samples*')
    assert len(cell_line_samples) > 0, "File cell_line_samples does not exist."
    
    cell_line_experiments = glob.glob('cell_line_experiments*')
    assert len(cell_line_experiments) > 0, "File cell_line_experiments does not exist."
    
    cell_line_mutations = glob.glob('cell_line_mutations*')
    assert len(cell_line_mutations) > 0, "File cell_line_mutations does not exist."
    
    cell_line_proteomics = glob.glob('cell_line_proteomics*')
    assert len(cell_line_proteomics) > 0, "File cell_line_proteomics does not exist."
    
    cell_line_transcriptomics = glob.glob('cell_line_transcriptomics*')
    assert len(cell_line_transcriptomics) > 0, "File cell_line_transcriptomics does not exist."
    
    cell_line_drugs = glob.glob('cell_line_drugs*')
    assert len(cell_line_drugs) > 0, "File cell_line_drugs does not exist."
    
    cell_line_copy_number = glob.glob('cell_line_copy_number*')
    assert len(cell_line_copy_number) > 0, "File cell_line_copy_number does not exist."
    
    
    dataset_type = 'cell_line'
    expected_data_types = ['mutations', 'samples', 'experiments', 'transcriptomics','proteomics','drugs','copy_number']
    
    # Initialize DatasetLoader with the temporary directory
    loader = DatasetLoader(dataset_type)

    # Check if the correct datasets are loaded
    for data_type in expected_data_types:
        assert hasattr(loader, data_type), f"{data_type} dataset not loaded for {dataset_type}"
        loaded_data = getattr(loader, data_type)
        assert isinstance(loaded_data, pd.DataFrame), f"{data_type} is not a DataFrame"
        assert not loaded_data.empty, f"{data_type} DataFrame is empty"