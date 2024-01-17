# tests/test_info.py

import pytest
from coderdata.load.loader import DatasetLoader
from io import StringIO
import sys
import pandas as pd

def test_info_function():
    # Create a mock DatasetLoader instance
    loader = DatasetLoader("hcmi")
    loader.transcriptomics = pd.DataFrame({'improve_sample_id': [1, 2], 'entrez_id': [3, 4], 'transcriptomics': [.1, .2]})
    loader.proteomics = pd.DataFrame({'improve_sample_id': [1, 2], 'entrez_id': [5, 6], 'proteomics': [.3, .4]})
    loader.data_format_params = {
        'transcriptomics': ('improve_sample_id', 'entrez_id', 'transcriptomics'),
        'proteomics': ('improve_sample_id', 'entrez_id', 'proteomics')
    }

    # Capture the output of the info function
    captured_output = StringIO()
    sys.stdout = captured_output
    loader.info()
    sys.stdout = sys.__stdout__

    # Define the expected output
    expected_output = (
        "Dataset Type: hcmi\n"
        "Human Cancer Models Initiative (HCMI) data was collected though the National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Portal.\n\n"
        "Available Datatypes and Their Formats:\n"
        "- copy_number: Format not specified\n"
        "- mutations: Format not specified\n"
        "- proteomics: long format\n"
        "- samples: Format not specified\n"
        "- transcriptomics: long format\n"
    )

    # Assert that the captured output matches the expected output
    assert captured_output.getvalue() == expected_output
