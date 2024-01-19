# tests/test_info.py

from coderdata.load.loader import DatasetLoader
from io import StringIO
import sys
import re
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

    # Get captured output as a string
    output = captured_output.getvalue()

    # Define expected patterns
    expected_patterns = [
        r"Dataset Type",
        r"Human Cancer Models Initiative",
        r"Available Datatypes",
        r"transcriptomics",
        r"proteomics"
    ]

    # Check for each pattern in output
    for pattern in expected_patterns:
        assert re.search(pattern, output), f"Pattern not found in output: {pattern}"
