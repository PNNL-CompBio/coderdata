# tests/test_reformat_datasets.py

from coderdata.load.loader import DatasetLoader
import pandas as pd

def test_reformat_dataset():
    # Initialize DatasetLoader with sample data
    loader = DatasetLoader("test")
    loader.transcriptomics = pd.DataFrame({
        'improve_sample_id': [1, 1, 2, 2],
        'entrez_id': ['gene1', 'gene2', 'gene1', 'gene2'],
        'transcriptomics': [10, 20, 30, 40]
    })

    # Convert to wide format and test
    loader.reformat_dataset('transcriptomics', 'wide')
    wide_data = loader.transcriptomics
    assert wide_data.shape == (2, 3), "Incorrect shape after converting to wide format"
    assert 'gene1' in wide_data.columns and 'gene2' in wide_data.columns, "Columns missing in wide format"

    # Convert back to long format and test
    loader.reformat_dataset('transcriptomics', 'long')
    long_data = loader.transcriptomics
    assert long_data.shape == (4, 3), "Incorrect shape after converting to long format"
    assert 'entrez_id' in long_data.columns, "Column 'entrez_id' missing in long format"

    # Test with a dataset that doesn't exist
    try:
        loader.reformat_dataset('nonexistent_dataset', 'wide')
    except ValueError as e:
        assert str(e) == "Dataset 'nonexistent_dataset' is empty or already in desired format."

    # Test with invalid format argument
    try:
        loader.reformat_dataset('transcriptomics', 'invalid_format')
    except ValueError as e:
        assert str(e) == "Format must be 'long' or 'wide'."

    # Ensure original data is preserved when attempting to reformat already formatted data
    loader.reformat_dataset('transcriptomics', 'wide')
    assert loader.transcriptomics.equals(wide_data), "Data altered when reformatting wide format data the second time"

    loader.reformat_dataset('transcriptomics', 'long')
    assert loader.transcriptomics.equals(long_data), "Data altered when reformatting long format data the second time"
