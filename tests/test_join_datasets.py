# tests/test_join_datasets.py

import pandas as pd
from coderdata.load.loader import DatasetLoader, join_datasets

def test_join_datasets():
    # Create mock data for two loaders
    data1 = pd.DataFrame({'improve_sample_id': [1, 2], 'entrez_id': [3, 4], 'transcriptomics': [.1, .2]})
    data2 = pd.DataFrame({'improve_sample_id': [1, 2], 'entrez_id': [5, 6], 'transcriptomics': [.3, .4]})

    # Create two DatasetLoader instances
    loader1 = DatasetLoader("dataset1")
    loader2 = DatasetLoader("dataset2")

    # Manually set data for loaders
    loader1.transcriptomics = data1
    loader2.transcriptomics = data2

    # Perform the join
    joined_loader = join_datasets(loader1, loader2)

    # Expected result after join
    expected_data = pd.DataFrame({
        'improve_sample_id': [1, 2, 1, 2],
        'entrez_id': [3, 4, 5, 6],
        'transcriptomics': [.1, .2, .3, .4]
    })

    # Assert that joined data is as expected
    pd.testing.assert_frame_equal(joined_loader.transcriptomics, expected_data)
