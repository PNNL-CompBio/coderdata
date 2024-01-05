# tests/test_download_cptac_cli.py

import subprocess
import tempfile
import os
import glob

def test_cli_download():

    # Run the command line tool
    subprocess.run(['coderdata', 'download', '--prefix', 'cptac'], check=True)

    # Check if the expected files are downloaded
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
    