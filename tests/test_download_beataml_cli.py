# tests/test_download_beataml_cli.py

import subprocess
import os
import glob

def test_cli_download():

    # Run the command line tool
    subprocess.run(['coderdata', 'download', '--prefix', 'beataml'], check=True)

    # Check if the expected files are downloaded
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
