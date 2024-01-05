# tests/test_download_cell_line_cli.py

import subprocess
import os
import glob

def test_cli_download():

    # Run the command line tool
    subprocess.run(['coderdata', 'download', '--prefix', 'cell_line'], check=True)

    # Check if the expected files are downloaded
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
    