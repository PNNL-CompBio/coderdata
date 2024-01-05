# tests/test_download_hcmi_cli.py

import subprocess
import tempfile
import os
import glob

def test_cli_download():
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Change the current working directory to the temporary directory
        os.chdir(temp_dir)

        # Run the command line tool
        subprocess.run(['coderdata', 'download', '--prefix', 'hcmi'], check=True)

        # Check if the expected files are downloaded
        hcmi_mutations = glob.glob('hcmi_mutations*')
        assert len(hcmi_mutations) > 0, "File hcmi_mutations does not exist."
        
        hcmi_samples = glob.glob('hcmi_samples*')
        assert len(hcmi_samples) > 0, "File hcmi_samples does not exist."
        
        hcmi_transcriptomics = glob.glob('hcmi_transcriptomics*')
        assert len(hcmi_transcriptomics) > 0, "File hcmi_transcriptomics does not exist."
        
        hcmi_copynum = glob.glob('hcmi_copy_number*')
        assert len(hcmi_copynum) > 0, "File hcmi_copy_number  does not exist."
