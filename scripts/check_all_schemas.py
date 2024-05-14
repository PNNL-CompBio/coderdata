import subprocess
import os
import argparse

def run_schema_checker(script_name):
    """
    Runs a schema checker shell script and returns the outcome.
    
    Parameters:
        script_name (str): The filename of the shell script to run.
    
    Returns:
        bool: True if the validation succeeded, False if it failed.
    """
    try:
        # Build the full path to the script
        script_path = os.path.join('schema', script_name)
        # Run the shell script
        result = subprocess.run(['bash', script_path], capture_output=True, text=True)
        # Print output and error from shell script
        print(result.stdout)
        if result.stderr:
            print("Error:", result.stderr)
        # Return True if the script executed successfully
        return result.returncode == 0
    except Exception as e:
        print(f"An error occurred while running {script_name}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Run schema validations for specified datasets.")
    parser.add_argument('-d', '--datasets', nargs='*', help='List of datasets to validate (e.g., beataml, cptac, depmap, hcmi)', default=None)
    args = parser.parse_args()

    # Mapping from dataset names to script names
    schema_mapping = {
        'beataml': 'check_beataml_linkml.sh',
        'cptac': 'check_cptac_linkml.sh',
        'depmap': 'check_depmap_linkml.sh',
        'hcmi': 'check_hcmi_linkml.sh',
        'mpnst': 'check_mpnst_linkml.sh'
    }

    all_passed = True
    scripts_to_run = schema_mapping.values() if not args.datasets else [schema_mapping[dataset] for dataset in args.datasets if dataset in schema_mapping]

    # Iterate over each script and run it
    for script_name in scripts_to_run:
        print(f"Running {script_name}...")
        if not run_schema_checker(script_name):
            all_passed = False
            print(f"Validation failed for {script_name}")

    if all_passed:
        print("All schema validations passed successfully.")
    else:
        print("Some schema validations failed.")

if __name__ == '__main__':
    main()
