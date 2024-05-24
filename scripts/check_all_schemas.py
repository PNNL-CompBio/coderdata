import subprocess
import os
import argparse
import concurrent.futures

def run_schema_checker(script_name):
    """
    Runs a schema checker shell script and returns the outcome.
    
    Parameters:
        script_name (str): The filename of the shell script to run.
    
    Returns:
        tuple: (script_name, bool) where bool is True if the validation succeeded, False if it failed.
    """
    try:
        # Build the full path to the script
        script_path = os.path.join('schema', script_name)
        # Run the shell script
        result = subprocess.run(['bash', script_path], capture_output=True, text=True)
        # Print output and error from shell script
        print(result.stdout)
        if result.stderr:
            print(f"Error in {script_name}:", result.stderr)
        # Return the script name and True if the script executed successfully
        return (script_name, result.returncode == 0)
    except Exception as e:
        print(f"An error occurred while running {script_name}: {e}")
        return (script_name, False)

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

    scripts_to_run = schema_mapping.values() if not args.datasets else [schema_mapping[dataset] for dataset in args.datasets if dataset in schema_mapping]

    all_passed = True
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {executor.submit(run_schema_checker, script): script for script in scripts_to_run}

        for future in concurrent.futures.as_completed(futures):
            script_name, result = future.result()
            if not result:
                all_passed = False
                print(f"Validation failed for {script_name}")

    if all_passed:
        print("All schema validations passed successfully.")
    else:
        print("Some schema validations failed.")

if __name__ == '__main__':
    main()
