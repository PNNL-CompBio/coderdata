import subprocess
import argparse
import concurrent.futures
import sys
import yaml
import os

def run_validations_for_dataset(dataset_name, validations):
    """
    Runs validations for a dataset and returns whether all validations passed.

    Parameters:
        dataset_name (str): The name of the dataset.
        validations (list): List of validations, each a dict with 'target_class' and 'file'.

    Returns:
        tuple: (dataset_name, bool) where bool is True if all validations passed, False otherwise.
    """
    validation_failed = False
    for validation in validations:
        target_class = validation['target_class']
        file_path = validation['file']
        print(f"Validating {target_class} in file {file_path} for dataset {dataset_name}...")

        # Run the validation command
        try:
            result = subprocess.run(
                ['linkml-validate', '--schema', 'schema/coderdata.yaml', '--target-class', target_class, file_path],
                stdout=sys.stdout, stderr=sys.stderr, text=True
            )
            
            # Print output and error from command
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(f"Error in validating {target_class} in {file_path}:", result.stderr)
            # Check the exit status
            if result.returncode != 0:
                print(f"Validation failed for {target_class} in file {file_path}.")
                validation_failed = True
            else:
                print(f"Validation succeeded for {target_class} in file {file_path}.")
        except Exception as e:
            print(f"An error occurred while validating {target_class} in {file_path}: {e}")
            validation_failed = True

    if validation_failed:
        print(f"One or more validations failed for dataset {dataset_name}.")
        return (dataset_name, False)
    else:
        print(f"All validations succeeded for dataset {dataset_name}.")
        return (dataset_name, True)

def main():
    parser = argparse.ArgumentParser(description="Run schema validations for specified datasets.")
    parser.add_argument('-d', '--datasets', nargs='*', help='List of datasets to validate (e.g., "beataml cptac ccle hcmi")', default=None)
    args = parser.parse_args()

    # Read the config file
    config_path = os.path.join('schema', 'expected_files.yaml')
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    available_datasets = config['datasets'].keys()
    datasets_to_validate = args.datasets if args.datasets else available_datasets
    datasets_to_validate = [dataset for dataset in datasets_to_validate if dataset in available_datasets]

    print(f"Datasets to validate: {datasets_to_validate}")

    all_passed = True
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(
                run_validations_for_dataset,
                dataset,
                config['datasets'][dataset]
            ): dataset for dataset in datasets_to_validate
        }

        for future in concurrent.futures.as_completed(futures):
            dataset_name, result = future.result()
            if not result:
                all_passed = False
                print(f"Validation failed for dataset {dataset_name}")
            else:
                print(f"Validation passed for dataset {dataset_name}")

    if all_passed:
        print("All schema validations passed successfully.")
        sys.exit(0)
    else:
        print("Some schema validations failed.")
        sys.exit(1)

if __name__ == '__main__':
    main()
