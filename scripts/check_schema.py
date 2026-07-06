import argparse
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import yaml

# Path to the LinkML schema, relative to the repo root (the working directory
# used both locally and inside the `upload` Docker image).
SCHEMA_PATH = os.path.join('schema', 'coderdata.yaml')

# One Validator per worker process. Constructing a Validator compiles the
# schema (and, lazily, the per-target-class JSON Schema) exactly once and
# caches it, so reusing a single instance across many files avoids the
# per-file process startup + schema recompilation that dominated runtime when
# this script shelled out to `linkml-validate` once per file.
#
# The plugin configuration below reproduces the exact defaults of the
# `linkml-validate` CLI (JsonschemaValidationPlugin with closed=True, and
# strict=False), so validation is identical to the previous implementation --
# same engine, same closed-world JSON Schema, same results.
_VALIDATOR = None


def _init_worker():
    """Initializer run once per worker process: build the shared Validator."""
    global _VALIDATOR
    # LinkML emits chatty INFO logs while importing the schema; silence them so
    # the validation output stays readable.
    logging.disable(logging.INFO)
    from linkml.validator import Validator
    from linkml.validator.plugins import JsonschemaValidationPlugin

    _VALIDATOR = Validator(
        SCHEMA_PATH,
        validation_plugins=[JsonschemaValidationPlugin(closed=True)],
    )


def _validate_one(task):
    """Validate a single file against its target class using the shared Validator.

    Returns (dataset_name, target_class, file_path, [error_messages]).
    """
    dataset_name, target_class, file_path = task
    from linkml.validator.loaders import default_loader_for_file

    loader = default_loader_for_file(file_path)
    errors = []
    for result in _VALIDATOR.iter_results_from_source(loader, target_class):
        severity = getattr(result.severity, 'value', str(result.severity))
        message = (
            f"[{severity}] [{loader.source}/{result.instance_index}] "
            f"{result.message}"
        )
        # linkml-validate exits non-zero only on ERROR (or worse); mirror that
        # so a WARNING/INFO does not fail the build.
        if str(severity).upper() in ('ERROR', 'FATAL', 'CRITICAL'):
            errors.append(message)
    return (dataset_name, target_class, file_path, errors)


def main():
    parser = argparse.ArgumentParser(
        description="Run schema validations for specified datasets."
    )
    parser.add_argument(
        '-d', '--datasets', nargs='*',
        help='List of datasets to validate (e.g., "beataml cptac ccle hcmi")',
        default=None,
    )
    args = parser.parse_args()

    config_path = os.path.join('schema', 'expected_files.yaml')
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    available_datasets = list(config['datasets'].keys())
    datasets_to_validate = args.datasets if args.datasets else available_datasets
    datasets_to_validate = [d for d in datasets_to_validate if d in available_datasets]

    print(f"Datasets to validate: {datasets_to_validate}")

    # Flatten to per-file tasks so all files validate in parallel across
    # processes (the previous version only parallelized across datasets).
    tasks = []
    for dataset in datasets_to_validate:
        for validation in config['datasets'][dataset]:
            tasks.append((dataset, validation['target_class'], validation['file']))

    # Track pass/fail per dataset so the summary matches the old output.
    dataset_failed = {d: False for d in datasets_to_validate}

    # Cap workers: each worker can load a multi-GB decompressed omics file into
    # memory, so bound concurrency to avoid exhausting RAM on large builds.
    max_workers = min(len(tasks), os.cpu_count() or 1, 8) or 1
    with ProcessPoolExecutor(max_workers=max_workers, initializer=_init_worker) as executor:
        for dataset_name, target_class, file_path, errors in executor.map(_validate_one, tasks):
            print(f"Validating {target_class} in file {file_path} for dataset {dataset_name}...")
            if errors:
                for message in errors:
                    print(message)
                print(f"Validation failed for {target_class} in file {file_path}.")
                dataset_failed[dataset_name] = True
            else:
                print(f"Validation succeeded for {target_class} in file {file_path}.")

    all_passed = True
    for dataset in datasets_to_validate:
        if dataset_failed[dataset]:
            all_passed = False
            print(f"One or more validations failed for dataset {dataset}.")
            print(f"Validation failed for dataset {dataset}")
        else:
            print(f"All validations succeeded for dataset {dataset}.")
            print(f"Validation passed for dataset {dataset}")

    if all_passed:
        print("All schema validations passed successfully.")
        sys.exit(0)
    else:
        print("Some schema validations failed.")
        sys.exit(1)


if __name__ == '__main__':
    main()
