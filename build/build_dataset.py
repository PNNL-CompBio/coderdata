"""
Script that builds a single dataset.
"""

import os
import argparse
import subprocess
import shutil
import gzip
from concurrent.futures import ThreadPoolExecutor
import glob 

def run_docker_cmd(cmd_arr, filename):
    '''
    Wrapper for 'docker run' command. Executes a Docker container with the specified command.
    '''
    print('Running...', filename)
    env = os.environ.copy()
    if 'SYNAPSE_AUTH_TOKEN' not in env:
        print('You need to set the SYNAPSE_AUTH_TOKEN to access the MPNST and beatAML datasets')
        docker_run = ['docker', 'run', '-v', f"{env['PWD']}/local/:/tmp/", '--platform=linux/amd64']
    else:
        docker_run = ['docker', 'run',  '-v', f"{env['PWD']}/local/:/tmp/", '-e', f"SYNAPSE_AUTH_TOKEN={env['SYNAPSE_AUTH_TOKEN']}", '--platform=linux/amd64']

    cmd = docker_run + cmd_arr
    print('Executing command:', ' '.join(cmd))
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stderr.decode())
        exit(f'{filename} failed')
    else:
        print(f'{filename} completed successfully')

def process_docker(dataset,validate):
    '''
    Build Docker images required for the specified dataset.
    '''
    compose_file = 'build/docker/docker-compose.yml'
    dataset_map = {
        'broad_sanger': ['broad_sanger_exp', 'broad_sanger_omics'],
        'hcmi': ['hcmi'],
        'beataml': ['beataml'],
        'mpnst': ['mpnst'],
        'cptac': ['cptac'],
        'genes': ['genes'],
        'upload': ['upload']
    }

    # Collect container names to build based on the dataset provided. Always build 'genes'.
    datasets_to_build = ['genes']
    # Append upload if validation step is included
    if validate is True:
        datasets_to_build.append('upload')
        
    datasets_to_build.extend(dataset_map.get(dataset, []))

    compose_command = ['docker-compose', '-f', compose_file, 'build'] + datasets_to_build

    log_file_path = 'local/docker.log'
    env = os.environ.copy()

    print(f"Docker-compose is building images for {', '.join(datasets_to_build)}. View output in {log_file_path}.")

    with open(log_file_path, 'w') as log_file:
        try:
            subprocess.run(compose_command, env=env, stdout=log_file, stderr=log_file, text=True, check=True)
            log_file.write("Docker images built successfully.\n")
            print(f"Docker images for {', '.join(datasets_to_build)} built successfully. Details logged in {log_file_path}.")
        except subprocess.CalledProcessError as e:
            log_file.write(f"Docker compose build failed with error: {e}\n")
            print(f"Docker compose build failed. See {log_file_path} for details.")
            raise

def process_genes(executor):
    '''
    Build the genes file if it does not exist.
    '''
    if not os.path.exists('local/genes.csv'):
        executor.submit(run_docker_cmd, ['genes', 'sh', 'build_genes.sh'], 'genes file')

def process_samples(executor, dataset, use_prev_dataset, should_continue):
    '''
    Build the samples file for the specified dataset.
    '''
    samples_file = f'local/{dataset}_samples.csv'
    if should_continue and os.path.exists(samples_file):
        print(f"Samples file for {dataset} already exists. Skipping samples build.")
        return

    prev_samples_file = f'/tmp/{use_prev_dataset}_samples.csv' if use_prev_dataset else ''
    di = 'broad_sanger_omics' if dataset == 'broad_sanger' else dataset
    filename = f'{dataset} samples'
    executor.submit(run_docker_cmd, [di, 'sh', 'build_samples.sh', prev_samples_file], filename)

def process_drugs(executor, dataset, use_prev_dataset, should_continue):
    '''
    Build the drugs file for the specified dataset.
    '''
    if dataset in ['cptac', 'hcmi']:
        return  # No drugs to process for these datasets

    drugs_file = f'local/{dataset}_drugs.tsv'
    if should_continue and os.path.exists(drugs_file):
        print(f"Drugs file for {dataset} already exists. Skipping drugs build.")
        return

    prev_drugs_file = f'/tmp/{use_prev_dataset}_drugs.tsv' if use_prev_dataset else ''
    dflist = [prev_drugs_file] if use_prev_dataset else []
    di = 'broad_sanger_exp' if dataset == 'broad_sanger' else dataset
    filename = f'{dataset} drugs'
    executor.submit(run_docker_cmd, [di, 'sh', 'build_drugs.sh', ','.join(dflist)], filename)


def process_omics(executor, dataset, should_continue):
    '''
    Build the omics files for the specified dataset.
    '''
    # Map datasets to their expected omics files
    dataset_omics_files = {
        'beataml': ['mutations', 'proteomics', 'transcriptomics'],
        'mpnst': ['copy_number', 'mutations', 'proteomics', 'transcriptomics'],
        'broad_sanger': ['copy_number', 'mutations', 'proteomics', 'transcriptomics'],
        'cptac': ['copy_number', 'mutations', 'proteomics', 'transcriptomics'],
        'hcmi': ['mutations', 'transcriptomics']
    }

    expected_omics = dataset_omics_files.get(dataset, [])

    if not expected_omics:
        print(f"No omics data expected for dataset {dataset}. Skipping omics build.")
        return

    # Check if all expected omics files exist
    omics_files_exist = True
    for omics_type in expected_omics:
        patterns = [
            f'local/{dataset}_{omics_type}.csv',
            f'local/{dataset}_{omics_type}.csv.gz',
            f'local/{dataset}_{omics_type}.tsv',
            f'local/{dataset}_{omics_type}.tsv.gz'
        ]
        file_found = False
        for pattern in patterns:
            matches = glob.glob(pattern)
            if matches:
                file_found = True
                break
        if not file_found:
            omics_files_exist = False
            break  # If any omics files are missing, just build / rebuild them all.

    if should_continue and omics_files_exist:
        print(f"Omics files for {dataset} already exist. Skipping omics build.")
        return

    di = 'broad_sanger_omics' if dataset == 'broad_sanger' else dataset
    filename = f'{dataset} omics'
    executor.submit(run_docker_cmd, [di, 'sh', 'build_omics.sh', '/tmp/genes.csv', f'/tmp/{dataset}_samples.csv'], filename)


def process_experiments(executor, dataset, should_continue):
    '''
    Build the experiments files for the specified dataset.
    '''
    if dataset in ['cptac', 'hcmi']:
        return  # No experiments to process for these datasets

    experiments_file = f'local/{dataset}_experiments.tsv'
    if should_continue and os.path.exists(experiments_file):
        print(f"Experiments file for {dataset} already exists. Skipping experiments build.")
        return

    di = 'broad_sanger_exp' if dataset == 'broad_sanger' else dataset
    filename = f'{dataset} experiments'
    executor.submit(run_docker_cmd, [di, 'sh', 'build_exp.sh', f'/tmp/{dataset}_samples.csv', f'/tmp/{dataset}_drugs.tsv'], filename)



def process_misc(executor, datasets, high_mem):
    '''
    Run all misc scripts concurrently or one at a time.
    '''
    last_misc_future = None
    #Currently this only applies to broad_sanger. Add others here if they need a final step.
    if "broad_sanger" in datasets:
        datasets = ["broad_sanger"]
    else:
        return
    for da in datasets:
        di = 'broad_sanger_omics' if da == 'broad_sanger' else da
        #Run all at once:
        if high_mem:
            executor.submit(run_docker_cmd, [di, 'sh', 'build_misc.sh'], f'{da} misc')
        #Run one at a time.
        else:
            if last_misc_future:
                last_misc_future.result() 
            last_misc_future = executor.submit(run_docker_cmd,  [di, 'sh', 'build_misc.sh'], f'{da} misc')
    


def decompress_file(file_path):
    """Decompress a gzip file and delete the original compressed file."""
    with gzip.open(file_path, 'rb') as f_in:
        decompressed_file_path = file_path[:-3]  # Remove '.gz' from the filename
        with open(decompressed_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)

def compress_file(file_path):
    """Compress a file using gzip and delete the original uncompressed file."""
    compressed_file_path = file_path + '.gz'
    with open(file_path, 'rb') as f_in:
        with gzip.open(compressed_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)

def run_docker_validate_cmd(cmd_arr, all_files_dir, name):
    '''
    Wrapper for 'docker run' command used during validation and uploads.
    '''
    env = os.environ.copy()
    docker_run = ['docker', 'run', '-v', f"{env['PWD']}/local/{all_files_dir}:/tmp"]
    docker_run.extend(['upload']) 
    docker_run.extend(cmd_arr)
    print('Executing:', ' '.join(docker_run))
    res = subprocess.run(docker_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stderr.decode())
        exit(f'{name} failed')
    else:
        print(f'{name} completed successfully')

def run_schema_checker(dataset):
    '''
    Run schema checker on the built files for the specified dataset.
    '''
    # Prepare the directory with the built files
    prefixes = ['genes', dataset]
    datasets = [dataset]
    broad_sanger_datasets = ["ccle","ctrpv2","fimm","gdscv1","gdscv2","gcsi","prism","nci60"]
    all_files_dir = 'all_files_dir'
    if "broad_sanger" == dataset:
            prefixes.extend(broad_sanger_datasets)
            datasets.extend(broad_sanger_datasets)
            datasets.remove("broad_sanger")
            prefixes.remove("broad_sanger")
    
    if not os.path.exists(f'local/{all_files_dir}'):
        os.makedirs(f'local/{all_files_dir}')

    # Move relevant files to all_files_dir
    for file in os.listdir('local'):
        if any(file.startswith(prefix) for prefix in prefixes):
            shutil.move(os.path.join('local', file), os.path.join('local', all_files_dir, file))

    # Decompress any compressed files
    for file in os.listdir(f'local/{all_files_dir}'):
        if file.endswith('.gz'):
            decompress_file(os.path.join('local', all_files_dir, file))

    # Run schema checker
    schema_check_command = ['python3', 'scripts/check_all_schemas.py', '--datasets'] + datasets
    run_docker_validate_cmd(schema_check_command, all_files_dir, 'Validation')

def main():
    parser = argparse.ArgumentParser(
        description="This script builds a single dataset."
    )
    parser.add_argument('--dataset', required=True, help='Name of the dataset to build')
    parser.add_argument('--use_prev_dataset', help='Prefix of the previous dataset for sample and drug ID assignment')
    parser.add_argument('--build', action='store_true', help='Run data build.')
    parser.add_argument('--validate', action='store_true', help='Run schema checker on the built files')
    parser.add_argument('--continue', dest='should_continue', action='store_true', help='Continue from where the build left off by skipping existing files')

    args = parser.parse_args()

    if not os.path.exists('local'):
        os.mkdir('local')
        
    # Build Docker Image
    process_docker(args.dataset,args.validate)

    if args.build:
        # Use ThreadPoolExecutor for parallel execution
        with ThreadPoolExecutor() as executor:
            # Always build genes file
            process_genes(executor)

            # Build samples and drugs
            samples_future = executor.submit(process_samples, executor, args.dataset, args.use_prev_dataset, args.should_continue)
            drugs_future = executor.submit(process_drugs, executor, args.dataset, args.use_prev_dataset, args.should_continue)

            samples_future.result()
            drugs_future.result()
            
        print("Samples and Drugs Files Completed.")

        with ThreadPoolExecutor() as executor:
            
            # Build omics and experiments
            omics_future = executor.submit(process_omics, executor, args.dataset, args.should_continue)
            experiments_future = executor.submit(process_experiments, executor, args.dataset, args.should_continue)

            omics_future.result()
            experiments_future.result()

        print("Experiments and Omics Files completed.")
        
        with ThreadPoolExecutor() as executor:
            
            if args.all:
                misc_thread = executor.submit(process_misc, executor, args.dataset, args.high_mem)
            if args.all:
                misc_thread.result()
                print("Final build step complete.")
                
    if args.validate:
        run_schema_checker(args.dataset)
        print("Validation completed.")

if __name__ == '__main__':
    main()
