"""
script that builds the coderdata package and stores locally
"""

import os
import argparse
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor
import shutil
import gzip
from glob import glob
import sys
    
def main():
    parser=argparse.ArgumentParser(
        description="This script initializes all docker containers, builds datasets, validates them, and uploads to Figshare.",
        epilog="""Examples of usage:

Build all datasets in a high memory environment, validate them, and upload to Figshare:
  python build/build_all.py --all --high_mem --validate --figshare --version 0.1.29

Build only experiment files. This assumes preceding steps (docker images, samples, omics, and drugs) have already been completed:
  python build/build_all.py --exp

Validate all local files without building or uploading. These files must be located in ./local. Includes compression/decompression steps.
  python build/build_all.py --validate

Upload the latest data to Figshare (ensure tokens are set in the local environment):
  python build/build_all.py --figshare --version 0.1.30
        """
    )
    parser.add_argument('--docker',dest='docker',default=False,action='store_true', help="Build all docker images.")
    parser.add_argument('--samples',dest='samples',default=False,action='store_true', help="Build all sample files.")
    parser.add_argument('--omics',dest='omics',default=False,action='store_true', help="Build all omics files.")
    parser.add_argument('--drugs',dest='drugs',default=False,action='store_true', help="Build all drug files")
    parser.add_argument('--exp',dest='exp',default=False,action='store_true', help="Build all experiment file.")
    parser.add_argument('--validate', action='store_true', help="Run schema checker on all local files. Note this will be run, whether specified or not, if figshare arguments are included.")
    parser.add_argument('--figshare', action='store_true', help="Upload all local data to Figshare. FIGSHARE_TOKEN must be set in local environment.")
    parser.add_argument('--all',dest='all',default=False,action='store_true', help="Run all data build commands. This includes docker, samples, omics, drugs, exp arguments. This does not run the validate or figshare commands")
    parser.add_argument('--high_mem',dest='high_mem',default=False,action='store_true',help = "If you have 32 or more CPUs, this option is recommended. It will run many code portions in parallel. If you don't have enough memory, this will cause a run failure.")
    parser.add_argument('--dataset',dest='datasets',default='broad_sanger,hcmi,beataml,cptac,mpnst,mpnstpdx',help='Datasets to process. Defaults to all available.')
    parser.add_argument('--version', type=str, required=False, help='Version number for the Figshare upload title (e.g., "0.1.29"). This is required for Figshare upload. This must be a higher version than previously published versions.')
    parser.add_argument('--github-username', type=str, required=False, help='GitHub username for the repository.')
    parser.add_argument('--github-email', type=str, required=False, help='GitHub email for the repository.')
    
    args = parser.parse_args()
                    
    # Simulation command for testing order of everything:
    # def run_docker_cmd(cmd_arr, filename):
    #     # Simulate running the command by printing what would be run
    #     print(f'Running: {filename} with command {" ".join(cmd_arr)}')
    #     # Simulate execution time with a random delay
    #     time.sleep(2)
    #     print(f'Completed: {filename}')
    
    def run_docker_cmd(cmd_arr,filename):
        '''
        Essentially a wrapper for 'docker run'. Also provides output.
        '''
        print('running...'+filename)
        env = os.environ.copy()
        if 'SYNAPSE_AUTH_TOKEN' not in env.keys():
            print('You need to set the SYNAPSE_AUTH_TOKEN to acess the MPNST and beatAML Datasets')
            docker_run = ['docker','run','--rm','-v',env['PWD']+'/local/:/tmp/','--platform=linux/amd64']
        else:
            docker_run = ['docker','run','--rm','-v',env['PWD']+'/local/:/tmp/','-e','SYNAPSE_AUTH_TOKEN='+env['SYNAPSE_AUTH_TOKEN'],'--platform=linux/amd64']
            
            
        cmd = docker_run+cmd_arr
        print(cmd)
        # res = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        res = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
        if res.returncode !=0:
            print(res.stderr)
            exit(filename+' file failed')
        else:
            print(filename+' retrieved')
        
    
    # def process_docker():
    #     '''
    #     Build all docker images using docker compose
    #     All output and errors are logged at local/docker.log
    #     '''
    #     compose_file = 'build/docker/docker-compose.yml'
    #     compose_command = ['docker-compose', '-f', compose_file, 'build', '--parallel']
    #     log_file_path = 'local/docker.log'
    #     env = os.environ.copy()
    #     print(f"Docker-compose is building all images. View output in {log_file_path}.")
    #     with open(log_file_path, 'w') as log_file:
    #         # Execute the docker-compose command
    #         res = subprocess.run(compose_command,  env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    #         # Log both stdout and stderr to the log file
    #         log_file.write(res.stdout)
    #         if res.returncode != 0:
    #             log_file.write("Docker compose build failed.\n")
    #             print(f"Docker compose build failed. See {log_file_path} for details.")
    #             exit(1)
    #         else:
    #             log_file.write("Docker images built successfully.\n")
    #             print(f"Docker images built successfully. Details logged in {log_file_path}")


    def process_docker(datasets):
        '''
        Build specific docker images using docker-compose based on the dataset argument.
        All output and errors are logged at local/docker.log.
        
        Parameters:
        - datasets: list of datasets to process (e.g., ['broad_sanger', 'hcmi', 'mpnst'])
        '''
        compose_file = 'build/docker/docker-compose.yml'
        
        # Map datasets to corresponding Docker Containers
        dataset_map = {
            'broad_sanger': ['broad_sanger_exp', 'broad_sanger_omics'],
            'hcmi': ['hcmi'],
            'beataml': ['beataml'],
            'mpnst': ['mpnst'],
            'mpnstpdx': ['mpnstpdx'],
            'cptac': ['cptac'],
            'genes': ['genes'],
            'upload': ['upload']
        }
        
        # Collect container names to build based on the datasets provided. Always build genes and upload.
        datasets_to_build = ['genes', 'upload']
        for dataset in datasets:
            datasets_to_build.extend(dataset_map.get(dataset, []))
        
        # Build the docker-compose command, adding specific datasets
        compose_command = ['docker-compose', '-f', compose_file, 'build', '--parallel'] + datasets_to_build
        
        log_file_path = 'local/docker.log'
        env = os.environ.copy()
        
        print(f"Docker-compose is building images for {', '.join(datasets_to_build)}. View output in {log_file_path}.")
        
        with open(log_file_path, 'w') as log_file:
            try:
                # Execute the docker-compose command
                res = subprocess.run(compose_command, env=env, stdout=log_file, stderr=log_file, text=True, check=True)
                log_file.write("Docker images built successfully.\n")
                print(f"Docker images for {', '.join(datasets_to_build)} built successfully. Details logged in {log_file_path}.")
            except subprocess.CalledProcessError as e:
                log_file.write(f"Docker compose build failed with error: {e}\n")
                print(f"Docker compose build failed. See {log_file_path} for details.")
                raise

    def process_drugs(executor, datasets):
        '''
        Build all drug files sequentially
        '''
        last_drug_future = None
        dflist = []  
            
        # WE NEED A METHOD TO CONFIRM THAT DRUG FILES ARE NOT INCOMPLETE
        ##THIS IS BUILT IN- always rerun drug code to check
        # Check for existing files and update dflist with processed files
        for da in datasets:
            if da not in ['cptac', 'hcmi']: 
                file_path = f'local/{da}_drugs.tsv'
                desc_path = f'local/{da}_drug_descriptor.tsv.gz'
                #if os.path.exists(file_path): ##always rerun drug process
                #    dflist.append(f'/tmp/{da}_drugs.tsv')  # Add to dflist if already processed

        for da in datasets:
            if da not in ['cptac', 'hcmi']:
                di = 'broad_sanger_exp' if da == 'broad_sanger' else da
                #if not os.path.exists(f'local/{da}_drugs.tsv'):
                if last_drug_future:
                    last_drug_future.result()  # Ensure the last drug process is completed before starting the next
                last_drug_future = executor.submit(run_docker_cmd, [di, 'bash', 'build_drugs.sh', ','.join(dflist)], f'{da} drugs')
                dflist.append(f'/tmp/{da}_drugs.tsv')
    
    def process_samples(executor, datasets):
        '''
        Build all samples files sequentially
        '''
        last_sample_future = None
        sf = ''
        for da in datasets:
            file_path = f'local/{da}_samples.csv'
            if os.path.exists(file_path):
                sf = f'/tmp/{da}_samples.csv'  # Set the most recent successfully processed dataset file
            else:
                break
        for da in datasets:
                di = 'broad_sanger_omics' if da == 'broad_sanger' else da
                if not os.path.exists(f'local/{da}_samples.csv'):
                    if last_sample_future:
                        last_sample_future.result() 
                    last_sample_future = executor.submit(run_docker_cmd, [di, 'bash', 'build_samples.sh', sf], f'{da} samples')
                    sf = f'/tmp/{da}_samples.csv'
                    
    def process_omics(executor, datasets,high_mem):
        '''
        Build all omics files concurrently
        '''
        last_omics_future = None
        for da in datasets:
            di = 'broad_sanger_omics' if da == 'broad_sanger' else da
            #Run all at once:
            if high_mem:
                executor.submit(run_docker_cmd, [di, 'bash', 'build_omics.sh', '/tmp/genes.csv', f'/tmp/{da}_samples.csv'], f'{da} omics')
            #Run one at a time.
            else:
                if last_omics_future:
                    last_omics_future.result() 
                last_omics_future = executor.submit(run_docker_cmd, [di, 'bash', 'build_omics.sh', '/tmp/genes.csv', f'/tmp/{da}_samples.csv'], f'{da} omics')
        
    def process_experiments(executor, datasets, high_mem):
        '''
        Build all experiments files concurrently
        '''
        last_experiments_future = None
        for da in datasets:
            if da not in ['cptac', 'hcmi']:
                di = 'broad_sanger_exp' if da == 'broad_sanger' else da
                if not os.path.exists(f'local/{da}_experiments.tsv'):
                    #Run all at once
                    if high_mem:
                        executor.submit(run_docker_cmd, [di, 'bash', 'build_exp.sh', f'/tmp/{da}_samples.csv', f'/tmp/{da}_drugs.tsv'], f'{da} experiments')
                    #Run one at a time
                    else:
                        if last_experiments_future:
                            last_experiments_future.result() 
                        last_experiments_future = executor.submit(run_docker_cmd, [di, 'bash', 'build_exp.sh', f'/tmp/{da}_samples.csv', f'/tmp/{da}_drugs.tsv'], f'{da} experiments')


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
            #Running the build_misc.sh in broad_sanger_omics
            di = 'broad_sanger_omics' if da == 'broad_sanger' else da
            #Run all at once:
            if high_mem:
                executor.submit(run_docker_cmd, [di, 'bash', 'build_misc.sh'], f'{da} misc')
            #Run one at a time.
            else:
                if last_misc_future:
                    last_misc_future.result() 
                last_misc_future = executor.submit(run_docker_cmd,  [di, 'bash', 'build_misc.sh'], f'{da} misc')
        

    def process_genes(executor):
        if not os.path.exists('/tmp/genes.csv'):
            executor.submit(run_docker_cmd,['genes','bash','build_genes.sh'],'gene file')
        
        
    def run_docker_upload_cmd(cmd_arr, all_files_dir, name, version):
        '''
        Wrapper for 'docker run'. This one is focused on uploads.
        '''
        env = os.environ.copy()
        docker_run = ['docker', 'run', '--rm', '-v', f"{env['PWD']}/local/{all_files_dir}:/tmp", '-e', f"VERSION={version}"]

        # Add Appropriate Environment Variables
        if 'FIGSHARE_TOKEN' in env and name == 'Figshare':
            docker_run.extend(['-e', f"FIGSHARE_TOKEN={env['FIGSHARE_TOKEN']}", 'upload'])
        if name == "validate":
            docker_run.extend(['upload'])
        if 'GITHUB_TOKEN' in env and name == "GitHub":
            docker_run.extend(['-e', f"GITHUB_TOKEN={env['GITHUB_TOKEN']}", 'upload'])

        # Full command to run including version update
        docker_run.extend(cmd_arr)
        print('Executing:', ' '.join(docker_run))
        # res = subprocess.run(docker_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res = subprocess.run(docker_run, stdout=sys.stdout, stderr=sys.stderr)
        if res.returncode != 0:
            print(res.stderr)
            exit(f'{name} failed')
        else:
            print(f'{name} successful')
            

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
            
    ######
    ### Pre-Build Environment Token Check
    #####

    figshare_token = os.getenv('FIGSHARE_TOKEN')
    synapse_auth_token = os.getenv('SYNAPSE_AUTH_TOKEN')
    github_token = os.getenv('GITHUB_TOKEN')


    # Error handling for required tokens
    if args.figshare and not figshare_token:
        raise ValueError("FIGSHARE_TOKEN environment variable is not set.")
    if ('beataml' in args.datasets or 'mpnst' in args.datasets) and not synapse_auth_token:
        if args.docker or args.samples or args.omics or args.drugs or args.exp or args.all: # Token only required if building data, not upload or validate.
            raise ValueError("SYNAPSE_AUTH_TOKEN is required for accessing MPNST and beatAML datasets.")
       
    ######
    ### Begin Pipeline
    #####
    
    # Make a 'local' directory for output
    if not os.path.exists('local'):
        os.mkdir('local')

    # Get dataset names - default is all.
    datasets = args.datasets.split(',')

    ### Build Docker Images. These are all built in Parallel. Nothing else can run until these are built.
    # Ouput is logged at local/docker.log
    if args.docker or args.all:
        process_docker(datasets)
        print("Docker image generation completed")
        

    ### Build Drugs files, Samples files, and Genes file. These two steps are run in Parallel. 
    ### Within each step, sequential running is required.
    with ThreadPoolExecutor() as executor:
        if args.samples or args.all:
            sample_thread = executor.submit(process_samples,executor, datasets)
        if args.drugs or args.all:
            drug_thread = executor.submit(process_drugs,executor, datasets)
        if args.samples or args.omics or args.exp or args.all:
            gene_thread = executor.submit(process_genes,executor)
        
        # Wait for both processes to complete before proceeding to omics and experiments
        if args.drugs or args.all:
            drug_thread.result()
        if args.samples or args.all:##need to wait for samples for all of these
            sample_thread.result()
        if args.samples or args.omics or args.exp or args.all:
            gene_thread.result()
            
    print("All samples, drugs files, and genes file completed or skipped")


    ### At this point in the pipeline, all samples and drugs files have been created. There are no blockers to proceed.
    ### Build Omics files and Experiments files. These two steps are run in Parallel. 
    ### Within each step, all datasets are run in Parallel.
    
    with ThreadPoolExecutor() as executor:
        if args.omics or args.all:
            omics_thread = executor.submit(process_omics, executor, datasets, args.high_mem)
        if args.exp or args.all:
            exp_thread = executor.submit(process_experiments, executor, datasets, args.high_mem)
            
        if args.omics or args.all:
            omics_thread.result()
            print("All omics files completed")
        if args.exp or args.all:
            exp_thread.result()
            print("All experiments files completed")


    ### Final Step, some datasets may need an additional post build step. Add this here
    # Currently only the cell line datasets need this. This seperates broad_sanger into all of its component datasets.
    
    with ThreadPoolExecutor() as executor:
        if args.all:
            misc_thread = executor.submit(process_misc, executor, datasets, args.high_mem)
        if args.all:
            misc_thread.result()
            print("Final build step complete.")


    ######
    ### Begin Upload and/or validation
    #####
    
    if args.figshare or args.validate:
        # FigShare File Prefixes:
        prefixes = ['beataml', 'hcmi', 'cptac', 'mpnst', 'genes', 'drugs']
        broad_sanger_datasets = ["ccle","ctrpv2","fimm","gdscv1","gdscv2","gcsi","prism","nci60"]
        if "broad_sanger" in datasets:
            prefixes.extend(broad_sanger_datasets)
            datasets.extend(broad_sanger_datasets)
            datasets.remove("broad_sanger")

        
        figshare_token = os.getenv('FIGSHARE_TOKEN')

        all_files_dir = 'local/all_files_dir'
        if not os.path.exists(all_files_dir):
            os.makedirs(all_files_dir)
        
        # Ensure figshare tokens are available
        if  args.figshare and not figshare_token:
            raise ValueError("Required tokens (FIGSHARE) are not set in environment variables.")
        
        # Ensure version is specified
        if args.figshare and not args.version:
            raise ValueError("Version must be specified when pushing to figshare")

        # Move relevant files to a designated directory
        for file in glob(os.path.join("local", '*.*')):
            if any(file.startswith(os.path.join("local", prefix)) for prefix in prefixes):
                shutil.move(file, os.path.join(all_files_dir, os.path.basename(file)))

        # Decompress all compressed files in the directory for schema checking
        for file in glob(os.path.join(all_files_dir, '*.gz')):
            decompress_file(file)

        # Run schema checker - This will always run if uploading data.
        schema_check_command = ['python3', 'scripts/check_schema.py', '--datasets'] + datasets
        run_docker_upload_cmd(schema_check_command, 'all_files_dir', 'validate', args.version)
        
        print("Validation complete. Proceeding with file compression/decompression adjustments")
        
        # Compress or decompress files based on specific conditions after checking
        for file in glob(os.path.join(all_files_dir, '*')):
            is_compressed = file.endswith('.gz')
            if ('samples' in file or 'figshare' in file) and is_compressed:
                decompress_file(file)
            elif not ('samples' in file or 'figshare' in file) and not is_compressed:
                compress_file(file)

        print("File compression and decompression adjustments are complete.")
    
    # Upload to Figshare using Docker
        if args.figshare and args.version and figshare_token:
            figshare_command = ['python3', 'scripts/push_to_figshare.py', '--directory', "/tmp", '--title', f"CODERData{args.version}", '--token', os.getenv('FIGSHARE_TOKEN'), '--project_id', '189342', '--publish']
            run_docker_upload_cmd(figshare_command, 'all_files_dir', 'Figshare', args.version)

            
            # Push changes to GitHub using Docker
        if args.version and args.figshare and figshare_token and github_token and args.github_username and args.github_email:
            git_command = [
                'bash', '-c', (
                    f'git config --global user.name "{args.github_username}" '
                    f'&& git config --global user.email "{args.github_email}" '
                    f'&& cp /tmp/figshare_latest.yml /usr/src/app/coderdata/docs/_data/figshare_latest.yml '
                    f'&& git add docs/_data/figshare_latest.yml '
                    f'&& git commit -m "Data Built and Uploaded. New Tag: {args.version}" '
                    f'&& git tag {args.version} '
                    f'&& git push https://{args.github_username}:{github_token}@github.com/PNNL-CompBio/coderdata.git main '
                    f'&& git push https://{args.github_username}:{github_token}@github.com/PNNL-CompBio/coderdata.git --tags'
                )
            ]
            run_docker_upload_cmd(git_command, 'all_files_dir', 'GitHub', args.version)
            
    
if __name__ == '__main__':
    main()
