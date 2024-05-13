"""
script that builds the coderdata package and stores locally
"""

import os
import argparse
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor

    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--docker',dest='docker',default=False,action='store_true')
    parser.add_argument('--samples',dest='samples',default=False,action='store_true')
    parser.add_argument('--omics',dest='omics',default=False,action='store_true')
    parser.add_argument('--drugs',dest='drugs',default=False,action='store_true')
    parser.add_argument('--exp',dest='exp',default=False,action='store_true')
    parser.add_argument('--all',dest='all',default=False,action='store_true')
    parser.add_argument('--dataset',dest='datasets',default='broad_sanger,hcmi,beataml,mpnst,cptac',help='Datasets to process. Defaults to all available, but if there are synapse issues, please remove beataml and mpnst')

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
            docker_run = ['docker','run','-v',env['PWD']+'/local/:/tmp/','--platform=linux/amd64']
        else:
            docker_run = ['docker','run','-v',env['PWD']+'/local/:/tmp/','-e','SYNAPSE_AUTH_TOKEN='+env['SYNAPSE_AUTH_TOKEN'],'--platform=linux/amd64']
        cmd = docker_run+cmd_arr
        print(cmd)
        res = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if res.returncode !=0:
            print(res.stderr)
            exit(filename+' file failed')
        else:
            print(filename+' retrieved')
        
    
    def process_docker():
        '''
        Build all docker images using docker compose
        All output and errors are logged at local/docker.log
        '''
        compose_file = 'build/docker/docker-compose.yml'
        compose_command = ['docker-compose', '-f', compose_file, 'build', '--parallel']
        log_file_path = 'local/docker.log'
        env = os.environ.copy()
        print(f"Docker-compose is building all images. View output in {log_file_path}.")
        with open(log_file_path, 'w') as log_file:
            # Execute the docker-compose command
            res = subprocess.run(compose_command,  env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            # Log both stdout and stderr to the log file
            log_file.write(res.stdout)
            if res.returncode != 0:
                log_file.write("Docker compose build failed.\n")
                print(f"Docker compose build failed. See {log_file_path} for details.")
                exit(1)
            else:
                log_file.write("Docker images built successfully.\n")
                print(f"Docker images built successfully. Details logged in {log_file_path}")


    def process_drugs(executor, datasets):
        '''
        Build all drug files sequentially
        '''
        last_drug_future = None
        dflist = []  
            
        # WE NEED A METHOD TO CONFIRM THAT DRUG FILES ARE NOT INCOMPLETE
            
        # Check for existing files and update dflist with processed files
        for da in datasets:
            if da not in ['cptac', 'hcmi']: 
                file_path = f'local/{da}_drugs.tsv'
                if os.path.exists(file_path):
                    dflist.append(f'/tmp/{da}_drugs.tsv')  # Add to dflist if already processed

        for da in datasets:
            if da not in ['cptac', 'hcmi']:
                di = 'broad_sanger_exp' if da == 'broad_sanger' else da
                if not os.path.exists(f'local/{da}_drugs.tsv'):
                    if last_drug_future:
                        last_drug_future.result()  # Ensure the last drug process is completed before starting the next
                    last_drug_future = executor.submit(run_docker_cmd, [di, 'sh', 'build_drugs.sh', ','.join(dflist)], f'{da} drugs')
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
                    last_sample_future = executor.submit(run_docker_cmd, [di, 'sh', 'build_samples.sh', sf], f'{da} samples')
                    sf = f'/tmp/{da}_samples.csv'
                    
    def process_omics(executor, datasets):
        '''
        Build all omics files concurrently
        '''
        last_sample_future = None
        for da in datasets:
            di = 'broad_sanger_omics' if da == 'broad_sanger' else da
            
            #Run all at once:
            # executor.submit(run_docker_cmd, [di, 'sh', 'build_omics.sh', '/tmp/genes.csv', f'/tmp/{da}_samples.csv'], f'{da} omics')

            #Run 1 at a time.
            if last_sample_future:
                last_sample_future.result() 
            last_sample_future = executor.submit(run_docker_cmd, [di, 'sh', 'build_omics.sh', '/tmp/genes.csv', f'/tmp/{da}_samples.csv'], f'{da} omics')
    
    def process_experiments(executor, datasets):
        '''
        Build all experiments files concurrently
        '''
        for da in datasets:
            if da not in ['cptac', 'hcmi']:
                di = 'broad_sanger_exp' if da == 'broad_sanger' else da
                if not os.path.exists(f'local/{da}_experiments.tsv'):
                    
                    #Run all at once
                    # executor.submit(run_docker_cmd, [di, 'sh', 'build_exp.sh', f'/tmp/{da}_samples.csv', f'/tmp/{da}_drugs.tsv'], f'{da} experiments')
                    
                    #Run one at a time
                    if last_sample_future:
                        last_sample_future.result() 
                    last_sample_future = executor.submit(run_docker_cmd, [di, 'sh', 'build_exp.sh', f'/tmp/{da}_samples.csv', f'/tmp/{da}_drugs.tsv'], f'{da} experiments')
        



    def process_genes(executor):
        if not os.path.exists('/tmp/genes.csv'):
            executor.submit(run_docker_cmd,['genes','sh','build_genes.sh'],'gene file')
        
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
        process_docker()
        

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
        if args.samples or args.all:
            sample_thread.result()
        if args.samples or args.omics or args.exp or args.all:
            gene_thread.result()


    ### At this point in the pipeline, all samples and drugs files have been created. There are no blockers to proceed.
    ### Build Omics files and Experiments files. These two steps are run in Parallel. 
    ### Within each step, all datasets are run in Parallel.
    
    with ThreadPoolExecutor() as executor:
        if args.omics or args.all:
            omics_thread = executor.submit(process_omics, executor, datasets)
        if args.exp or args.all:
            exp_thread = executor.submit(process_experiments, executor, datasets)
            
        if args.omics or args.all:
            omics_thread.result()
        if args.exp or args.all:
            exp_thread.result()
    
    with ThreadPoolExecutor() as executor:
        if args.omics or args.all:
            omics_thread = executor.submit(process_omics, executor, datasets)
        if args.exp or args.all:
            exp_thread = executor.submit(process_experiments, executor, datasets)
            
        if args.omics or args.all:
            omics_thread.result()
        if args.exp or args.all:
            exp_thread.result()


if __name__ == '__main__':
    main()