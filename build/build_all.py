"""
script that builds the coderdata package and stores locally
"""


import os
import argparse
import subprocess

def run_cmd(cmd_arr,filename):
    '''
    short command that runs command and collates output
    '''
    print('running...'+filename)
    env = os.environ.copy()
    docker_run = ['docker','run','-v',env['PWD']+'/local/:/tmp/','-e','SYNAPSE_AUTH_TOKEN='+env['SYNAPSE_AUTH_TOKEN']+'--platform=linux/amd64']
    cmd = docker_run+cmd_arr
    print(cmd)
    res = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if res.returncode !=0:
        print(res.stderr)
        exit(filename+' file failed')
    else:
        print(filename+' retrieved')
            

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--docker',dest='docker',default=False,action='store_true')
    parser.add_argument('--samples',dest='samples',default=False,action='store_true')
    parser.add_argument('--omics',dest='omics',default=False,action='store_true')
    parser.add_argument('--drugs',dest='drugs',default=False,action='store_true')
    parser.add_argument('--exp',dest='exp',default=False,action='store_true')
    parser.add_argument('--all',dest='all',default=False,action='store_true')

    args = parser.parse_args()
                    
    ##make a 'local' directory for output
    if not os.path.exists('local'):
        os.mkdir('local')


    #docker_run='docker run -v $PWD/local/:/tmp/ --platform=linux/amd64 
    env = os.environ.copy()

    ###build docker images
    ###
    if args.docker or args.all:
        dlist = []
        for fn in os.listdir("build/docker"):
            dsname=fn.split('.')[1]
            print(dsname)
            dlist.append(dsname)
            #cmd = 'docker build --platform=linux/amd64 -t '+dsname+' .  -f build/docker/'+fn+ ' --build-arg HTTPS_PROXY=$HTTPS_PROXY'
            #print(cmd)
            res = subprocess.run(['docker','build','-t',dsname,'.','-f','build/docker/'+fn,'--build-arg','HTTPS_PROXY='+env['HTTPS_PROXY'],'--platform','linux/amd64'])
            #os.system(cmd)

    ### Any new sample creation must happened here.
    ### Each sample file requires the previous one to be created
    ### current order is : DepMap, Sanger, CPTAC, HCMI, BeatAML, MPNST
    ## can be run independently but first before omics/experiemnts
    if args.samples or args.all:
        ### build gene file
        run_cmd(['genes','sh','build_genes.sh'],'gene file')
        
        ###build sample files
        sf=''
        for di in ['broad_sanger_omics','cptac','hcmi','beataml','mpnst']:
            if di=='broad_sanger_omics':
                run_cmd([di,'sh','build_samples.sh'],di+' samples')
                sf='/tmp/broad_sanger_samples.csv'
            else:
                run_cmd([di,'sh','build_samples.sh',sf],di+' samples')
                sf = '/tmp/'+di+'_samples.csv' ##TODO make this into a list
                

     ### Drug matching scripts take a while
    ### they are their own step and can be run independentyly, before others, or alongside sample/omics
    ### DepMap/Sanger, MPNST, LINCS
    df=''
    if args.drugs or args.all:
        ###build drug data
        for di in ['broad_sanger_exp','beataml','mpnst']:
            if di=='broad_sanger_omics':
                run_cmd([di,'sh','build_drugs.sh'],di+' drugs')
                df = '/tmp/broad_sanger_drugs.csv'
            else:
                run_cmd([di,'sh','build_drugs.sh',df],di+' drugs')
                df = '/tmp/'+di+'_drugs.csv'

    #### Any new omics files are created here.
    ## depends on samples!
    ### these are not order dependent but require gene and sample files
    if args.omics or args.all:
        ###depmap cell line
        for di in ['broad_sanger_omics','beataml','mpsnst','cptac','hcmi']:
            if di=='broad_sanger_omics':
                df='broad_sanger'
            else:
                df = di
            run_cmd([di,'sh','build_omics.sh','/tmp/genes.csv','/tmp/'+df+'_samples.csv'],di+' omics')


    ### drug response data
    ## requires samplesa nd drugs to complete
    if args.exp or args.all:
        for di in ['broad_sanger_omics','beataml','mpsnst']:
            if di=='broad_sanger_omics':
                df='broad_sanger'
            else:
                df = di
            run_cmd([di,'sh','build_omics.sh','/tmp/'+df+'_samples.csv','/tmp/'+df+_'drugs.tsv'],di+' experiments')
    




main()
