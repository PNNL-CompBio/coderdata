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
            if 'HTTPS_PROXY' in env.keys():
                res = subprocess.run(['docker','build','-t',dsname,'.','-f','build/docker/'+fn,'--build-arg','HTTPS_PROXY='+env['HTTPS_PROXY'],'--platform','linux/amd64'])
            else:
                res = subprocess.run(['docker','build','-t',dsname,'.','-f','build/docker/'+fn,'--platform','linux/amd64'])
#os.system(cmd)

    datasets = args.datasets.split(',')


    ### Any new sample creation must happened here.
    ### Each sample file requires the previous one to be created
    ### current order is : DepMap, Sanger, CPTAC, HCMI, BeatAML, MPNST
    ## can be run independently but first before omics/experiemnts
    if args.samples or args.all:
        ### build gene file
        if not os.path.exists('/tmp/genes.csv'):
            run_cmd(['genes','sh','build_genes.sh'],'gene file')
        
        ###build sample files
        sf=''
        for da in datasets:#['broad_sanger_omics','cptac','hcmi','beataml','mpnst']:
            if da=='broad_sanger':
                di = 'broad_sanger_omics'
            else:
                di = da
            if not os.path.exists('local/'+da+'_samples.csv'):
                run_cmd([di,'sh','build_samples.sh',sf],da+' samples')
            sf = '/tmp/'+da+'_samples.csv' ##TODO make this into a list

                

     ### Drug matching scripts take a while
    ### they are their own step and can be run independentyly, before others, or alongside sample/omics
    ### DepMap/Sanger, MPNST, LINCS
    dflist=[]
    if args.drugs or args.all:
        ###build drug data
        for da in [a for a in datasets if a not in ['cptac','hcmi']]:
            if da == 'broad_sanger':
                di = 'broad_sanger_exp'
            else:
                di = da
                
            if not os.path.exists('local/'+da+'_drugs.tsv'):
                run_cmd([di,'sh','build_drugs.sh',','.join(dflist)],da+' drugs')
            dflist.append('/tmp/'+da+'_drugs.tsv')

    #### Any new omics files are created here.
    ## depends on samples!
    ### these are not order dependent but require gene and sample files
    if args.omics or args.all:
        ###depmap cell line
        for da in datasets:
            if da == 'broad_sanger':
                di = 'broad_sanger_omics'
            else:
                di = da
            run_cmd([di,'sh','build_omics.sh','/tmp/genes.csv','/tmp/'+da+'_samples.csv'],da+' omics')


    ### drug response data
    ## requires samplesa nd drugs to complete
    if args.exp or args.all:
        for da in [a for a in datasets if a not in ['cptac','hcmi']]:
            if da == 'broad_sanger':
                di = 'broad_sanger_exp'
            else:
                di = da
            if not os.path.exists('local/'+da+'_experiments.tsv'):
                run_cmd([di,'sh','build_exp.sh','/tmp/'+da+'_samples.csv','/tmp/'+da+'_drugs.tsv'],da+' experiments')
    




main()
