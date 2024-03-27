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
    docker_run = ['docker','run','-v',env['PWD']+'/local/:/tmp/','--platform=linux/amd64']
    cmd = docker_run+cmd_arr
    #print(cmd)
    res = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if res.returncode !=0:
        print(res.stderr)
        exit(filename+' samples file failed')
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
        run_cmd(['depmap_sanger','Rscript','00-buildGeneFile.R'],'gene file')
        
        ###build sample files
        run_cmd(['depmap_sanger','Rscript','01-depmap_sangerSamples.R'],'DepMap Sanger Samples')
        run_cmd(['depmap_sanger','Rscript','01a-pullSamples_LINCS.R','/tmp/depmap_sanger_samples.csv'],'LINCS samples')
        run_cmd(['cptac','--geneFile=/tmp/genes.csv','--prevSampleFile=/tmp/depmap_sanger_samples.csv'],'cptac samples')
        run_cmd(['hcmi','python','01-createHCMISamplesFile.py','--samples','/tmp/cptac_samples.csv'],'hcmi samples')
        run_cmd(['beataml','python','GetBeatAML.py','--token',env['SYNAPSE_AUTH_TOKEN'],'--samples', '--prevSamples','/tmp/hcmi_samples.csv'],'beatAML samples')
        run_cmd(['mpnst','Rscript','00_sample_gen.R','/tmp/beataml_samples.csv',env['SYNAPSE_AUTH_TOKEN']],'mpnst samples')

     ### Drug matching scripts take a while
    ### they are their own step and can be run independentyly, before others, or alongside sample/omics
    ### DepMap/Sanger, MPNST, LINCS
    if args.drugs or args.all:
        ###build drug data
        run_cmd(['depmap_sanger','Rscript','03-createDrugFile.R','CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM,NCI60'],'cell line drugs')
        run_cmd(['mpnst','Rscript','02_get_drug_data.R',env['SYNAPSE_AUTH_TOKEN'],'/tmp/depmap_sanger_drugs.tsv','/tmp/mpnst_drugs.tsv'],'mpnst drugs')
        run_cmd(['lincs','/opt/venv/bin/python','01b-pullDrugs_LINCS.py','--drugFile','/tmp/depmap_sanger_drugs.tsv'],'LINCS drugs')
        run_cmd(['beataml','python','GetBeatAML.py','--token',env['SYNAPSE_AUTH_TOKEN'], '--drugs','--drugFile','/tmp/depmap_sanger_drugs.tsv'],'BeatAML Drugs')

    #### Any new omics files are created here.
    ## depends on samples!
    ### these are not order dependent but require gene and sample files
    if args.omics or args.all:
        ###depmap cell line
        run_cmd(['depmap_sanger','Rscript','02-depmap-sanger-omics.R','/tmp/genes.csv','/tmp/depmap_sanger_samples.csv'],'depmap sanger omics')
        ###beataml
        run_cmd(['beataml','python','GetBeatAML.py','--token' ,env['SYNAPSE_AUTH_TOKEN'],'--omics','--curSamples','/tmp/beataml_samples.csv','--genes','/tmp/genes.csv'],'beatAML omics')
        ###mpnst
        run_cmd(['mpnst','Rscript','01_mpnst_get_omics.R',env['SYNAPSE_AUTH_TOKEN'],'/tmp/MPNST_samples.csv','/tmp/genes.csv'],'MPNST omics')
        ###HCMI - the folowing three steps are all required?        run_cmd(['depmap','/opt/venv/bin/python3','02a-depMapProts.py','--gene','/tmp/genes.csv','--sample','/tmp/depmap_samples.csv'],'depmap proteomics')
        ###cptac
        run_cmd(['cptac','--geneFile','/tmp/genes.csv','--curSampleFile','/tmp/cptac_samples.csv'],'cptac omics')
                ##beataml
        ###HCMI - the folowing three steps are all required?
        for dt in ['transcriptomics','copy_number','mutations']:
            run_cmd(['hcmi','python','02-getHCMIData.py','-m','full_manifest.txt','-t',dt,'-o','/tmp/hcmi_'+dt+'.csv'], 'hcmi '+dt+' omics')




    ### drug response data
    ## requires samplesa nd drugs to complete
    if args.exp or args.all:
        run_cmd(['mpnst','Rscript','03_get_drug_response_data.R',env['SYNAPSE_AUTH_TOKEN'],'/tmp/MPNST_samples.csv','/tmp/mpnst_drugs.tsv'],'MPNST experiments')
        run_cmd(['depmap_sanger','/opt/venv/bin/python','04-drug_dosage_and_curves.py','--drugfile','/tmp/depmap_sanger_drugs.tsv','--curSampleFile','/tmp/depmap_sanger_samples.csv'],'cell line experiments')
        run_cmd(['beataml','python','GetBeatAML.py','--exp','--token',env['SYNAPSE_AUTH_TOKEN'],'--curSamples','/tmp/beataml_samples.csv','--drugFile','/tmp/beataml_drugs.tsv'],'BeatAML experiments')
        run_cmd(['lincs','Rscript','05-LINCS_perturbations.R','/tmp/genes.csv','/tmp/lincs_drugs.tsv','/tmp/depmap_sanger_samples.csv'],'LINCS perturbations')
        




main()
