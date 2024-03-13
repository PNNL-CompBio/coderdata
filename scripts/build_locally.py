"""
script that builds the coderdata package and stores locally
"""


import os
import argparse
import subprocess

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

    env = os.environ.copy()
    #docker_run='docker run -v $PWD/local/:/tmp/ --platform=linux/amd64 '
    docker_run = ['docker','run','-v',env['PWD']+'/local/:/tmp/','--platform=linux/amd64']
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
            res = subprocess.run(['docker','build','--platform','linux/amd64','-t',dsname,'.','-f','build/docker/'+fn,'--build-arg','HTTPS_PROXY='+env['HTTPS_PROXY']])
            #os.system(cmd)

    ### Any new sample creation must happened here.
    ### Each sample file requires the previous one to be created
    ### current order is : DepMap, Sanger, CPTAC, HCMI, BeatAML, MPNST
    ## can be run independently but first before omics/experiemnts
    if args.samples or args.all:
        ### build gene file
        cmd = docker_run+['depmap','Rscript','00-buildGeneFile.R']
        print(cmd)
        res = subprocess.run(cmd)
        
        ###build sample files
        cmd = docker_run+['depmap','Rscript','01-depmapSamples.R']
        print(cmd)
        res = subprocess.run(cmd)
        # os.system(cmd)

        cmd = docker_run+['depmap','Rscript','01a-pullSamples_LINCS.R','/tmp/depmap_samples.csv']
        print(cmd)
        res = subprocess.run(cmd)
        
        cmd = docker_run+['cptac','--geneFile=/tmp/genes.csv','--prevSampleFile=/tmp/depmap_samples.csv']
        print(cmd)
        res = subprocess.run(cmd)
        
        cmd =docker_run+['hcmi','python','01-createHCMISamplesFile.py','--samples=/tmp/cptac_samples.csv']
        print(cmd)
        res = subprocess.run(cmd)
        
        cmd=docker_run+['beataml','python','GetBeatAML.py','--token','$SYNAPSE_AUTH_TOKEN','--samples=/tmp/hcmi_samples.csv']
        print(cmd)
        res = subprocess.run(cmd)
        
        cmd = docker_run+['mpnst','Rscript','00_sample_gen.R','/tmp/beataml_samples.csv',env['SYNAPSE_AUTH_TOKEN']]
        print(cmd)
        res = subprocess.run(cmd)

     ### Drug matching scripts take a while
    ### they are their own step and can be run independentyly, before others, or alongside sample/omics
    ### DepMap/Sanger, MPNST, LINCS
    if args.drugs or args.all:
        ###build drug data
        dcmd=docker_run+['depmap','Rscript','03-createDrugFile.R','CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM,NCI60']
        print(dcmd)
        res = subprocess.run(dcmd)

        emd = docker_run+['mpnst','Rscript','02_get_drug_data.R',env['SYNAPSE_AUTH_TOKEN'],'/tmp/depmap/drugs.tsv']
        res = subprocess.run(emd)
        
        ##build experimental data
        ecmd= docker_run+['depmap','/opt/venv/bin/python','01b-pullDrugs_LINCS.py','--drugFile','/tmp/drugs.tsv.gz']
        print(ecmd)
        res = subprocess.run(ecmd)
        

    #### Any new omics files are created here.
    ## depends on samples!
    ### these are not order dependent but require gene and sample files
    if args.omics or args.all:
        ##omics data

        ###depmap cell line
        ocmd = docker_run+['depmap','Rscript','02-pullDepMap.R','/tmp/genes.csv','/tmp/depmap_samples.csv']
        print(ocmd)
        res = subprocess.run(ocmd)
        
        ocmd = docker_run+['depmap','Rscript','02b-pullSanger.R','/tmp/genes.csv','/tmp/depmap_samples.csv']
        print(ocmd)
        res = subprocess.run(ocmd)

        ocmd = docker_run+['depmap','/opt/venv/bin/python3','02a-depMapProts.py','--gene','/tmp/genes.csv','--sample','/tmp/depmap_samples.csv']
        print(ocmd)
        res = subprocess.run(ocmd)
        
        ###cptac
        ocmd = docker_run+['cptac','--geneFile','/tmp/genes.csv','--curSampleFile','/tmp/cptac_samples.csv']
        print(ocmd)
        res = subprocess.run(ocmd)
        
        ###HCMI - the folowing three steps are all required?
        for dt in ['transcriptomics','copy_number','mutations']:
            ocmd = docker_run+['hcmi','/opt/venv/bin/python','02-getHCMIData.py','-m','full_manifest.txt ','-t',dt,'-o','hcmi_'+dt+'.csv']
            print(ocmd)
            res = subprocess.run(ocmd)
        

        ##beataml
        ocmd=docker_run+['beataml','python','GetBeatAML.py','--token' ,env['SYNAPSE_AUTH_TOKEN'],'--samples','/tmp/beataml_samples.csv']
        print(ocmd)
        res = subprocess.run(ocmd)
        

        ###mpnst
        ocmd=docker_run+['mpnst','Rscript','01_mpnst_get_omics.R',env['SYNAPSE_AUTH_TOKEN'],'/tmp/MPNST_samples.csv','/tmp/genes.csv']
        print(ocmd)
        res = subprocess.run(ocmd)


    ### drug response data
    ## requires samplesa nd drugs to complete
    if args.exp or args.all:
            
        dcmd=docker_run+['depmap','/opt/venv/bin/python','04-drug_dosage_and_curves.py','--drugfile','/tmp/drugs.tsv','--curSampleFile','/tmp/depmap_samples.csv']
        print(dcmd)
        res = subprocess.run(dcmd)

        dcmd = docker_run+['mpnst','Rscript','03_get_drug_response_data.R',env['SYNAPSE_AUTH_TOKEN'],'/tmp/MPNST_samples.csv','/tmp/depmap/drugs.tsv']
        print(dcmd)
        res = subprocess.run(dcmd)
        

        ecmd = docker_run+['cell-line Rscript','05-LINCS_perturbations.R','/tmp/genes.csv','/tmp/drugs.tsv.gz','/tmp/depmap_samples.csv']
        res = subprocess.run(ecmd)




main()
