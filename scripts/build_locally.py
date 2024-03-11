"""
script that builds the coderdata package and stores locally
"""


import os
import argparse


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
    
    docker_run='docker run -v $PWD/local:/tmp '
    ###build docker images
    ###
    if args.docker or args.all:
        dlist = []
        for fn in os.listdir("build/docker"):
            dsname=fn.split('.')[1]
            print(dsname)
            dlist.append(dsname)
            cmd = 'docker build -t '+dsname+' .  -f build/docker/'+fn+ ' --build-arg HTTPS_PROXY=$HTTPS_PROXY'
            print(cmd)
            os.system(cmd)

    ### Any new sample creation must happened here.
    ### Each sample file requires the previous one to be created
    ### current order is : DepMap, Sanger, CPTAC, HCMI, BeatAML, MPNST
    if args.samples or args.all:
        ### build gene file
        cmd = docker_run+'depmap Rscript 00-buildGeneFile.R'
        print(cmd)
        os.system(cmd)
        
        ###build sample files
        cmd = docker_run+'depmap Rscript 01-depmapSamples.R'
        print(cmd)
       # os.system(cmd)

        cmd = docker_run+'depmap Rscript 01a-pullSamples_LINCS.R /tmp/depmap_samples.csv'
        print(cmd)
        os.system(cmd)
        
        cmd = docker_run+'cptac --geneFile=/tmp/genes.csv --prevSampleFile=/tmp/depmap_samples.csv'
        print(cmd)
        os.system(cmd)
        
        cmd =docker_run+'hcmi python 01-createHCMISamplesFile.py --samples=/tmp/cptac_samples.csv'
        print(cmd)
        os.system(cmd)
        
        cmd=docker_run+'beataml python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --samples=/tmp/hcmi_samples.csv'
        print(cmd)
        os.system(cmd)
        
        cmd = docker_run+'mpnst Rscript 00_sample_gen.R /tmp/beataml_samples.csv $SYNAPSE_AUTH_TOKEN'
        print(cmd)
        os.system(cmd)

    #### Any new omics files are created here.
    ### these are not order dependent but require gene and sample files
    if args.omics or args.all:
        ##omics data

        ###depmap cell line
        ocmd = docker_run+' depmap Rscript 02-pullDepMap.R /tmp/genes.csv /tmp/depmap_samples.csv'
        print(ocmd)
        os.system(ocmd)
        
        ocmd = docker_run+' depmap Rscript 02b-pullSanger.R /tmp/genes.csv /tmp/depmap_samples.csv'
        print(ocmd)
        os.system(ocmd)

        ocmd = docker_run+' depmap python3 depMapProts.py /tmp/genes.csv /tmp/depmap_samples.csv'
        print(ocmd)
        os.system(ocmd)
        
        ###cptac
        ocmd = docker_run+' cptac python3 getCptacData.py --geneFile=/tmp/genes.csv --curSampleFile=/tmp/cptac_samples.csv'
        print(ocmd)
        os.system(ocmd)
        
        ###HCMI - the folowing three steps are all required?
        ocmd = docker_run+' hcmi python 02-getHCMIData.py -m transcriptomics_gdc_manifest.txt -t transcriptomics -o hcmi_transcriptomics.csv'
        print(ocmd)
        os.system(ocmd)
        
        ocmd = docker_run+' hcmi /opt/venv/bin/python 02-getHCMData.py -m transcriptomics_gdc_manifest.txt -t transcriptomics -o hcmi_transcriptomics.csv'
        print(ocmd)
        os.system(ocmd)
        
        ocmd = docker_run+' hcmi /opt/venv/bin/python 02-getHCMData.py -m transcriptomics_gdc_manifest.txt -t transcriptomics -o hcmi_transcriptomics.csv'
        print(ocmd)
        os.system(ocmd)
        
        ##beataml
        ocmd=docker_run+' beataml python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --samples /tmp/beataml_samples.csv'
        print(ocmd)
        os.system(ocmd)
        

        ###mpnst
        ocmd=docker_run+' mpnst  Rscript 01_mpnst_get_omics.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/depmap/genes.csv'
        print(ocmd)
        os.system(ocmd)

    ### Drug matching scripts take a while
    ### they are their own step and also order dependent. Current order is:
    ### DepMap/Sanger, MPNST, LINCS
    if args.drugs or args.all:
        ###build drug data
        dcmd=docker_run+' depmap Rscript 03-createDrugFile.R CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM,NCI60'
        print(dcmd)
        os.system(dcmd)

        emd = docker_run+' mpnst Rscript  02_get_drug_data.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/depmap/drugs.tsv.gz'
        print(emd)
        os.system(emd)
        
        ##build experimental data
        ecmd= docker_run+' depmap /opt/venv/bin/python 04-cellLineDrugs_LINCS.py --drugFile /tmp/drugs.tsv.gz'
        print(ecmd)
        os.system(ecmd)
        
    if args.exp or args.all:
            
        dcmd=docker_run+' depmap /opt/venv/bin/python 04-drug_dosage_and_curves.py --drugfile=/tmp/drugs.tsv.gz --curSampleFile=/tmp/depmap_samples.csv'
        print(dcmd)
        os.system(dcmd)

        dcmd = docker_run+' mpnst Rscript 03_get_drug_response_data.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/depmap/drugs.tsv.gz'
        print(dcmd)
        os.system(dcmd)

        ecmd = docker_run+' cell-line Rscript 05-LINCS_perturbations.R /tmp/genes.csv /tmp/drugs.tsv.gz /tmp/depmap_samples.csv'




main()
