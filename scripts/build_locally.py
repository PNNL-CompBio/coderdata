"""
script that builds the coderdata package and stores locally
"""


import os
if not os.path.exists('local'):
    os.mkdir('local')
docker_run='docker run -v $PWD/local:/tmp '
###build docker images
###
dlist = []
for fn in os.listdir("build/docker"):
    dsname=fn.split('.')[1]
    print(dsname)
    dlist.append(dsname)
    cmd = 'docker build -t '+dsname+' .  -f build/docker/'+fn
    print(cmd)
    os.system(cmd)

### build gene file
cmd = docker_run+'cell_line Rscript 00-buildGeneFile.R'
print(cmd)
os.system(cmd)

###build sample files
cmd = docker_run+'/ cell_line Rscript 01-cellLineSamples.R'
print(cmd)
os.system(cmd)

cmd = docker_run+' cptac python3 getCptacData.py --geneFile=/tmp/genes.csv --prevSampleFile=/tmp/cell_line_samples.csv'
print(cmd)
os.system(cmd)

cmd =docker_run+' hcmi python 01-createHCMISamplesFile.py'
print(cmd)
os.system(cmd)

cmd=docker_run+' beataml python GetBeatAML.py --token $SYNAPSE_AUTH_TOKEN --samples=/tmp/hcmi_samples.csv'
print(cmd)
os.system(cmd)

cmd = docker_run+' mpnst Rscript 00_sample_gen.R /tmp/beataml_samples.csv $SYNAPSE_AUTH_TOKEN'
print(cmd)
os.system(cmd)

##omics data

###depmap cell line
ocmd = docker_run+'/ cell_line Rscript 02-cellLineDepMap.R /tmp/genes.csv /tmp/cell_line_samples.csv'
print(ocmd)
os.system(ocmd)

ocmd = docker_run+'/ cell_line Rscript 02b-cellLineSanger.R /tmp/genes.csv /tmp/cell_line_samples.csv'
print(ocmd)
os.system(ocmd)

###cptac
ocmd = docker_run+' cptac python3 getCptacData.py --geneFile=/tmp/genes.csv --curSampleFile=/tmp/cptac_samples.csv'
print(ocmd)
os.system(ocmd)

###HCMI
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
ocmd=docker_run+' mpnst  Rscript 01_mpnst_get_omics.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/cell_line/genes.csv'
print(ocmd)
os.system(ocmd)

###build drug data
dcmd=docker_run+' cell_line Rscript 03-createDrugFile.R CTRPv2,GDSC,gCSI,PRISM,CCLE,FIMM,NCI60'
print(dcmd)
os.system(dcmd)

dcmd=docker_run+' cell_line /opt/venv/bin/python 04-drug_dosage_and_curves.py --drugfile=/tmp/drugs.tsv.gz --curSampleFile=/tmp/cell_line_samples.csv'
print(dcmd)
os.system(dcmd)

##build experimental data
ecmd= docker_run+' cell-line /opt/venv/bin/python 04-cellLineDrugs_LINCS.py --drugFile /tmp/drugs.tsv.gz'
ecmd = docker_run+' cell-line Rscript 05-LINCS_perturbations.R /tmp/genes.csv /tmp/drugs.tsv.gz /tmp/cell_line_samples.csv'
emd = docker_run+' mpnst Rscript  02_get_drug_data.R $SYNAPSE_AUTH_TOKEN /tmp/MPNST_samples.csv /tmp/cell_line/drugs.tsv.gz' 





