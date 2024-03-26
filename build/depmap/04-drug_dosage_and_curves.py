'''
This script represnets the end to end steps required to build the benchmark data module.
It requires R and some libraries.
TODO: add in R or dockerize requirements

4b- run fit_curve.py to fit doseRep.tsv files
5- concatenate doseRep files and compress
6- compress drugs.tsv
7- store data?

'''

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--curSampleFile',dest='samplefile',default=None,help='DepMap sample file')

parser.add_argument('--drugfile',dest='dfile',default=None,help='Drug database')

opts = parser.parse_args()

samplefile = opts.samplefile
drugfile = opts.dfile

####step 4a - get dose response data
cmd = 'Rscript 04a-drugResponseData.R '+samplefile+' '+drugfile+' CTRPv2,FIMM,GDSC'
print(cmd)
os.system(cmd)

cmd = 'Rscript 04a-drugResponseData.R '+samplefile+' '+drugfile+' gCSI,PRISM,CCLE'
print(cmd)
os.system(cmd)

cmd = 'Rscript 04a-drugResponseData.R '+samplefile+' '+drugfile+' NCI60'
print(cmd)
os.system(cmd)

########Step 4b fit curves
allfiles=[a for a in os.listdir('./') if 'DoseResponse' in a]
print(allfiles)
for a in allfiles:
    os.system('/opt/venv/bin/python fit_curve.py --input '+a+' --output '+a)

###step 4c concatenate all files

os.system('cat *.0 > /tmp/depmap_experiments.tsv')
#os.system('gzip -f /tmp/experiments.tsv')

