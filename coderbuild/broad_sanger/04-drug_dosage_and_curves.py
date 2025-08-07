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
import pandas as pd
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--curSampleFile',dest='samplefile',default=None,help='DepMap sample file')

parser.add_argument('--drugfile',dest='dfile',default=None,help='Drug database')

opts = parser.parse_args()

samplefile = opts.samplefile
drugfile = opts.dfile

cmd = ['/opt/venv/bin/python','04b-nci60-updated.py','--sampleFile='+samplefile,'--drugFile='+drugfile]
print(cmd)
subprocess.run(cmd)

####step 4a - get dose response data
cmd = ['Rscript','04a-drugResponseData.R',samplefile,drugfile,'CTRPv2,FIMM,GDSC']
print(cmd)
subprocess.run(cmd)

cmd = ['Rscript','04a-drugResponseData.R',samplefile,drugfile,'gCSI,PRISM,CCLE']
print(cmd)
subprocess.run(cmd)


#cmd = 'Rscript 04a-drugResponseData.R '+samplefile+' '+drugfile+' NCI60'
#print(cmd)
#os.system(cmd)

########Step 4b fit curves
allfiles=[a for a in os.listdir('./') if 'DoseResponse' in a]
print(allfiles)
for a in allfiles:
    subprocess.run(['/opt/venv/bin/python','fit_curve.py','--input='+a,'--output='+a])

###step 4c concatenate all files
outfiles = [a for a in os.listdir("./") if ".0" in a]
final_file = []
for of in outfiles:
    final_file.append(pd.read_csv(of,sep='\t'))

pd.concat(final_file).to_csv('/tmp/broad_sanger_experiments.tsv',index=False,sep='\t')
#os.system('cat *.0 > /tmp/broad_sanger_experiments.tsv')
#os.system('gzip -f /tmp/experiments.tsv')

