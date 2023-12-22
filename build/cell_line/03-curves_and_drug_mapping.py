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

####step 3a - get dose response data
cmd = 'Rscript 03a-drugAndResponseData.R'
os.system(cmd)

########Step 3b fit curves
allfiles=[a for a in os.listdir('./') if 'DoseResponse' in a]
print(allfiles)
for a in allfiles:
    os.system('python fit_curve.py --input '+a+' --output '+a)

###step 3c concatenate all files
os.system('cat *.0 > experiments_orig.tsv')

##now fix drug identifiers in experiments and drug files

### step 3d store files (figshare? ftp?
os.system('Rscript remapDrugsToSmiles.R drugs.tsv.gz experiments_orig.tsv')