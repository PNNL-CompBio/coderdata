'''
This script represnets the end to end steps required to build the benchmark data module.
It requires R and some libraries.
TODO: add in R or dockerize requirements

1- load up Gene names, create data/genes.csv
2- load up Sample names, create data/samples.csv
[3- use pharmacoGX to build expression.csv based on 1 and 2] - skipping now due to Priya's data
4a- use pharmacoGX to build *doseResponse files 
4b- run fit_curve.py to fit doseRep.tsv files
5- concatenate doseRep files and compress
6- compress drugs.tsv
7- store data?

'''


import os


###########Step 1, build genes
os.system('Rscript initialGeneDB.R') ###requires data/genes.csv to exist for step 3


###########Step 2, build sample names
os.system('Rscript initialSampleDB.R') ###requires data/samples.csv to exist for step 3

#########Step 3, get expression data
#os.system('Rscript -e "setwd(".."); source("pgx/pgxToImprove.R"); getCellLineExpData()') ##makes expression csv
os.system('Rscript getLatestCCLE.R') ###gets data from Priya's formatted data
##move it to dta andnow zip it up


#########Step 4a get dose response data from PGX
os.system('Rscript ../pgx/pgxToImprove.R')

########Step 4b fit curves
allfiles=[a for a in os.listdir('./') if 'DoseResponse' in a]
print(allfiles)
for a in allfiles:
    os.system('python fit_curve.py --input '+a+' --output '+a)

###step 5 concatenate all files
os.system('cat *.0 > experiments.tsv')

##step 6 zip up all files, upload them to fighsare

##now fix drug identifiers


### step 7 store files (figshare? ftp?
