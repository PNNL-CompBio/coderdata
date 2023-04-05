'''
This script represnets the end to end steps required to build the benchmark data module.
It requires R and some libraries.
TODO: add in R or dockerize requirements

1- load up Gene names, create data/genes.csv
2- load up Sample names, create data/samples.csv
3- use pharmacoGX to build expression.csv based on 1 and 2
4a- use pharmacoGX to build doseRep.tsv files (too big for 1) and drugs.tsv and 2
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
os.system('Rscript -e "setwd(".."); source("pgx/pgxToImprove.R"); getCellLineExpData()') ##makes expression csv
##move it to dta andnow zip it up


#########Step 4a get dose response data from PGX
os.system('')
os.system('')


########Step 4b fit curves
