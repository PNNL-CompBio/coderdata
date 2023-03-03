'''
This script represnets the end to end steps required to build the benchmark data module.
It requires R and some libraries.
TODO: add in R or dockerize requirements

1- load up Gene names, create data/genes.csv
2- load up Sample names, create data/samples.csv
3- use pharmacoGX to build doseRep.tsv files (too big for 1) and drugs.tsv and expression.csv files
4- run fit_curve.py to fit doseRep.tsv files
5- concatenate doseRep files and compress
6- compress drugs.tsv
7- store data?

'''


import os


###########Step 1, build genes
os.system('Rscript data/')


###########Step 2, build sample names
