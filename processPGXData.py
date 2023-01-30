'''
This is a script that leverages PHarmacogx R package as well as the curve.py script
to collect all the data into a single file for CANDLE
_author_='Sara Gosline'
_email_ = 'sara.gosline@pnnl.gov'
'''

import pandas as pd
import os

datasets = ['prism','tavor','fimm','beatAML','UHNBreast','gbm','gray','pdtx']

datasets = ['prism','beatAML']
#empty list of files
dose_files = []
gex_files = []

'''
we iterate through every dataset
'''
for dat in datasets:
  print("Analyzing "+dat+' dataset')
  
  rcmd = 'Rscript pgx/'+dat+'dataFormat.R'
  
  dose_rep = 'python ../curve/fit_curve.py --input '+dat+'doseResponse.tsv --output '+dat+'drugSummaryStats.tsv'
  gex = dat+'geneExpression.tsv'

  ##run the commands
  os.system(rcmd)
  os.system(dose_rep)
  
  ##add the files to the list
  dose_files = dose_files.append(dat+'drugSummaryStats.tsv')
  gex_files = gex_files.append(dat+'geneExpression.tsv')


print('No we can append all the files into one')
