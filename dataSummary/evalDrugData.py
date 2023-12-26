
'''
evalDrugData
Evaluates drug data and identifies gaps

'''

import matplotlib as mp
import pandas as pd

drugs = pd.read_csv('drugs.tsv.gz',sep='\t')

##which drugs are missing names
nonames = drugs[drugs['chem_name'].isna()]

##which drugs are missing smiles
nosmiles = drugs[drugs['canSMILES'].isna()].set_index('improve_chem_id')

##read in experiment data, figure out where those drugs came from
exp = pd.read_csv('experiments.tsv.gz',sep='\t')
drug_counts = exp.groupby('study')['improve_chem_id'].nunique()

exp = exp.set_index('improve_chem_id')
#join experimetns to drugs, figure out where missing came from


missing_exp = nosmiles.join(exp,on='improve_chem_id',how='left')

##count how many missing drugs per study
per_study = missing_exp.reset_index().groupby('study')['improve_chem_id'].nunique()

##count how many studies per missing drug
per_drug = missing_exp.groupby('improve_chem_id')['study'].nunique()
##now plot

combined = pd.DataFrame(drug_counts).rename(columns={'improve_chem_id':'total'}).join(pd.DataFrame(per_study)).rename(columns={'improve_chem_id':'missing'}).fillna(0.1)

combined = combined.reset_index().melt(value_vars=['total','missing'],id_vars='study').rename(columns={'variable':'Smiles','value':'Drugs'})

from plotnine import ggplot, aes, geom_bar, scales

p = ggplot(combined) + aes(x="study",fill='Smiles',y='Drugs') + geom_bar(stat='identity',position='dodge')+scales.scale_y_log10()
p.save(filename='drugNumbers.png')
