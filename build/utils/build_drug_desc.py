'''
build drug descriptor table from drug table


'''


import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from mordred import Calculator, descriptors

def smiles_to_fingerprint(smiles):
    '''
    takes smiles nad create morgan fingerprint
    '''
    fdict = []
    ##get morgan fingerprint
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)  # update these parameters
        fingerprint_array = np.array(fingerprint)
        fdict.append({'smile':s,'descriptor_value':str(fingerprint_array),'structural_descriptor':'morgan fingerprint'})
        
    return pd.DataFrame(fdict)#fingerprint_array


def smiles_to_mordred(smiles):
    '''
    get descriptors - which ones?
    '''
    

def main:
    parser = argparse.ArgumentParser('Build drug descriptor table')
    parser.add_arg('--drugtable',dest='drugtable')
    parser.add_arg('--desctable',dest='outtable')

    args  = parser.parse_args()

    tab = pd.read_csv(args.drugtable,sep='\t')

    cansmiles = list(set(tab.canSMILES))
    #    isosmiles = list(set(tab.isoSMILES))
    morgs = smiles_to_fingerprint(cansmiles)

    ids = pd.DataFrame(tab[['improve_drug_id','canSMILES']]).drop_duplicates()
    ids = ids.rename({"canSMILES":'smile'},axis=1).merge(morgs)[['improve_drug_id','structural_descriptor','descriptor_value']]

    ids.to_csv(args.outtable,sep='\t',index=False)
