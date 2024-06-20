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
    print('Computing morgan fingerprints for '+str(len(smiles))+' SMILES')
    for s in smiles:
       # print(s)
        mol = Chem.MolFromSmiles(s)
        try:
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)  # update these parameters
        except Error:
            print('Cannot compute fingerprint for '+s)
            continue
        fingerprint_array = np.array(fingerprint)
        fstr = ''.join([str(a) for a in fingerprint_array])
        fdict.append({'smile':s,'descriptor_value':fstr,'structural_descriptor':'morgan fingerprint'})
        
    return pd.DataFrame(fdict)#fingerprint_array


def smiles_to_mordred(smiles):
    '''
    get descriptors - which ones?
    '''
    print('Computing mordred descriptors for '+str(len(smiles))+' SMILES')

    mols = [Chem.MolFromSmiles(s) for s in smiles]

    calc = Calculator(descriptors, ignore_3D=True)
    dd = calc.pandas( mols, nmols=None, quiet=False, ipynb=False )
    values = dd.columns
    dd['smile'] = smiles
    ##reformat here
    longtab = pd.melt(dd,id_vars='smile',value_vars=values)
    longtab = longtab.rename({'variable':'structural_descriptor','value':'descriptor_value'},axis=1)
    return longtab

def main():
    parser = argparse.ArgumentParser('Build drug descriptor table')
    parser.add_argument('--drugtable',dest='drugtable')
    parser.add_argument('--desctable',dest='outtable')

    args  = parser.parse_args()

    print('Adding drug table for '+args.drugtable)
    tab = pd.read_csv(args.drugtable,sep='\t')

    cansmiles = [a for a in set(tab.canSMILES) if str(a)!='nan']
    #    isosmiles = list(set(tab.isoSMILES))
    morgs = smiles_to_fingerprint(cansmiles)

    ids = pd.DataFrame(tab[['improve_drug_id','canSMILES']]).drop_duplicates()
    id_morg = ids.rename({"canSMILES":'smile'},axis=1).merge(morgs)[['improve_drug_id','structural_descriptor','descriptor_value']]

    mords = smiles_to_mordred(cansmiles)
    
    id_mord = ids.rename({'canSMILES':'smile'},axis=1).merge(mords)[['improve_drug_id','structural_descriptor','descriptor_value']]
    
    full = pd.concat([id_morg,id_mord],axis=0)                     
    full.to_csv(args.outtable,sep='\t',index=False,compression='gzip')

if __name__=='__main__':
    main()
