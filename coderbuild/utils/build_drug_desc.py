'''
build drug descriptor table from drug table


'''


import argparse
import os
import tempfile

# ---------------------------------------------------------------------------
# Keep multiprocessing's scratch files OFF the bind mount.
#
# This MUST run before multiprocessing creates any worker, so it sits above the
# imports that pull it in.
#
# build_all.py bind-mounts the host's local/ directory at /tmp inside every
# container. Python 3.14 changed the default multiprocessing start method on
# Linux from "fork" to "forkserver", and forkserver opens a Unix domain socket
# under tempfile.gettempdir() -- i.e. /tmp, i.e. the bind mount. Docker's macOS
# file sharing does not support chmod() on a socket there, so worker startup
# died with:
#
#   OSError: [Errno 22] Invalid argument: '/tmp/pymp-8ai8jp5x/sock-2c162e8922a7'
#
# Mordred computes descriptors in parallel, so this killed the drugs step --
# and it did so at the very END of build_drugs.sh, after nci60, the PSets and
# the join had all completed, wasting the whole multi-hour drug build on every
# attempt. It is deterministic, so retries do not help.
#
# Pointing TMPDIR at container-local storage fixes it at the source. Verified:
# with /tmp as a bind mount a Pool raises EINVAL; with this set it succeeds.
# Moving the socket is preferred over forcing start_method="fork", which would
# reintroduce the fork-with-threads hazard that motivated the upstream change.
_MP_TMPDIR = os.environ.get("CODERDATA_MP_TMPDIR", "/opt/mp_tmp")
os.makedirs(_MP_TMPDIR, exist_ok=True)
os.environ["TMPDIR"] = _MP_TMPDIR
tempfile.tempdir = _MP_TMPDIR
# ---------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import ConvertToNumpyArray
import pandas as pd
import numpy as np
from mordred import Calculator, descriptors
import multiprocessing
import gzip

# Remove all of the Deprecation warnings. There is a github issue to update the code and 50k lines of warnings in the build log so I'm hiding them for now. 
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
# If this script suddently stops working, it is likely due to a change in rdkit. Unhide the warnings to see the error.


def smiles_to_fingerprint(smiles):
    '''
    Takes all SMILES and create morgan fingerprints for them
    '''
    fdict = []
    ##get morgan fingerprint
    print('Computing morgan fingerprints for '+str(len(smiles))+' SMILES')
 #   morgan_fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024, useCountSimulation=False)
    for s in smiles:
       # print(s)
        mol = Chem.MolFromSmiles(s)
        try:
            #this has been depracated despite being in Alex's original script
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)  # update these parameters
  #          fingerprint = morgan_fp_gen.GetFingerprint(mol)
            #            vec2 = np.array(fp2)
        except:
            print('Cannot compute fingerprint for '+s)
            continue
        fingerprint_array = np.array(fingerprint)
        fstr = ''.join([str(a) for a in fingerprint_array])
        fdict.append({'smile':s,'descriptor_value':fstr,'structural_descriptor':'morgan fingerprint'})
        
    return pd.DataFrame(fdict)#fingerprint_array


def smiles_to_mordred(smiles,nproc=2):
    '''
    get descriptors - which ones?
    '''
    print('Computing mordred descriptors for '+str(len(smiles))+' SMILES')

    
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    smols = []
    ssmil = []
    for i in range(len(mols)):
        m = mols[i]
        if m is not None:
            smols.append(m)
            ssmil.append(smiles[i])

    calc = Calculator(descriptors, ignore_3D=True)
    dd = calc.pandas(mols=smols, nproc=nproc, quiet=False, ipynb=False )
    values = dd.columns
    dd['smile'] = ssmil
    ##reformat here
    longtab = pd.melt(dd,id_vars='smile',value_vars=values)
    longtab = longtab.rename({'variable':'structural_descriptor','value':'descriptor_value'},axis=1)
    
    return longtab

## Number of SMILES whose descriptors are held in memory at once.
##
## Mordred returns a WIDE frame of ~1,600 object-dtype descriptors per molecule,
## which pd.melt() then explodes into rows: at 56,890 SMILES that is roughly 92
## MILLION rows, before the merge, the concat and two more full-frame copies in
## the cleaning step. Build v23 was OOM-killed doing exactly that
## ("build_drugs.sh: line 21: Killed") on a 23.4GiB Docker VM.
##
## v19 completed the same step at 46,319 SMILES, so this was never comfortable
## -- a 23% growth in the drug table was enough to cross the limit. Chunking
## makes peak memory independent of how many drugs the release contains, so it
## does not silently creep back up to the ceiling next time.
SMILES_CHUNK = int(os.environ.get("DRUG_DESC_SMILES_CHUNK", "2000"))


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def _clean_block(block):
    """Coerce bad values and drop malformed ids for one block of rows."""
    # Mordred puts the exception object into a cell when a descriptor cannot be
    # computed for a molecule, so str() of the error lands in the table (e.g. the
    # "module 'numpy' has no attribute 'float'" AttributeError from np.float being
    # removed in NumPy >= 1.24, or messages containing "invalid"/"missing"/"min...").
    # Every Mordred descriptor is numeric, so coerce any non-numeric value to
    # "NaN" -- this robustly catches ANY such error string, unlike the previous
    # fixed substring blocklist which missed the numpy AttributeError. Morgan
    # fingerprints are legitimate bit strings and are left untouched.
    block['descriptor_value'] = block['descriptor_value'].astype(str)
    is_fingerprint = block['structural_descriptor'] == 'morgan fingerprint'
    numeric = pd.to_numeric(block['descriptor_value'], errors='coerce')

    # Render every Mordred value through float, so the text written does not
    # depend on where the chunk boundaries fall. Mordred returns ints for some
    # descriptors; whether a given column also picks up a NaN decides if pandas
    # holds it as int64 ("0") or float64 ("0.0"). That makes the SAME descriptor
    # serialise differently depending on how many molecules share a chunk, which
    # would make the output non-deterministic. Fingerprints are bit strings and
    # are left untouched.
    block.loc[~is_fingerprint, 'descriptor_value'] = (
        numeric[~is_fingerprint]
        .map(lambda v: "NaN" if pd.isna(v) else repr(float(v)))
    )

    # Remove Data that is incorrectly written by mordred or rdkit. - Very rare bug, but it happens.
    block['improve_drug_id'] = block['improve_drug_id'].astype(str).str.strip()
    return block[block['improve_drug_id'].str.match(r'^SMI_\d+$')]


def main():
    parser = argparse.ArgumentParser('Build drug descriptor table')
    parser.add_argument('--drugtable',dest='drugtable')
    parser.add_argument('--desctable',dest='outtable')

    args = parser.parse_args()

    cores = multiprocessing.cpu_count()
    ncors = cores-1
    tab = pd.read_csv(args.drugtable,sep='\t')

    cansmiles = [a for a in set(tab.canSMILES) if str(a)!='nan']

    ids = pd.DataFrame(tab[['improve_drug_id','canSMILES']]).drop_duplicates()
    ids = ids.rename({"canSMILES":'smile'},axis=1)
    del tab

    cols = ['improve_drug_id','structural_descriptor','descriptor_value']
    total = len(cansmiles)
    print(f'Building descriptors for {total} SMILES in chunks of {SMILES_CHUNK}',
          flush=True)

    # Stream each block straight to the gzip output. Morgan fingerprints are
    # emitted first, then Mordred descriptors, preserving the original row order.
    with gzip.open(args.outtable, 'wt', newline='') as out:
        header = True
        for label, fn in (('morgan', smiles_to_fingerprint),
                          ('mordred', lambda c: smiles_to_mordred(c, nproc=ncors))):
            for i, chunk in enumerate(_chunks(cansmiles, SMILES_CHUNK), start=1):
                desc = fn(chunk)
                if desc is None or desc.empty:
                    continue
                block = ids.merge(desc)[cols]
                block = _clean_block(block)
                block.to_csv(out, sep='\t', index=False, header=header)
                header = False
                done = min(i * SMILES_CHUNK, total)
                print(f'  [{label}] {done}/{total} SMILES -> {len(block)} rows',
                      flush=True)
                del desc, block

if __name__=='__main__':
    main()
