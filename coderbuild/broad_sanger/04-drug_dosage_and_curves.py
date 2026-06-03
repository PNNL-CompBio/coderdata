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
import multiprocessing
from datetime import datetime

def _log(msg):
    print(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}] {msg}', flush=True)

parser = argparse.ArgumentParser()
parser.add_argument('--curSampleFile',dest='samplefile',default=None,help='DepMap sample file')
parser.add_argument('--drugfile',dest='dfile',default=None,help='Drug database')
parser.add_argument('--workers', type=int, default=max(1, round(multiprocessing.cpu_count() * 3 / 4)), help='Number of parallel worker processes for curve fitting (default: 3/4 of available CPUs)')

opts = parser.parse_args()

samplefile = opts.samplefile
drugfile = opts.dfile

def run(cmd):
    _log('running: ' + ' '.join(cmd))
    res = subprocess.run(cmd)
    if res.returncode != 0:
        _log(f'WARNING: command exited {res.returncode}: {" ".join(cmd)}')
    return res.returncode

run(['/opt/venv/bin/python','04b-nci60-updated.py','--sampleFile='+samplefile,'--drugFile='+drugfile])

####step 4a - get dose response data
run(['Rscript','04a-drugResponseData.R',samplefile,drugfile,'CTRPv2,FIMM,GDSC'])
run(['Rscript','04a-drugResponseData.R',samplefile,drugfile,'gCSI,PRISM,CCLE'])


#cmd = 'Rscript 04a-drugResponseData.R '+samplefile+' '+drugfile+' NCI60'
#print(cmd)
#os.system(cmd)

########Step 4b fit curves
allfiles=[a for a in os.listdir('./') if 'DoseResponse' in a]
if not allfiles:
    raise RuntimeError('No DoseResponse files found — data generation step likely failed')
_log(f'fitting {len(allfiles)} DoseResponse files with {opts.workers} workers: {allfiles}')
fit_failures = []
for a in allfiles:
    _log(f'starting fit_curve: {a}')
    cmd = ['/opt/venv/bin/python', 'fit_curve.py',
           '--input=' + a, '--output=' + a, '--workers=' + str(opts.workers)]
    if 'nci60' in a.lower():
        cmd.append('--chunk_size=25000')
    rc = run(cmd)
    if rc != 0:
        _log(f'WARNING: fit_curve failed for {a} (exit {rc}) — will be excluded from output')
        fit_failures.append(a)
    _log(f'finished fit_curve: {a}')

###step 4c concatenate all files
outfiles = [a for a in os.listdir("./") if a.endswith('.0')]
if not outfiles:
    raise RuntimeError('No .0 output files found after curve fitting — all fit_curve steps failed')
final_file = []
for of in outfiles:
    df = pd.read_csv(of, sep='\t')
    if df.empty:
        _log(f'WARNING: {of} is empty, skipping')
        continue
    _log(f'loaded {of}: {len(df)} rows')
    final_file.append(df)

if not final_file:
    raise RuntimeError('All .0 output files were empty — cannot produce experiments output')
if fit_failures:
    _log(f'WARNING: {len(fit_failures)} dataset(s) missing from output: {fit_failures}')

pd.concat(final_file).drop_duplicates().to_csv('/tmp/broad_sanger_experiments.tsv',index=False,sep='\t')
_log(f'wrote /tmp/broad_sanger_experiments.tsv')
#os.system('cat *.0 > /tmp/broad_sanger_experiments.tsv')
#os.system('gzip -f /tmp/experiments.tsv')

