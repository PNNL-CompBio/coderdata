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
import glob
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

# A failure here is FATAL, not a warning.
#
# In build v32 this exited 1 after DOSERESP.zip truncated, run() merely logged
# "WARNING: command exited 1", and the build carried on. nci60DoseResponse was
# never written, so the fitting below saw 7 studies instead of 8 and nci60 was
# published as eight empty files -- which pass schema validation, because a
# file with no rows cannot violate a schema.
_rc = run(['/opt/venv/bin/python','04b-nci60-updated.py','--sampleFile='+samplefile,'--drugFile='+drugfile])
if _rc != 0:
    raise RuntimeError(
        f'04b-nci60-updated.py failed (exit {_rc}). It produces the NCI60 '
        f'dose-response input, and continuing without it publishes nci60 as '
        f'empty files that still validate.')

####step 4a - get dose response data
run(['Rscript','04a-drugResponseData.R',samplefile,drugfile,'CTRPv2,FIMM,GDSC'])
run(['Rscript','04a-drugResponseData.R',samplefile,drugfile,'gCSI,PRISM,CCLE'])


#cmd = 'Rscript 04a-drugResponseData.R '+samplefile+' '+drugfile+' NCI60'
#print(cmd)
#os.system(cmd)

########Step 4b fit curves
# Collect dose-response inputs from BOTH places they are written.
#
# 04b-nci60-updated.py writes ./nci60DoseResponse into the working directory,
# but 04a-drugResponseData.R writes through dose_response_cache_path(), i.e.
# $DOSE_RESPONSE_CACHE_DIR/<lowercase study>_doseResponse.tsv, so that a failed
# attempt can resume instead of regenerating every study. This scan only looked
# in './' for 'DoseResponse', which matches neither the directory nor the
# lowercase 'd', so it silently found NCI60 alone.
#
# The consequence was severe and completely silent: broad_sanger_experiments.tsv
# came out NCI60-only, and 05b_separate_datasets.py derives each split dataset's
# sample and drug ids FROM the experiments rows, so ccle, ctrpv2, fimm, gcsi,
# gdscv1, gdscv2 and prism were published with every file present but empty --
# samples, omics, drugs and experiments alike. Header-only files are schema
# valid, so validate passed. Build v27 shipped exactly that.
#
# Excluding '.0' matters too: fit_curve writes '<input>.0', which itself
# contains 'DoseResponse', so a retry would re-fit its own output.
CACHE_DIR = os.environ.get('DOSE_RESPONSE_CACHE_DIR', '/tmp/dose_response_cache')
jobs = []   # (input path, output stem -- kept local so '<stem>.0' lands in ./)
for a in sorted(os.listdir('./')):
    if 'DoseResponse' in a and not a.endswith('.0'):
        jobs.append((a, a))
for pth in sorted(glob.glob(os.path.join(CACHE_DIR, '*_doseResponse.tsv'))):
    stem = os.path.basename(pth)[:-len('.tsv')]
    if any(stem == j[1] for j in jobs):
        continue
    jobs.append((pth, stem))

if not jobs:
    raise RuntimeError('No DoseResponse files found — data generation step likely failed')
_log(f'fitting {len(jobs)} DoseResponse files with {opts.workers} workers: '
     f'{[j[1] for j in jobs]}')
fit_failures = []
for src, stem in jobs:
    _log(f'starting fit_curve: {stem} (from {src})')
    cmd = ['/opt/venv/bin/python', 'fit_curve.py',
           '--input=' + src, '--output=' + stem, '--workers=' + str(opts.workers)]
    if 'nci60' in stem.lower():
        cmd.append('--chunk_size=25000')
    rc = run(cmd)
    if rc != 0:
        _log(f'WARNING: fit_curve failed for {stem} (exit {rc}) — will be excluded from output')
        fit_failures.append(stem)
    _log(f'finished fit_curve: {stem}')

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

combined = pd.concat(final_file).drop_duplicates()

# Every study that was fitted must appear in the output.
#
# 05b_separate_datasets.py derives each split dataset's sample and drug ids from
# these rows, so a study missing here is not a partial result -- it publishes
# that whole dataset as header-only files: samples, omics, drugs, experiments.
# Those pass schema validation (no rows, nothing to violate), which is why v27
# shipped seven empty datasets and still reported "validate successful".
# Check against the studies this step is SUPPOSED to produce, not merely the
# ones that happened to yield an input file.
#
# The previous version derived the expectation from `jobs`, so a study whose
# dose-response file was never generated simply was not expected, and its
# absence passed silently. That is precisely how v32 shipped an empty nci60:
# the input never existed, so nothing noticed it was gone.
REQUIRED_STUDIES = {'nci60', 'ccle', 'ctrpv2', 'fimm', 'gcsi',
                    'gdscv1', 'gdscv2', 'prism'}
fitted_studies = {s.replace('_doseResponse', '').replace('DoseResponse', '').lower()
                  for s in (j[1] for j in jobs)} - {''}
not_fitted = sorted(REQUIRED_STUDIES - fitted_studies)
if not_fitted:
    raise RuntimeError(
        f"No dose-response input was produced for: {', '.join(not_fitted)}. "
        f"Fitted only {sorted(fitted_studies)}. Each missing study becomes an "
        f"entirely empty published dataset, so the build is stopping here.")
expected_studies = fitted_studies | REQUIRED_STUDIES
got_studies = {str(v).lower() for v in combined['study'].unique()} if 'study' in combined.columns else set()
missing_studies = sorted(st for st in expected_studies if st not in got_studies)
if missing_studies:
    raise RuntimeError(
        f"broad_sanger_experiments is missing every row for: {', '.join(missing_studies)}. "
        f"Present: {sorted(got_studies)}. Each missing study becomes an entirely empty "
        f"published dataset, so the build is stopping rather than shipping header-only files.")
_log(f'studies in output: {sorted(got_studies)} ({len(combined):,} rows)')

combined.to_csv('/tmp/broad_sanger_experiments.tsv',index=False,sep='\t')
_log(f'wrote /tmp/broad_sanger_experiments.tsv')
#os.system('cat *.0 > /tmp/broad_sanger_experiments.tsv')
#os.system('gzip -f /tmp/experiments.tsv')

