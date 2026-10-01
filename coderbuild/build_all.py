"""
script that builds the coderdata package and stores locally
"""

import os
import argparse
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor
import shutil
import gzip
from glob import glob
import sys
import requests
import threading
import atexit
import signal
import tempfile
import yaml
from datetime import datetime

def _log(*args, **kwargs):
    ts = datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    print(ts, *args, **kwargs)

def _redact(cmd):
    """Render a command for logging with secret values masked.

    Tokens are passed to containers as `-e NAME=value`, so logging the command
    verbatim wrote them in cleartext. build_all_v19.log contained the live
    SYNAPSE_AUTH_TOKEN 47 times, and the build logs sit untracked in the repo
    root where a `git add -A` would commit them.
    """
    out = []
    for part in cmd:
        if '=' in part:
            key, _, val = part.partition('=')
            if val and any(m in key.upper() for m in ('TOKEN', 'SECRET', 'PASSWORD', 'KEY')):
                part = f'{key}=***REDACTED***'
        out.append(part)
    return ' '.join(out)

def main():
    parser=argparse.ArgumentParser(
        description="This script initializes all docker containers, builds datasets, validates them, and uploads to Figshare.",
        epilog="""Examples of usage:

Build all datasets in a high memory environment, validate them, and upload to Figshare:
  python coderbuild/build_all.py --all --high_mem --validate --figshare --version 0.1.29

Build only experiment files. This assumes preceding steps (docker images, samples, omics, and drugs) have already been completed:
  python coderbuild/build_all.py --exp

Validate all local files without building or uploading. These files must be located in ./local. Includes compression/decompression steps.
  python coderbuild/build_all.py --validate

Upload the latest data to Figshare (ensure tokens are set in the local environment):
  python coderbuild/build_all.py --figshare --version 0.1.30
        """
    )
    parser.add_argument('--docker',dest='docker',default=False,action='store_true', help="Build all docker images.")
    parser.add_argument('--samples',dest='samples',default=False,action='store_true', help="Build all sample files.")
    parser.add_argument('--omics',dest='omics',default=False,action='store_true', help="Build all omics files.")
    parser.add_argument('--drugs',dest='drugs',default=False,action='store_true', help="Build all drug files")
    parser.add_argument('--prev_drugs', type=str, default='', help='Comma-separated list of previously built drug files (e.g. local_old/broad_sanger_drugs.tsv,local_old/beataml_drugs.tsv). Each matching dataset is still rebuilt, but its previous file seeds pubchem_retrieval so only genuinely new drugs are queried. Existing drugs are preserved; removed drugs fall off naturally.')
    parser.add_argument('--misc', action='store_true', help="Run the final misc post-build step (e.g., split broad_sanger datasets).")
    parser.add_argument('--exp',dest='exp',default=False,action='store_true', help="Build all experiment file.")
    parser.add_argument('--validate', action='store_true', help="Run schema checker on all local files. Note this will be run, whether specified or not, if figshare arguments are included.")
    parser.add_argument('--figshare', action='store_true', help="Upload all local data to Figshare. FIGSHARE_TOKEN must be set in local environment.")
    parser.add_argument('--all',dest='all',default=False,action='store_true', help="Run all data build commands. This includes docker, samples, omics, drugs, exp arguments. This does not run the validate or figshare commands")
    parser.add_argument('--high_mem',dest='high_mem',default=False,action='store_true',help = "If you have 32 or more CPUs, this option is recommended. It will run many code portions in parallel. If you don't have enough memory, this will cause a run failure.")
    parser.add_argument('--dataset',dest='datasets',default='broad_sanger,beataml,pancreatic,bladder,sarcoma,liver,novartis,colorectal,mpnst,cptac,hcmi',help='Datasets to process. Defaults to all available.')
    parser.add_argument('--version', type=str, required=False, help='Version number for the Figshare upload title (e.g., "0.1.29"). This is required for Figshare upload. This must be a higher version than previously published versions.')
    parser.add_argument('--github-username', type=str, required=False, help='GitHub username for the repository.')
    parser.add_argument('--github-email', type=str, required=False, help='GitHub email for the repository.')
    parser.add_argument('--depmap-ready', dest='depmap_ready', default=False, action='store_true',
                        help='Confirm that the static DepMap files on Synapse are up to date for '
                             'the current DepMap release. Required for any build that includes the '
                             'broad_sanger dataset, because DepMap files can no longer be downloaded '
                             'programmatically and must be refreshed by hand.')
    parser.add_argument('--continue', dest='resume', default=False, action='store_true',
                        help='Resume a partial build: skip steps whose sentinel output already '
                             'exists in local/, rebuild Docker images only for datasets that '
                             'failed in the previous run, and log all activity to '
                             'local/build_progress.log.')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # Resume / progress tracking helpers
    # -------------------------------------------------------------------------
    PROGRESS_LOG = 'local/build_progress.log'
    _progress_lock = threading.Lock()

    def _step_sentinel(step_name):
        """Return the file/dir whose existence proves this step completed."""
        if step_name == 'genes file':
            return 'local/genes.csv'
        if step_name == 'phosphosites file':
            return 'local/phosphosites.csv'
        # All other steps are "{dataset} {type}" e.g. "colorectal omics"
        da, _, stype = step_name.rpartition(' ')
        if not stype:
            return None
        if stype == 'samples':
            return f'local/{da}_samples.csv'
        if stype == 'drugs':
            return f'local/{da}_drugs.tsv'
        if stype == 'omics':
            # transcriptomics is the last large file written — proves full completion
            for ext in ('.csv.gz', '.csv'):
                candidate = f'local/{da}_transcriptomics{ext}'
                if os.path.exists(candidate):
                    return candidate
            # default expectation for first-time check
            gzip_ds = {'broad_sanger', 'beataml', 'pancreatic', 'hcmi'}
            ext = '.csv.gz' if da in gzip_ds else '.csv'
            return f'local/{da}_transcriptomics{ext}'
        if stype == 'experiments':
            return f'local/{da}_experiments.tsv'
        if stype == 'misc':
            return 'local/all_files_dir'
        return None

    def _sentinel_done(step_name):
        """Return True if the sentinel for this step exists and is non-trivial."""
        s = _step_sentinel(step_name)
        if s is None:
            return False
        if os.path.isdir(s):
            return bool(os.listdir(s))
        return os.path.isfile(s) and os.path.getsize(s) > 100

    def _step_cleanup_globs(step_name):
        """Glob patterns (relative to cwd) to remove when a step fails."""
        if step_name == 'genes file':
            return ['local/genes.csv']
        if step_name == 'phosphosites file':
            return ['local/phosphosites.csv']
        da, _, stype = step_name.rpartition(' ')
        if not stype:
            return []
        if stype == 'samples':
            return [f'local/{da}_samples.csv']
        if stype == 'drugs':
            return [f'local/{da}_drugs.tsv', f'local/{da}_drug_descriptors.tsv.gz']
        if stype == 'omics':
            if da == 'broad_sanger':
                # broad_sanger omics produces both broad_* and sanger_* files
                return [
                    'local/broad_mutations*', 'local/broad_transcriptomics*',
                    'local/broad_copy_number*', 'local/broad_proteomics*',
                    'local/sanger_mutations*', 'local/sanger_transcriptomics*',
                    'local/sanger_copy_number*', 'local/sanger_proteomics*',
                ]
            return [
                f'local/{da}_mutations*', f'local/{da}_transcriptomics*',
                f'local/{da}_copy_number*', f'local/{da}_proteomics*',
                f'local/{da}_phosphoproteomics*',
            ]
        if stype == 'experiments':
            return [f'local/{da}_experiments.tsv']
        if stype == 'misc':
            return ['local/all_files_dir']
        return []

    def _cleanup_step(step_name):
        """Remove partial outputs from a failed step; log what was removed."""
        removed = []
        for pattern in _step_cleanup_globs(step_name):
            for f in glob(pattern):
                try:
                    if os.path.isdir(f):
                        shutil.rmtree(f)
                    else:
                        os.remove(f)
                    removed.append(f)
                except Exception as exc:
                    _log(f'[cleanup] could not remove {f}: {exc}')
        if removed:
            _log(f'[cleanup] {step_name}: removed {len(removed)} file(s): {", ".join(removed)}')
            _write_progress(step_name, 'CLEANUP', f'removed: {", ".join(removed)}')
        else:
            _log(f'[cleanup] {step_name}: nothing to remove')

    def _clean_local_intermediates(reason='', protect_prefixes=()):
        """Move leftover intermediate files in local/ into a quarantine
        directory (local/unexpected_intermediates/) to keep the working set
        small (raw downloads like rnaseq_all_data_*.csv, cell.xml,
        proteomics_all_*.csv, *.xlsx, etc. are never needed after the step that
        consumed them and otherwise accumulate until a later step hits ENOSPC).

        Files are MOVED, not deleted, so a mis-declared build output is never
        silently lost -- if something you expected to be uploaded ends up in
        local/unexpected_intermediates/, it means it was not declared in
        schema/expected_files.yaml. NOTE: because the quarantine directory is on
        the same filesystem, moving does NOT reclaim disk space during the
        build; clear local/unexpected_intermediates/ (or point it at another
        volume) if you hit disk pressure.

        Safety: a file is KEPT in place if it is a schema-declared output (per
        schema/expected_files.yaml -- these may still be sitting in local/ before
        the move to all_files_dir), a reference file, an incremental-build cache,
        build infrastructure, or matches one of protect_prefixes (e.g.
        broad*/sanger*, which build_misc.sh still needs to split broad_sanger).
        Only regular files matching local/*.* are considered (mirrors the move
        step's own glob), so directories like all_files_dir are never touched.
        """
        try:
            with open('schema/expected_files.yaml') as _f:
                _schema = yaml.safe_load(_f)
        except Exception as exc:
            _log(f'[cleanup-intermediates] {reason}: skipped, cannot read schema: {exc}')
            return
        keep = set()
        for _entries in _schema.get('datasets', {}).values():
            for _entry in _entries:
                _bn = os.path.basename(_entry['file'])
                keep.add(_bn)
                keep.add(_bn + '.gz')
                if _bn.endswith('.gz'):
                    keep.add(_bn[:-3])
        protected_exact = {'genes.csv', 'phosphosites.csv', 'prev_drug_files.txt'}
        quarantine = os.path.join('local', 'unexpected_intermediates')
        moved, moved_bytes = 0, 0
        for f in glob(os.path.join('local', '*.*')):
            if os.path.isdir(f):
                continue
            bn = os.path.basename(f)
            if bn in keep or bn in protected_exact:
                continue
            if '_prev.' in bn or bn.endswith(('.cid', '.log', '.jsonl')):
                continue
            if any(bn.startswith(p) for p in protect_prefixes):
                continue
            try:
                os.makedirs(quarantine, exist_ok=True)
                sz = os.path.getsize(f)
                dest = os.path.join(quarantine, bn)
                if os.path.exists(dest):
                    os.remove(dest)  # replace a stale copy from an earlier pass
                shutil.move(f, dest)
                moved += 1
                moved_bytes += sz
            except Exception as exc:
                _log(f'[cleanup-intermediates] could not move {f}: {exc}')
        if moved:
            _log(f'[cleanup-intermediates] {reason}: moved {moved} file(s) '
                 f'({moved_bytes / 1e9:.2f} GB) to {quarantine}/')
        else:
            _log(f'[cleanup-intermediates] {reason}: nothing to move')

    def _write_progress(step, status, note=''):
        """Append one line to local/build_progress.log (thread-safe)."""
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        line = f'{ts} | {status:<8} | {step}'
        if note:
            line += f'  ({note})'
        line += '\n'
        with _progress_lock:
            try:
                with open(PROGRESS_LOG, 'a') as f:
                    f.write(line)
            except Exception:
                pass  # never let logging crash the build

    def _step_to_docker_dataset(step_name):
        """Map a step name to the dataset name used by process_docker."""
        if step_name in ('genes file', 'gene file'):
            return 'genes'
        if step_name == 'phosphosites file':
            return 'phosphosites'
        da, _, stype = step_name.rpartition(' ')
        if not stype:
            return None
        return da  # e.g. "colorectal", "broad_sanger", "liver"

    def _last_run_failures():
        """
        Parse build_progress.log and return the set of step names that
        FAILED in the most recent run (since the last START line) and were
        not subsequently fixed (no SUCCESS line after the FAILURE).
        """
        if not os.path.exists(PROGRESS_LOG):
            return set()
        current_failures = set()
        in_run = False
        with open(PROGRESS_LOG) as f:
            for line in f:
                parts = [p.strip() for p in line.split('|', 2)]
                if len(parts) < 3:
                    continue
                _ts, status, rest = parts
                step = rest.split('(')[0].strip()  # strip any trailing note
                if status == 'START':
                    in_run = True
                    current_failures = set()
                elif in_run and status == 'FAILURE':
                    current_failures.add(step)
                elif in_run and status == 'SUCCESS':
                    current_failures.discard(step)
        return current_failures

    # -------------------------------------------------------------------------
    # Simulation command for testing order of everything:
    # def run_docker_cmd(cmd_arr, filename):
    #     # Simulate running the command by printing what would be run
    #     print(f'Running: {filename} with command {" ".join(cmd_arr)}')
    #     # Simulate execution time with a random delay
    #     time.sleep(2)
    #     print(f'Completed: {filename}')
    
    # --- Container cleanup on exit/interrupt ---
    _active_cid_files = []
    _cid_lock = threading.Lock()

    def _kill_all_containers():
        with _cid_lock:
            cids = list(_active_cid_files)
        for cid_file in cids:
            try:
                with open(cid_file) as f:
                    cid = f.read().strip()
                if cid:
                    _log(f"Cleaning up container {cid[:12]}...")
                    subprocess.run(['docker', 'kill', cid],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except Exception:
                pass

    def _signal_handler(_signum, _frame):
        _log("\nBuild interrupted — killing all running containers...")
        _kill_all_containers()
        sys.exit(1)

    atexit.register(_kill_all_containers)
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)

    def run_docker_cmd(cmd_arr, filename):
        '''
        Essentially a wrapper for 'docker run'. Also provides output.
        Tracks each container via --cidfile so they can be killed on exit/interrupt.
        With --continue: skips steps whose sentinel file already exists, and
        cleans up partial outputs after all retries are exhausted.
        '''
        # --- Resume: skip if already done ---
        if args.resume and _sentinel_done(filename):
            sentinel = _step_sentinel(filename)
            _log(f'[{filename}] skipping — already complete ({sentinel})')
            _write_progress(filename, 'SKIPPED', f'sentinel: {sentinel}')
            return

        retries = 3
        delays = [3 * 60, 10 * 60]  # 3 minutes, 10 minutes
        _log('running...' + filename)
        env = os.environ.copy()
        if 'SYNAPSE_AUTH_TOKEN' not in env.keys():
            _log('You need to set the SYNAPSE_AUTH_TOKEN to acess the MPNST, beatAML, bladder, pancreatic, liver, sarcoma, cnf datasets')
            docker_run = ['docker', 'run', '--rm', '-v', env['PWD']+'/local/:/tmp/', '--platform=linux/amd64']
        else:
            docker_run = ['docker', 'run', '--rm', '-v', env['PWD']+'/local/:/tmp/', '-e', 'SYNAPSE_AUTH_TOKEN='+env['SYNAPSE_AUTH_TOKEN'], '--platform=linux/amd64']

        attempt = 1
        while attempt <= retries:
            with tempfile.NamedTemporaryFile(suffix='.cid', delete=False) as _tf:
                cid_file = _tf.name
            os.remove(cid_file)  # Docker requires the file to not exist before writing
            cmd = docker_run + ['--cidfile', cid_file] + cmd_arr
            with _cid_lock:
                _active_cid_files.append(cid_file)
            try:
                _log(f"[{filename}] Attempt {attempt}/{retries}: {_redact(cmd)}")
                res = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
                if res.returncode == 0:
                    _log(f"[{filename}] succeeded on attempt {attempt}.")
                    _write_progress(filename, 'SUCCESS')
                    return
                else:
                    _log(f"[{filename}] failed (exit {res.returncode}).")
                    if attempt < retries:
                        delay = delays[attempt - 1]
                        _log(f"[{filename}] waiting {delay//60} minutes before retrying...")
                        time.sleep(delay)
            finally:
                with _cid_lock:
                    if cid_file in _active_cid_files:
                        _active_cid_files.remove(cid_file)
                try:
                    os.remove(cid_file)
                except OSError:
                    pass
            attempt += 1

        # All retries exhausted — log failure and remove partial outputs
        _write_progress(filename, 'FAILURE', f'all {retries} attempts failed')
        _cleanup_step(filename)
        raise RuntimeError(f"{filename} failed after {retries} attempts")

        
    
    def process_docker(datasets):
        '''
        Build specific docker images using docker-compose based on the dataset argument.
        All output and errors are logged at local/docker.log.
        
        Parameters:
        - datasets: list of datasets to process (e.g., ['broad_sanger', 'hcmi', 'mpnst'])
        '''
        compose_file = 'coderbuild/docker/docker-compose.yml'
        
        # Map datasets to corresponding Docker Containers
        dataset_map = {
            'broad_sanger': ['broad_sanger_exp', 'broad_sanger_omics'],
            'hcmi': ['hcmi'],
            'beataml': ['beataml'],
            'mpnst': ['mpnst'],
            'pancreatic': ['pancreatic'],
            'bladder': ['bladder'],
            'sarcoma': ['sarcoma'],
            'colorectal': ['colorectal'],
            'cptac': ['cptac'],
            'genes': ['genes'],
            'phosphosites': ['phosphosites'],
            'upload': ['upload'],
            'liver': ['liver'],
            'novartis': ['novartis'],
            'cnf': ['cnf']
        }
        
        # Collect container names to build based on the datasets provided. Always build genes and upload.
        datasets_to_build = ['genes', 'upload']
        if 'cnf' in datasets:
            datasets_to_build.append('phosphosites')
        for dataset in datasets:
            datasets_to_build.extend(dataset_map.get(dataset, []))
        
        # Build the docker-compose command, adding specific datasets
        compose_command = ['docker', 'compose', '-f', compose_file, 'build', '--parallel'] + datasets_to_build
        
        log_file_path = 'local/docker.log'
        env = os.environ.copy()
        
        _log(f"Docker-compose is building images for {', '.join(datasets_to_build)}. View output in {log_file_path}.")
        
        with open(log_file_path, 'w') as log_file:
            try:
                # Execute the docker-compose command
                res = subprocess.run(compose_command, env=env, stdout=log_file, stderr=log_file, text=True, check=True)
                log_file.write("Docker images built successfully.\n")
                _log(f"Docker images for {', '.join(datasets_to_build)} built successfully. Details logged in {log_file_path}.")
            except subprocess.CalledProcessError as e:
                log_file.write(f"Docker compose build failed with error: {e}\n")
                _log(f"Docker compose build failed. See {log_file_path} for details.")
                raise


    # target_class values in schema/expected_files.yaml, grouped by the build
    # step that produces them.
    _STEP_CLASSES = {
        'samples':     {'Sample'},
        'omics':       {'Transcriptomics', 'Proteomics', 'Mutations', 'Copy Number',
                        'Phosphoproteomics'},
        'drugs':       {'Drug', 'Drug Descriptor'},
        'experiments': {'Experiments', 'Combinations'},
    }

    def _verify_step_outputs(dataset, step):
        """Fail a step that exited 0 without writing the files it declares.

        build_all only checks a container's exit status, so a script that
        returns 0 after doing nothing counts as success. beataml omics did
        exactly that on 2026-09-13: it ran for 88 seconds (against 4m24s in a
        known-good build), printed "Starting Transcriptomics Data", never
        reached proteomics, exited 0, and wrote none of its three files. The
        build carried on through experiments, misc and three upload stages and
        only fell over ~11 hours later in validate.

        schema/expected_files.yaml already declares what each dataset must
        produce, so check it here and fail at the step that went wrong, while
        the log still shows why.
        """
        want = _STEP_CLASSES.get(step)
        if not want:
            return
        try:
            with open('schema/expected_files.yaml') as fh:
                declared = yaml.safe_load(fh)['datasets'].get(dataset, [])
        except Exception as exc:
            _log(f'[verify] could not read expected_files.yaml ({exc}); skipping output check.')
            return

        missing = []
        for entry in declared:
            if entry.get('target_class') not in want:
                continue
            base = os.path.basename(entry['file'])
            stem = base[:-3] if base.endswith('.gz') else base
            # the build writes either the plain or the gzipped form
            if not any(os.path.exists(os.path.join('local', c))
                       for c in (stem, stem + '.gz')):
                missing.append(stem)

        if missing:
            raise RuntimeError(
                f"{dataset} {step} reported success but produced none of: "
                f"{', '.join(missing)}. The step exited 0 without writing the "
                f"files schema/expected_files.yaml declares for it, so the build "
                f"is stopping here rather than failing later in validate.")

    def _await_all(futures, step_label):
        """Retrieve every future so no step failure can be silently dropped.

        This is load-bearing. Each process_* function used to await only the
        PREVIOUS future before submitting the next, so the LAST dataset's
        exception was never retrieved -- and an unretrieved exception in a
        ThreadPoolExecutor is discarded, not raised. In --high_mem mode nothing
        was awaited at all.

        On 2026-09-04 that let `hcmi omics` fail all three attempts while the
        build logged "All omics files completed" and carried on to publish and
        validate a release with hcmi transcriptomics missing entirely.

        Raises the first failure after collecting the rest, so one broken step
        cannot hide behind a later success.
        """
        errors = []
        for fut in futures:
            try:
                fut.result()
            except Exception as exc:
                errors.append(exc)
        if errors:
            if len(errors) > 1:
                _log(f'{step_label}: {len(errors)} steps failed; raising the first.')
            raise errors[0]

    def process_drugs(executor, datasets, prev_drugs=None):
        '''
        Build all drug files sequentially.

        prev_drugs: list of local paths (e.g. ['local_old/broad_sanger_drugs.tsv', ...])
            Files are copied into local/ and listed in local/prev_drug_files.txt.
            pubchem_retrieval.py auto-loads this file on import and uses it as a
            local cache: known drugs skip HTTP requests entirely, but ID assignment
            proceeds exactly as in a fresh run (dflist chain is unchanged).
        '''
        last_drug_future = None
        _futures = []
        dflist = []

        if prev_drugs:
            docker_paths = []
            for p in prev_drugs:
                base = os.path.basename(p)
                if base.endswith('.tsv.gz'):
                    stem = base[:-len('.tsv.gz')]
                    prev_name = f'{stem}_prev.tsv.gz'
                elif base.endswith('.tsv'):
                    stem = base[:-len('.tsv')]
                    prev_name = f'{stem}_prev.tsv'
                else:
                    stem = os.path.splitext(base)[0]
                    prev_name = f'{stem}_prev.tsv'
                shutil.copy(p, f'local/{prev_name}')
                docker_paths.append(f'/tmp/{prev_name}')
                _log(f'Copied prev drug cache: {p} → local/{prev_name}')
            with open('local/prev_drug_files.txt', 'w') as f:
                f.write('\n'.join(docker_paths) + '\n')
            _log(f'Wrote local/prev_drug_files.txt with {len(docker_paths)} cache file(s).')

        for da in datasets:
            if da not in ['cptac', 'hcmi']:
                di = 'broad_sanger_exp' if da == 'broad_sanger' else da
                if last_drug_future:
                    last_drug_future.result()  # Ensure the last drug process is completed before starting the next
                def _drugs_step(image=di, arg=','.join(dflist), ds=da):
                    run_docker_cmd([image, 'bash', 'build_drugs.sh', arg], f'{ds} drugs')
                    _verify_step_outputs(ds, 'drugs')
                last_drug_future = executor.submit(_drugs_step)
                _futures.append(last_drug_future)
                dflist.append(f'/tmp/{da}_drugs.tsv')
        _await_all(_futures, 'drugs')

    def process_samples(executor, datasets):
        '''
        Build all samples files sequentially
        '''
        last_sample_future = None
        _futures = []
        sf = ''
        for da in datasets:
            file_path = f'local/{da}_samples.csv'
            if os.path.exists(file_path):
                sf = f'/tmp/{da}_samples.csv'  # Set the most recent successfully processed dataset file
            else:
                break
        for da in datasets:
            di = 'broad_sanger_omics' if da == 'broad_sanger' else da
            if not os.path.exists(f'local/{da}_samples.csv'):
                if last_sample_future:
                    last_sample_future.result() 
                def _samples_step(image=di, arg=sf, ds=da):
                    run_docker_cmd([image, 'bash', 'build_samples.sh', arg], f'{ds} samples')
                    _verify_step_outputs(ds, 'samples')
                last_sample_future = executor.submit(_samples_step)
                _futures.append(last_sample_future)
                sf = f'/tmp/{da}_samples.csv'
        _await_all(_futures, 'samples')

    def process_phosphosites(executor):
        '''
        Build the phosphosites reference file if it does not exist.
        Caller must ensure genes.csv exists first. Returns a Future (or None).
        '''
        if not os.path.exists('local/phosphosites.csv'):
            return executor.submit(
                run_docker_cmd,
                ['phosphosites', 'bash', 'build_phosphosites.sh', '/tmp/genes.csv'],
                'phosphosites file',
            )
        return None

    def process_omics(executor, datasets, high_mem):
        '''
        Build all omics files concurrently
        '''
        last_omics_future = None
        _futures = []
        for da in datasets:
            di = 'broad_sanger_omics' if da == 'broad_sanger' else da
            omics_cmd = [di, 'bash', 'build_omics.sh', '/tmp/genes.csv', f'/tmp/{da}_samples.csv']
            if da == 'cnf':
                omics_cmd.append('/tmp/phosphosites.csv')
            def _omics_step(cmd=omics_cmd, ds=da):
                run_docker_cmd(cmd, f'{ds} omics')
                _verify_step_outputs(ds, 'omics')
            if high_mem:
                _futures.append(executor.submit(_omics_step))
            else:
                if last_omics_future:
                    last_omics_future.result()
                last_omics_future = executor.submit(_omics_step)
                _futures.append(last_omics_future)
        _await_all(_futures, 'omics')

    def process_experiments(executor, datasets, high_mem):
        '''
        Build all experiments files concurrently
        '''
        last_experiments_future = None
        _futures = []
        for da in datasets:
            if da not in ['cptac', 'hcmi']:
                di = 'broad_sanger_exp' if da == 'broad_sanger' else da
                if not os.path.exists(f'local/{da}_experiments.tsv'):
                    #Run all at once
                    if high_mem:
                        _futures.append(executor.submit(run_docker_cmd, [di, 'bash', 'build_exp.sh', f'/tmp/{da}_samples.csv', f'/tmp/{da}_drugs.tsv'], f'{da} experiments'))
                    #Run one at a time
                    else:
                        if last_experiments_future:
                            last_experiments_future.result() 
                        last_experiments_future = executor.submit(run_docker_cmd, [di, 'bash', 'build_exp.sh', f'/tmp/{da}_samples.csv', f'/tmp/{da}_drugs.tsv'], f'{da} experiments')
                        _futures.append(last_experiments_future)
        _await_all(_futures, 'experiments')

    def process_misc(executor, datasets, high_mem):
        '''
        Run all misc scripts concurrently or one at a time.
        '''
        last_misc_future = None
        _futures = []
        #Currently this only applies to broad_sanger. Add others here if they need a final step.
        if "broad_sanger" in datasets:
            datasets = ["broad_sanger"]
        else:
            return
        for da in datasets:
            #Running the build_misc.sh in broad_sanger_omics
            di = 'broad_sanger_omics' if da == 'broad_sanger' else da
            #Run all at once:
            if high_mem:
                _futures.append(executor.submit(run_docker_cmd, [di, 'bash', 'build_misc.sh'], f'{da} misc'))
            #Run one at a time.
            else:
                if last_misc_future:
                    last_misc_future.result() 
                last_misc_future = executor.submit(run_docker_cmd,  [di, 'bash', 'build_misc.sh'], f'{da} misc')
                _futures.append(last_misc_future)
        _await_all(_futures, 'misc')

    def process_genes(executor):
        if not os.path.exists('local/genes.csv'):
            return executor.submit(run_docker_cmd,['genes','bash','build_genes.sh'],'gene file')
        return None
        
        
    def run_docker_upload_cmd(cmd_arr, all_files_dir, name, version):
        '''
        Wrapper for 'docker run'. This one is focused on uploads.
        '''
        env = os.environ.copy()
        docker_run = ['docker', 'run', '--rm', '-v', f"{env['PWD']}/local/{all_files_dir}:/tmp", '-e', f"VERSION={version}"]

        # Add Appropriate Environment Variables
        if name == "validate":
            docker_run.extend(['upload'])
        if 'FIGSHARE_TOKEN' in env and name == 'Figshare':
            docker_run.extend(['-e', f"FIGSHARE_TOKEN={env['FIGSHARE_TOKEN']}", 'upload'])
        if name in ["Map_Drugs", "Map_Samples", "Align_Drug_Descriptors"]:
            docker_run.extend(['upload'])
        if 'GITHUB_TOKEN' in env and name == "GitHub":
            docker_run.extend(['-e', f"GITHUB_TOKEN={env['GITHUB_TOKEN']}", 'upload'])

        # Full command to run including version update
        docker_run.extend(cmd_arr)
        _log('Executing:', _redact(docker_run))
        # res = subprocess.run(docker_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res = subprocess.run(docker_run, stdout=sys.stdout, stderr=sys.stderr)
        if res.returncode != 0:
            # Do NOT log res.stderr here. stderr is passed through to sys.stderr
            # above, so subprocess does not capture it and res.stderr is always
            # None -- which is exactly what this line used to print, hiding the
            # real cause of a `validate failed` behind the word "None".
            # The container's own output is already interleaved above.
            _log(f'{name} failed: exited with code {res.returncode}.'
                 f' See the container output immediately above for the cause.')
            exit(f'{name} failed')
        else:
            _log(f'{name} successful')
            

    def decompress_file(file_path):
        """Decompress a gzip file and delete the original compressed file."""
        with gzip.open(file_path, 'rb') as f_in:
            decompressed_file_path = file_path[:-3]  # Remove '.gz' from the filename
            with open(decompressed_file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_path)

    def compress_file(file_path):
        """Compress a file using gzip and delete the original uncompressed file."""
        compressed_file_path = file_path + '.gz'
        with open(file_path, 'rb') as f_in:
            with gzip.open(compressed_file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_path)
        
    def get_latest_commit_hash(owner, repo, branch='main'):
        """
        Returns the SHA of the latest commit on the specified branch.
        """
        url = f"https://api.github.com/repos/{owner}/{repo}/commits/{branch}"
        response = requests.get(url)
        response.raise_for_status()
        
        # The commit data is in JSON format; the 'sha' field is the full commit hash.
        commit_data = response.json()
        return commit_data['sha']
            
    ######
    ### Pre-Build Environment Token Check
    #####

    figshare_token = os.getenv('FIGSHARE_TOKEN')
    synapse_auth_token = os.getenv('SYNAPSE_AUTH_TOKEN')
    github_token = os.getenv('GITHUB_TOKEN')


    # Error handling for required tokens
    if args.figshare and not figshare_token:
        raise ValueError("FIGSHARE_TOKEN environment variable is not set.")
    # NOTE: broad_sanger now needs a Synapse token too -- its DepMap inputs are
    # read from Synapse (syn75028495) because the DepMap portal added a
    # Cloudflare challenge. See the DepMap gate below.
    if any(dataset in args.datasets for dataset in ['beataml', 'mpnst', 'bladder', 'pancreatic','sarcoma','liver','novartis','colorectal','cnf','broad_sanger']) and not synapse_auth_token:
        if args.docker or args.samples or args.omics or args.drugs or args.exp or args.all: # Token only required if building data, not upload or validate.
            raise ValueError("SYNAPSE_AUTH_TOKEN is required for accessing MPNST, beatAML, bladder, pancreatic, liver, novartis, colorectal, sarcoma, cnf, broad_sanger (DepMap files) datasets.")

    # -------------------------------------------------------------------------
    # DepMap freshness gate.
    #
    # DepMap put a Cloudflare Turnstile ("verify you are a person") challenge in
    # front of https://depmap.org/portal/api/download/files, so the broad_sanger
    # samples and omics steps can no longer fetch DepMap data from the portal.
    # coderdata instead reads a static copy from the "DepMap Raw" Synapse folder
    # (syn75028495); coderbuild/utils/fetch_depmap.py downloads it inside the
    # container at the start of build_samples.sh and build_omics.sh.
    #
    # That copy must be refreshed BY HAND for each new DepMap release. Because a
    # stale copy fails silently -- the build succeeds and quietly publishes the
    # previous release's cell lines -- we require an explicit confirmation flag
    # rather than trusting that someone remembered.
    #
    # Only the steps that actually read the DepMap files are gated (samples and
    # omics); --docker on its own still works for image debugging.
    # -------------------------------------------------------------------------
    DEPMAP_SYNAPSE_FOLDER = 'syn75028495'
    DEPMAP_REQUIRED_FILES = [
        'Model.csv',
        'OmicsSomaticMutations.csv',
        'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv',
        'PortalOmicsCNGeneLog2.csv',
    ]
    if 'broad_sanger' in args.datasets and (args.samples or args.omics or args.all):
        if not args.depmap_ready:
            raise SystemExit(
                "\n"
                "=====================================================================\n"
                " STOP: DepMap files must be confirmed up to date before building\n"
                "=====================================================================\n"
                " This build includes broad_sanger, which requires DepMap data.\n"
                "\n"
                " DepMap files can NO LONGER be downloaded from the DepMap portal --\n"
                " it now serves a Cloudflare 'verify you are a person' challenge.\n"
                " coderdata reads a STATIC COPY from Synapse instead, which must be\n"
                " updated MANUALLY for each new DepMap release.\n"
                "\n"
                f"   Synapse folder: https://www.synapse.org/Synapse:{DEPMAP_SYNAPSE_FOLDER}\n"
                "\n"
                " To refresh the files for a new DepMap release:\n"
                "   1. Download them from the DepMap download page:\n"
                "        https://depmap.org/portal/data_page/?tab=allData\n"
                "      (accept the terms, then use the download button on each file)\n"
                f"   2. Upload them to {DEPMAP_SYNAPSE_FOLDER} as NEW VERSIONS of the\n"
                "      existing entities (this keeps the syn IDs stable)\n"
                "   3. Bump DEPMAP_RELEASE in coderbuild/utils/fetch_depmap.py\n"
                "\n"
                " Files tracked (4 total, ~1.29 GB):\n"
                + ''.join(f"   - {f}\n" for f in DEPMAP_REQUIRED_FILES) +
                "\n"
                " Repository: https://github.com/PNNL-CompBio/coderdata\n"
                "\n"
                " Once you have confirmed the DepMap files on Synapse are current,\n"
                " re-run this command with:\n"
                "\n"
                "     --depmap-ready\n"
                "\n"
                "=====================================================================\n"
            )
        _log(f"DepMap files confirmed current on Synapse ({DEPMAP_SYNAPSE_FOLDER}) "
             f"via --depmap-ready; fetch_depmap.py will download them in-container.")


    ######
    ### Begin Pipeline
    #####
    
    # Make a 'local' directory for output
    if not os.path.exists('local'):
        os.mkdir('local')

    # Record this run in the progress log
    _write_progress('BUILD START', 'START', f'args={sys.argv[1:]}')

    # Get dataset names - default is all.
    datasets = args.datasets.split(',')

    ### Build Docker Images. These are all built in Parallel. Nothing else can run until these are built.
    # Ouput is logged at local/docker.log
    if args.docker or args.all:
        process_docker(datasets)
        _log("Docker image generation completed")

    # --continue: rebuild Docker images only for datasets that failed last run
    if args.resume:
        failed_steps = _last_run_failures()
        if failed_steps:
            failed_datasets = {
                _step_to_docker_dataset(s) for s in failed_steps
                if _step_to_docker_dataset(s)
            }
            # Only rebuild images for datasets we're actually processing
            failed_datasets &= set(datasets) | {'genes', 'phosphosites'}
            if failed_datasets:
                _log(f'--continue: rebuilding Docker images for failed datasets: {sorted(failed_datasets)}')
                process_docker(list(failed_datasets))
            else:
                _log('--continue: no failed datasets from last run match current dataset list — skipping Docker rebuild')
        else:
            _log('--continue: no failures found in build_progress.log — skipping Docker rebuild')
        

    ### Build Drugs files, Samples files, and Genes file. These two steps are run in Parallel.
    ### Within each step, sequential running is required.
    with ThreadPoolExecutor() as executor:
        if args.samples or args.all:
            sample_thread = executor.submit(process_samples,executor, datasets)
        if args.drugs or args.all:
            prev_drugs = [p.strip() for p in args.prev_drugs.split(',') if p.strip()]
            drug_thread = executor.submit(process_drugs, executor, datasets, prev_drugs or None)

        # Genes must finish before phosphosites can start (phosphosites reads genes.csv).
        # Run genes now and wait; then submit phosphosites (which can overlap with samples/drugs).
        if args.samples or args.omics or args.exp or args.all:
            genes_future = process_genes(executor)
            if genes_future is not None:
                genes_future.result()

        phosphosite_future = None
        if (args.omics or args.all) and 'cnf' in datasets:
            phosphosite_future = process_phosphosites(executor)

        # Wait for all remaining tasks to complete before proceeding to omics and experiments
        if args.drugs or args.all:
            drug_thread.result()
        if args.samples or args.all:
            sample_thread.result()
        if phosphosite_future is not None:
            phosphosite_future.result()

    _log("All samples, drugs files, genes, and phosphosites files completed or skipped")


    ### At this point in the pipeline, all samples and drugs files have been created. There are no blockers to proceed.
    ### Build Omics files and Experiments files. These two steps are run in Parallel. 
    ### Within each step, all datasets are run in Parallel.
    
    with ThreadPoolExecutor() as executor:
        if args.omics or args.all:
            omics_thread = executor.submit(process_omics, executor, datasets, args.high_mem)
        if args.exp or args.all:
            exp_thread = executor.submit(process_experiments, executor, datasets, args.high_mem)
            
        if args.omics or args.all:
            omics_thread.result()
            _log("All omics files completed")
        if args.exp or args.all:
            exp_thread.result()
            _log("All experiments files completed")

    # Omics + experiments are done. Drop their large raw intermediates now
    # (e.g. rnaseq_all_data_*.csv ~5GB, cell.xml, proteomics_all_*.csv), but keep
    # broad*/sanger* — build_misc.sh still needs them to split broad_sanger.
    _clean_local_intermediates(reason='after omics/experiments',
                               protect_prefixes=('broad', 'sanger'))


    ### Final Step, some datasets may need an additional post build step. Add this here
    # Currently only the cell line datasets need this. This seperates broad_sanger into all of its component datasets.
    
    with ThreadPoolExecutor() as executor:
        if args.misc or args.all:
            misc_thread = executor.submit(process_misc, executor, datasets, args.high_mem)
        if args.misc or args.all:
            misc_thread.result()
            _log("Final build step complete.")

    # misc has split broad_sanger into its per-dataset outputs, so the broad*/
    # sanger* intermediates are no longer needed. Clear all remaining
    # intermediates from local/ before the disk-heavy map/align/upload steps
    # (this is where the previous build hit "No space left on device").
    _clean_local_intermediates(reason='after misc (pre-upload)')


    ######
    ### Begin Upload and/or validation
    #####
    if args.figshare or args.validate or github_token:
    # if args.figshare or args.validate:
        # FigShare File Prefixes:
        
        broad_sanger_datasets = ["ccle","ctrpv2","fimm","gdscv1","gdscv2","gcsi","prism","nci60"]
        if "broad_sanger" in datasets:
            datasets.extend(broad_sanger_datasets)
            datasets.remove("broad_sanger")

        figshare_token = os.getenv('FIGSHARE_TOKEN')

        all_files_dir = 'local/all_files_dir'
        if not os.path.exists(all_files_dir):
            os.makedirs(all_files_dir)
        
        # Ensure figshare tokens are available
        if  args.figshare and not figshare_token:
            raise ValueError("Required tokens (FIGSHARE) are not set in environment variables.")
        
        # Ensure version is specified
        if args.figshare and not args.version:
            raise ValueError("Version must be specified when pushing to figshare")

        # Build exact allowlist from schema/expected_files.yaml so only
        # schema-declared output files (not _prev caches, docker.log, etc.)
        # end up in all_files_dir.
        with open('schema/expected_files.yaml') as _f:
            _schema = yaml.safe_load(_f)
        expected_basenames = set()
        for _entries in _schema.get('datasets', {}).values():
            for _entry in _entries:
                _bn = os.path.basename(_entry['file'])
                expected_basenames.add(_bn)
                expected_basenames.add(_bn + '.gz')

        # Move only expected files to a designated directory
        for file in glob(os.path.join("local", '*.*')):
            if os.path.basename(file) in expected_basenames:
                shutil.move(file, os.path.join(all_files_dir, os.path.basename(file)))

        # Decompress all compressed files in the directory for schema checking
        for file in glob(os.path.join(all_files_dir, '*.gz')):
            decompress_file(file)

        ### These should be done before schema checking.
        version_args = ['--version', args.version] if args.version is not None else []
        sample_mapping_command = ['python3', 'scripts/map_improve_sample_ids.py', '--local_dir', "/tmp"] + version_args
        run_docker_upload_cmd(sample_mapping_command, 'all_files_dir', 'Map_Samples', args.version)

        drug_mapping_command = ['python3', 'scripts/map_improve_drug_ids.py', '--local_dir', "/tmp"] + version_args
        run_docker_upload_cmd(drug_mapping_command, 'all_files_dir', 'Map_Drugs', args.version)

        drug_mapping_command_2 = ['python3', 'scripts/align_drug_descriptors.py', '--local_dir', "/tmp"] + version_args
        run_docker_upload_cmd(drug_mapping_command_2, 'all_files_dir', 'Align_Drug_Descriptors', args.version)

        # Run schema checker - This will always run if uploading data.
        schema_check_command = ['python3', 'scripts/check_schema.py', '--datasets'] + datasets
        run_docker_upload_cmd(schema_check_command, 'all_files_dir', 'validate', args.version)
        
        _log("Validation complete. Proceeding with file compression/decompression adjustments")
        
        # Compress or decompress files based on specific conditions after checking
        for file in glob(os.path.join(all_files_dir, '*')):
            is_compressed = file.endswith('.gz')
            if ('samples' in file or 'figshare' in file) and is_compressed:
                decompress_file(file)
            elif not ('samples' in file or 'figshare' in file) and not is_compressed:
                compress_file(file)

        _log("File compression and decompression adjustments are complete.")
    
        ### Upload to Figshare using Docker
        if args.figshare and args.version and figshare_token:
            figshare_command = ['python3', 'scripts/push_to_figshare.py', '--directory', "/tmp", '--title', f"CODERData{args.version}", '--token', os.getenv('FIGSHARE_TOKEN'), '--project_id', '189342', '--version', args.version, '--publish']
            run_docker_upload_cmd(figshare_command, 'all_files_dir', 'Figshare', args.version)

            ### Push changes to GitHub using Docker
            # if args.version and args.figshare and figshare_token and github_token and args.github_username and args.github_email:
            
            # You can only upload to Github after Figshare upload is completed - otherwise figshare_latest.yml and dataset.yml won't be available.
            if args.version and github_token and args.github_username and args.github_email:

                git_command = [
                    'bash', '-c', (
                        f'git config --global user.name "{args.github_username}" '
                        f'&& git config --global user.email "{args.github_email}" '
                        
                        # Checkout a new branch
                        f'&& git checkout -b testing-auto-build-pr-{args.version} '
                        
                        # Copy and add the necessary files
                        f'&& cp /tmp/improve_sample_mapping.json.gz /usr/src/app/coderdata/coderbuild/improve_sample_mapping.json.gz '
                        f'&& cp /tmp/improve_drug_mapping.json.gz /usr/src/app/coderdata/coderbuild/improve_drug_mapping.json.gz '
                        f'&& gunzip /usr/src/app/coderdata/coderbuild/*.gz '
                        f'&& git add -f coderbuild/improve_sample_mapping.json coderbuild/improve_drug_mapping.json '
                        f'&& cp /tmp/figshare_latest.yml /usr/src/app/coderdata/docs/_data/figshare_latest.yml '
                        f'&& cp /tmp/dataset.yml /usr/src/app/coderdata/coderdata/dataset.yml '
                        f'&& git add -f docs/_data/figshare_latest.yml coderdata/dataset.yml'
                        
                        # Tag and push
                        f'&& git commit -m "Data Built and Uploaded. New Tag: {args.version}" '
                        f'&& git tag {args.version} '
                        f'&& git push https://{args.github_username}:{github_token}@github.com/PNNL-CompBio/coderdata.git testing-auto-build-pr-{args.version} '
                        
                        # Create a PR using GitHub CLI
                        f'&& gh pr create --title "Testing Auto PR instead of auto Merge {args.version}" '
                        f'--body "This PR was automatically generated by the build process." '
                        f'--base main --head testing-auto-build-pr-{args.version}'
                    )
                ]
            
                run_docker_upload_cmd(git_command, 'all_files_dir', 'GitHub', args.version)
            
if __name__ == '__main__':
    main()
