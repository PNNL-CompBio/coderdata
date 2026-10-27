#!/usr/bin/env python3
"""
Fetch the DepMap files required by the broad_sanger build from Synapse.

WHY THIS EXISTS
---------------
DepMap put a Cloudflare Turnstile ("verify you are a person") challenge in front
of https://depmap.org/portal/api/download/files. That endpoint now returns an
HTML challenge page with HTTP 200, so the previous programmatic download silently
produced garbage: read_csv() parsed the HTML, and the build failed downstream with
a confusing "object 'filename' not found". DepMap also explicitly asks that the
portal not be scraped, so this is not something to work around.

coderdata therefore keeps a STATIC COPY of the required DepMap files on Synapse.
This script downloads them into a container-local scratch dir (see
DEPMAP_WORKDIR) before the broad_sanger samples and omics R scripts run. Those
R scripts read them from there, delete each one once consumed, and never touch
the DepMap portal.

UPDATING FOR A NEW DEPMAP RELEASE (must be done by hand)
--------------------------------------------------------
1. Download the files from the DepMap download page:
       https://depmap.org/portal/data_page/?tab=allData
   Accept the terms, then use the download button on each file below.
2. Upload them into the "DepMap Raw" Synapse folder
       https://www.synapse.org/Synapse:syn75028495
   as NEW VERSIONS of the existing entities listed in DEPMAP_SYN_IDS.
   Uploading as a new version keeps the syn ID stable, so nothing in this file
   needs to change. (If you create brand-new Synapse entities instead, update
   the IDs in DEPMAP_SYN_IDS below.)
3. Bump DEPMAP_RELEASE below to the new release name.
4. Run the build with --depmap-ready to confirm the Synapse copy is current.

Repository: https://github.com/PNNL-CompBio/coderdata

STANDALONE USE
--------------
Run this directly to verify Synapse access and the configured IDs before
launching a long build:

    SYNAPSE_AUTH_TOKEN=... python3 fetch_depmap.py --dest /tmp/depmap_check
"""

import argparse
import os
import sys

# ---------------------------------------------------------------------------
# Configuration -- edit DEPMAP_RELEASE and DEPMAP_SYN_IDS when DepMap ships a
# new release (see "UPDATING FOR A NEW DEPMAP RELEASE" above).
#
# Every value can also be overridden by environment variable, which is how
# other coderdata builders (cnf, liver, bladder, ...) parameterise Synapse IDs.
# ---------------------------------------------------------------------------
DEPMAP_RELEASE = os.environ.get('DEPMAP_RELEASE', '26Q1')

# Synapse entity IDs for the four files the broad_sanger build consumes.
# These live in the "DepMap Raw" folder: https://www.synapse.org/Synapse:syn75028495
#
# Upload each new DepMap release as a NEW VERSION of these same entities and
# the IDs below never need to change -- syn.get() always returns the latest
# version.
DEPMAP_SYNAPSE_FOLDER = 'syn75028495'  # "DepMap Raw" -- for reference/uploads

DEPMAP_SYN_IDS = {
    'Model.csv':
        os.environ.get('DEPMAP_MODEL_SYN_ID', 'syn75028594'),
    'OmicsSomaticMutations.csv':
        os.environ.get('DEPMAP_MUTATIONS_SYN_ID', 'syn77171195'),
    'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv':
        os.environ.get('DEPMAP_TRANSCRIPTOMICS_SYN_ID', 'syn77171175'),
    'PortalOmicsCNGeneLog2.csv':
        os.environ.get('DEPMAP_COPY_NUMBER_SYN_ID', 'syn77171184'),
}

DOWNLOAD_PAGE = 'https://depmap.org/portal/data_page/?tab=allData'
REPO_URL = 'https://github.com/PNNL-CompBio/coderdata'

# Container-local scratch directory for the downloaded DepMap files.
#
# Deliberately NOT under /tmp: build_all.py bind-mounts the host's local/
# directory at /tmp inside every container, so anything written there lands on
# the host (and, here, inside a OneDrive-synced folder) and survives the run.
# These files are transient build inputs -- they are consumed by the R scripts
# and deleted -- so they belong in the container's own filesystem, which is
# discarded when the `docker run --rm` exits.
DEPMAP_WORKDIR = os.environ.get('DEPMAP_DIR', '/opt/depmap_data')


def default_dest():
    """Container-local directory the R scripts read the DepMap files from."""
    return DEPMAP_WORKDIR


def _unconfigured():
    """Names of files whose Synapse ID has not been set yet."""
    return [name for name, sid in DEPMAP_SYN_IDS.items() if not sid.strip()]


def main():
    parser = argparse.ArgumentParser(
        description='Download the static DepMap files for broad_sanger from Synapse.'
    )
    parser.add_argument('--dest', default=None,
                        help='Directory to download into '
                             '(default: $DEPMAP_DIR or /opt/depmap_data, container-local)')
    parser.add_argument('--print-dest', action='store_true',
                        help='Print the download directory and exit. Lets the '
                             'build shell scripts derive DEPMAP_DIR without '
                             'duplicating the release name.')
    args = parser.parse_args()

    dest = args.dest or default_dest()

    if args.print_dest:
        print(dest)
        return

    missing_ids = _unconfigured()
    if missing_ids:
        sys.exit(
            '\n'
            '=====================================================================\n'
            ' DepMap Synapse IDs are not configured\n'
            '=====================================================================\n'
            f' No Synapse entity ID is set for:\n'
            + ''.join(f'   - {n}\n' for n in missing_ids) +
            '\n'
            ' DepMap files can no longer be downloaded programmatically from the\n'
            ' DepMap portal, so coderdata reads a static copy from Synapse.\n'
            '\n'
            ' To fix this:\n'
            f'   1. Download the files from {DOWNLOAD_PAGE}\n'
            '   2. Upload them to Synapse\n'
            '   3. Put their syn IDs in DEPMAP_SYN_IDS in\n'
            '      coderbuild/utils/fetch_depmap.py (or set the matching\n'
            '      DEPMAP_*_SYN_ID environment variables)\n'
            '\n'
            f' Repository: {REPO_URL}\n'
            '=====================================================================\n'
        )

    if not os.environ.get('SYNAPSE_AUTH_TOKEN'):
        sys.exit('SYNAPSE_AUTH_TOKEN is not set; cannot download DepMap files '
                 'from Synapse.')

    try:
        import synapseclient
    except ImportError:
        sys.exit("synapseclient is not installed in this image; it is required "
                 "to fetch the DepMap files from Synapse. Add 'synapseclient' "
                 "to coderbuild/broad_sanger/requirements.txt and rebuild.")

    os.makedirs(dest, exist_ok=True)

    syn = synapseclient.Synapse()
    syn.login()

    print(f'Fetching DepMap {DEPMAP_RELEASE} files from Synapse into {dest}')
    for fname, syn_id in DEPMAP_SYN_IDS.items():
        print(f'  {fname} <- {syn_id}')
        # Always re-download so the copy used by the build provably matches
        # Synapse. `dest` is container-local, so this costs no host disk and is
        # discarded when the container exits.
        entity = syn.get(
            syn_id,
            downloadLocation=dest,
            ifcollision='overwrite.local',
        )

        # Synapse names the downloaded file after the entity, which may differ
        # from the name the R scripts expect. Normalise it.
        got = getattr(entity, 'path', None)
        if not got or not os.path.exists(got):
            sys.exit(f'Synapse returned no downloaded file for {fname} ({syn_id}). '
                     f'Check that the entity exists and that this account has '
                     f'download permission on {DEPMAP_SYNAPSE_FOLDER}.')
        want = os.path.join(dest, fname)
        if os.path.abspath(got) != os.path.abspath(want):
            os.replace(got, want)

        size = os.path.getsize(want)
        if size == 0:
            sys.exit(f'Downloaded {fname} from {syn_id} but it is empty.')
        print(f'    -> {want} ({size / 1e6:.1f} MB)')

    print(f'All {len(DEPMAP_SYN_IDS)} DepMap {DEPMAP_RELEASE} files ready in {dest}')


if __name__ == '__main__':
    main()
