import os
import time
import pandas as pd
import argparse
from zipfile import ZipFile
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

# Waits between whole-download attempts: 1, 3, 10 and 15 minutes, matching
# coderbuild/utils/retry_utils.R. Long enough to ride out a real incident
# rather than only a momentary blip; if all attempts fail the build fails.
RETRY_SLEEPS = (60, 180, 600, 900)


def robust_download(url, dest_path, sleeps=RETRY_SLEEPS):
    """Download `url` to `dest_path`, resuming an interrupted transfer.

    The previous implementation mounted a urllib3 Retry adapter and assumed that
    made it robust. It did not: Retry only covers connection setup and HTTP
    status codes. Once the response body is streaming, a truncated body raises
    ChunkedEncodingError/ProtocolError out of iter_content(), which went straight
    to the except clause and became a fatal RuntimeError -- no retry, no resume.

    That is exactly how build v24 died, 12MB into a 112MB file:

      Failed to download .../Proteomics_20221214.zip:
      ('Connection broken: IncompleteRead(12353536 bytes read,
        100318054 more expected)')

    So: retry the whole transfer, and RESUME it with a Range request rather than
    starting over, since restarting a large file on a flaky link tends to fail
    again at a similar point. Verified that this host sends accept-ranges: bytes
    and answers a Range request with 206.

    Bytes accumulate in "<dest>.part", which is renamed onto dest_path only once
    the finished size matches the server's Content-Length. An unverified or
    partial file therefore never appears at the destination.
    """
    part = dest_path + ".part"

    retry_strategy = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET", "HEAD"},
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    expected = None
    try:
        head = session.head(url, timeout=(5, 30), allow_redirects=True)
        if head.ok and head.headers.get("Content-Length"):
            expected = int(head.headers["Content-Length"])
    except requests.exceptions.RequestException:
        pass  # size unknown; completeness is then judged by a clean stream end

    last_error = None
    for attempt in range(len(sleeps) + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        try:
            headers = {"Range": f"bytes={have}-"} if have else {}
            with session.get(url, stream=True, timeout=(5, 60), headers=headers) as r:
                if have and r.status_code == 206:
                    mode = "ab"                     # server honoured the resume
                else:
                    if have:
                        print(f"  server ignored Range for {url}; restarting", flush=True)
                    r.raise_for_status()
                    mode, have = "wb", 0
                with open(part, mode) as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:                   # filter out keep-alive chunks
                            f.write(chunk)

            got = os.path.getsize(part)
            if expected is not None and got != expected:
                raise IOError(f"incomplete: {got} of {expected} bytes")

            os.replace(part, dest_path)
            return

        except (requests.exceptions.RequestException, IOError) as e:
            last_error = e
            got = os.path.getsize(part) if os.path.exists(part) else 0
            # Keep the .part file: the next attempt resumes from here. Deleting
            # it would mean a flaky link could never accumulate a large file.
            if attempt < len(sleeps):
                wait = sleeps[attempt]
                print(f"  download of {url} failed at {got} bytes "
                      f"({e}); attempt {attempt + 1}/{len(sleeps) + 1}, "
                      f"resuming in {wait // 60} min", flush=True)
                time.sleep(wait)

    raise RuntimeError(
        f"Failed to download {url} after {len(sleeps) + 1} attempts: {last_error}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', dest='samplefile', default=None, help='DepMap sample file')
    parser.add_argument('--gene', dest='genefile', default=None, help='Gene file')

    opts = parser.parse_args()

    samplefile = opts.samplefile
    gfile = opts.genefile

    samps = pd.read_csv(samplefile)
    print(samps)
    genes = pd.read_csv(gfile)[['gene_symbol','entrez_id']]
    genes = genes.drop_duplicates()

    print(genes)
    protfile='https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx'


    
    prots = pd.read_excel(protfile,'Normalized Protein Expression')
    print(prots)


    vvars = [col for col in prots.columns if '_Ten' in col]

    prot2 = prots[['Gene_Symbol']+vvars]

    ##we can just do the metlt here
    plong = pd.melt(prot2,id_vars='Gene_Symbol',value_vars=vvars,var_name='cellline',value_name='proteomics')

    ##rename gene symbol column
    plong = plong.rename({'Gene_Symbol':'gene_symbol'},axis=1)
    print(plong)
    
    ##split cell lin
    plong['other_id'] = [a.split('_Ten')[0] for a in plong.cellline]

    full = plong.merge(genes,on='gene_symbol')
    full = full.merge(samps,on='other_id')

    full = full.loc[:,['entrez_id','proteomics','improve_sample_id']].drop_duplicates().dropna()

    full[['study']] = 'DepMap'
    full[['source']] = 'Broad'
    ##now save to separate files
    full.dropna(axis=0)
    full.to_csv('/tmp/broad_proteomics.csv.gz', index=False, compression='gzip')


    #old download, (was failing too much)
    # sanger_protfile='https://cog.sanger.ac.uk/cmp/download/Proteomics_20221214.zip'
    # r = requests.get(sanger_protfile)
    # sanger_loc ='/tmp/sp.zip'
    # open(sanger_loc , 'wb').write(r.content)
    # zf = ZipFile(sanger_loc,'r')
    # zf.extractall(path='/tmp/')
    # pdat = pd.read_csv('/tmp/Protein_matrix_averaged_zscore_20221214.tsv',sep='\t',skiprows=[0])
    
    ##now get sanger
    sanger_protfile = "https://cog.sanger.ac.uk/cmp/download/Proteomics_20221214.zip"
    sanger_loc = "/tmp/sp.zip"
    robust_download(sanger_protfile, sanger_loc)
    with ZipFile(sanger_loc, "r") as zf:
        zf.extractall(path="/tmp/")
    pdat = pd.read_csv(
        "/tmp/Protein_matrix_averaged_zscore_20221214.tsv",
        sep="\t",
        skiprows=[0],
    )
    
    vv=pdat.columns[2:]
    plong = pd.melt(pdat,id_vars='symbol',value_vars=vv)
    pres = plong.rename({'symbol':'other_names','variable':'gene_symbol','value':'proteomics'},axis=1)
    pres = pres.merge(genes,on='gene_symbol')
    pres = pres.merge(samps,on='other_names')

    full2 = pres.loc[:,['entrez_id','improve_sample_id','proteomics']].drop_duplicates().dropna()
    full2.loc[:,['study']] = 'Sanger'
    full2.loc[:,['source']] = 'Sanger'
    
    #full3 = pd.concat([full,full2])
    #print(full3)
    full2.dropna(axis=0)
    full2.to_csv('/tmp/sanger_proteomics.csv.gz',index=False, compression='gzip')
    
main()
