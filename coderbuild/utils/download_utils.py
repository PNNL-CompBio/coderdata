"""Resumable, verified HTTP downloads for the build scripts.

The NCI wiki serves DOSERESP.zip at ~346MB. A plain urlopen().read() of that
either succeeds or raises

    http.client.IncompleteRead: IncompleteRead(32601897 bytes read,
                                               313702904 more expected)

which is exactly how build v32 lost the whole nci60 dataset: the download
truncated at 32MB, 04b-nci60-updated.py exited 1, the caller only warned, and
nci60DoseResponse was never written -- so curve fitting saw 7 studies instead
of 8 and nci60 shipped as eight empty files that still validated clean.

Retrying from zero is not enough for a file this size on a link that drops
part-way; it tends to fail again at a similar point. Each attempt therefore
RESUMES with a Range request. Verified that wiki.nci.nih.gov sends
accept-ranges: bytes and answers a Range request with 206.
"""

import os
import time
from urllib import request, error

# 1, 3, 10 and 15 minutes, matching coderbuild/utils/retry_utils.R.
RETRY_SLEEPS = (60, 180, 600, 900)

# The NCI wiki blocks Python's default User-Agent.
BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)


def _remote_size(url, headers):
    try:
        req = request.Request(url, headers=headers, method="HEAD")
        with request.urlopen(req, timeout=60) as resp:
            n = resp.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:
        return None


def retrieve_url(url, dest, sleeps=RETRY_SLEEPS, headers=None):
    """Download `url` to `dest`, resuming an interrupted transfer.

    Bytes accumulate in "<dest>.part", which is renamed onto `dest` only once
    the finished size matches the server's Content-Length, so a truncated file
    never appears at the destination. The ".part" is kept between attempts --
    deleting it would mean a flaky link could never accumulate a large file.
    """
    headers = dict(headers or {"User-Agent": BROWSER_UA})
    part = dest + ".part"
    expected = _remote_size(url, headers)

    last_error = None
    for attempt in range(len(sleeps) + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        try:
            hdrs = dict(headers)
            if have:
                hdrs["Range"] = f"bytes={have}-"
            req = request.Request(url, headers=hdrs)
            with request.urlopen(req, timeout=120) as resp:
                resumed = resp.status == 206
                if have and not resumed:
                    # Server ignored the Range; start over rather than append
                    # onto bytes it is about to send again.
                    have = 0
                mode = "ab" if (have and resumed) else "wb"
                with open(part, mode) as f:
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)

            got = os.path.getsize(part)
            if expected is not None and got != expected:
                raise IOError(f"incomplete: {got} of {expected} bytes")
            os.replace(part, dest)
            return dest

        except (error.URLError, error.HTTPError, IOError, OSError) as e:
            last_error = e
            got = os.path.getsize(part) if os.path.exists(part) else 0
            if attempt < len(sleeps):
                wait = sleeps[attempt]
                print(f"  download of {url} failed at {got:,} bytes ({e}); "
                      f"attempt {attempt + 1}/{len(sleeps) + 1}, resuming in "
                      f"{wait // 60} min", flush=True)
                time.sleep(wait)

    raise RuntimeError(
        f"Failed to download {url} after {len(sleeps) + 1} attempts: {last_error}")
