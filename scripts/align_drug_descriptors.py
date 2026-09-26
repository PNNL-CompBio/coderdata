#!/usr/bin/env python3
"""Harmonize drug_descriptor values across datasets.

For each (improve_drug_id, structural_descriptor) seen across the per-dataset
`*_drug_descriptors.tsv[.gz]` files, the FIRST value encountered (files processed
in sorted order) is treated as canonical, and any later file whose value differs
is rewritten to match.

Memory-bounded implementation
-----------------------------
The reference map is kept in an on-disk SQLite database, not a Python dict, so the
step scales to very large descriptor files (e.g. nci60 ~89M rows / 2.7 GB). The
previous in-memory dict held every (drug_id, descriptor) -> value pair at once and
was OOM-killed by the container on large builds.

Fingerprint descriptors are skipped from harmonization: a fingerprint is a
deterministic function of chemical structure, so the same drug always yields the
same value -- there is nothing to harmonize, and these rows dominate both memory
and file size. They are passed through unchanged.
"""
import os
import sys
import gzip
import shutil
import csv
import argparse
import sqlite3
import tempfile

# Descriptor values (e.g. Morgan fingerprints) can be long strings.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

# Descriptor names (case-insensitive substring match) excluded from harmonization.
SKIP_DESCRIPTOR_SUBSTRINGS = ("fingerprint",)

_INSERT_BATCH = 100_000


def _skip_descriptor(name):
    n = name.lower()
    return any(s in n for s in SKIP_DESCRIPTOR_SUBSTRINGS)


# Helper scripts
def decompress_gz_if_needed(path):
    """If path ends with .gz, decompress to a temp file and return its name plus True."""
    if path.endswith('.gz'):
        out = path[:-3]
        with gzip.open(path, 'rb') as f_in, open(out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return out, True
    return path, False


def recompress_if_needed(decompressed, was_gz, original):
    """If was_gz, recompress decompressed back to original and remove decompressed."""
    if was_gz:
        with open(decompressed, 'rb') as f_in, gzip.open(original, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(decompressed)


def find_descriptor_files(directory):
    """Find all *_drug_descriptors.tsv and .tsv.gz files in directory."""
    files = []
    for fn in os.listdir(directory):
        if fn.endswith('_drug_descriptors.tsv') or fn.endswith('_drug_descriptors.tsv.gz'):
            files.append(os.path.join(directory, fn))
    return sorted(files)


# Actual work
def build_reference_db(files, conn):
    """Populate SQLite `ref` with the first-seen value per (drug_id, descriptor).

    INSERT OR IGNORE keeps the first value encountered, matching the original
    'first found wins' semantics. Fingerprint descriptors are not stored.
    """
    cur = conn.cursor()
    cur.execute(
        "CREATE TABLE ref (drug_id TEXT, descriptor TEXT, value TEXT, "
        "PRIMARY KEY (drug_id, descriptor)) WITHOUT ROWID"
    )
    for fp in files:
        path, gz = decompress_gz_if_needed(fp)
        with open(path, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            batch = []
            for row in reader:
                desc = row['structural_descriptor']
                if _skip_descriptor(desc):
                    continue
                batch.append((row['improve_drug_id'], desc, row['descriptor_value']))
                if len(batch) >= _INSERT_BATCH:
                    cur.executemany("INSERT OR IGNORE INTO ref VALUES (?,?,?)", batch)
                    batch.clear()
            if batch:
                cur.executemany("INSERT OR IGNORE INTO ref VALUES (?,?,?)", batch)
        conn.commit()
        recompress_if_needed(path, gz, fp)


def rewrite_files(files, conn):
    """Rewrite any mismatched descriptor_value in-place to the canonical value."""
    cur = conn.cursor()
    for fp in files:
        path, gz = decompress_gz_if_needed(fp)
        tmp = path + '.tmp'
        changed = False

        with open(path, newline='', encoding='utf-8') as fin, \
             open(tmp, 'w', newline='', encoding='utf-8') as fout:

            reader = csv.DictReader(fin, delimiter='\t')
            writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter='\t')
            writer.writeheader()

            for row in reader:
                desc = row['structural_descriptor']
                if not _skip_descriptor(desc):
                    cur.execute(
                        "SELECT value FROM ref WHERE drug_id=? AND descriptor=?",
                        (row['improve_drug_id'], desc),
                    )
                    res = cur.fetchone()
                    if res is not None and row['descriptor_value'] != res[0]:
                        row['descriptor_value'] = res[0]
                        changed = True
                writer.writerow(row)

        if changed:
            os.replace(tmp, path)
        else:
            os.remove(tmp)

        recompress_if_needed(path, gz, fp)


def main():
    parser = argparse.ArgumentParser(
        description="Harmonize drug_descriptor values across multiple files."
    )
    parser.add_argument('--local_dir', default='.', help='Folder containing *_drug_descriptors.tsv[.gz]')
    parser.add_argument('--version', help=argparse.SUPPRESS)  # ignore the version input
    args = parser.parse_args()
    files = find_descriptor_files(args.local_dir)
    if not files:
        print("No drug_descriptor files found in", args.local_dir)
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        conn = sqlite3.connect(os.path.join(tmpdir, "descriptors.sqlite"))
        # Temp DB: durability is irrelevant, so favour speed.
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        try:
            build_reference_db(files, conn)
            rewrite_files(files, conn)
        finally:
            conn.close()
    print("Done.")


if __name__ == '__main__':
    main()
