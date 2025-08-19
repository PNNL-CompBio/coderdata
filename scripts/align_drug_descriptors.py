#!/usr/bin/env python3
import os
import gzip
import shutil
import csv
import argparse

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
def build_reference_map(files):
    """Build map from (drug_id, descriptor) with first found value."""
    ref = {}
    for fp in files:
        path, gz = decompress_gz_if_needed(fp)
        with open(path, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                key = (row['improve_drug_id'], row['structural_descriptor'])
                val = row['descriptor_value']
                if key not in ref:
                    ref[key] = val
        recompress_if_needed(path, gz, fp)
    return ref

def rewrite_files(files, ref):
    """Go back and rewrite any mismatches in-place."""
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
                key = (row['improve_drug_id'], row['structural_descriptor'])
                correct = ref.get(key)
                if correct is not None and row['descriptor_value'] != correct:
                    # print(f"Fixing {key} in {os.path.basename(fp)}: "
                    #       f"{row['descriptor_value']} to {correct}")
                    row['descriptor_value'] = correct
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
    parser.add_argument('--version',   help=argparse.SUPPRESS)  # ignore the version input
    args = parser.parse_args()
    files = find_descriptor_files(args.local_dir)
    if not files:
        print("No drug_descriptor files found in", args.local_dir)
        return

    ref = build_reference_map(files)
    rewrite_files(files, ref)
    print("Done.")

if __name__ == '__main__':
    main()
