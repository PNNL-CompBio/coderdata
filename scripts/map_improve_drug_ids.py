#!/usr/bin/env python3
import json
import os
import csv
import argparse
from datetime import datetime
import gzip
import shutil
from collections import defaultdict


##### Load, make structure, save functions for improve_drug_mapping.json

def load_mapping(mapping_file='improve_drug_mapping.json'):
    """
    Loads an existing improve_drug_mapping.json if available.
    Otherwise returns an empty base structure with metadata and drugs.
    """
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            return json.load(f), True
    else:
        return {
            "metadata": {
                "builds": []
            },
            "drugs": []
        }, False

def save_mapping(mapping_data):
    """Saves mapping data to disk as JSON."""
    mapping_file='/tmp/improve_drug_mapping.json'
    with open(mapping_file, 'w') as f:
        json.dump(mapping_data, f, indent=2)

def get_current_build_metadata(build_date, version):
    """Returns dict describing current build metadata."""
    return {
        "build_date": build_date,
        # "git_commit": git_commit,
        # "description": description,
        "version": version
    }

#### (De)Compression Functions

def decompress_gz_if_needed(file_path):
    """
    If file_path ends with .gz, decompress into a temp file and return that path,
    plus a bool indicating it was gz compressed. Otherwise return original path, False.
    """
    if file_path.endswith('.gz'):
        decompressed_path = file_path[:-3]
        try:
            with gzip.open(file_path,'rb') as f_in, open(decompressed_path,'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
            return decompressed_path, True
        except Exception as e:
            print(f"Error decompressing {file_path}: {e}")
            return file_path, False
    return file_path, False

def recompress_if_needed(original_path, decompressed_path, was_gz):
    """
    If was_gz is True, re-compress decompressed_path back into original_path
    and remove decompressed_path. Otherwise do nothing.
    """
    if was_gz:
        try:
            with open(decompressed_path,'rb') as f_in, gzip.open(original_path,'wb',compresslevel=5) as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(decompressed_path)
        except Exception as e:
            print(f"Error recompressing {decompressed_path} to {original_path}: {e}")


##### Read {dataset}_drugs.tsv files

def read_drugs_tsv(file_path):
    """
    Reads a {dataset}_drugs.tsv or .tsv.gz and returns a list of row dicts.
    """
    if not os.path.exists(file_path):
        gz_path = file_path + '.gz'
        if os.path.exists(gz_path):
            file_path = gz_path
        else:
            print(f"File not found: {file_path} or {gz_path}")
            return []

    file_path, was_gz = decompress_gz_if_needed(file_path)
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        recompress_if_needed(file_path, file_path, was_gz)
        print(f"Empty or missing file after decompression: {file_path}")
        return []

    rows = []
    print(f"Reading drugs file: {file_path}")
    with open(file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)

    recompress_if_needed(file_path, file_path, was_gz)
    return rows


##### Create or pick stable IDs for drugs

def generate_new_stable_id(drugs_list):
    """
    Generate a new stable ID in the form SMI_{n} by examining existing stable_ids that 
    either match 'SMI_{#}' or are purely numeric, then picking the next integer.
    """
    existing_nums = []
    for drug in drugs_list:
        sid = drug["stable_id"]
        # If stable_id is numeric or in the form SMI_{number}, parse out the number portion
        if sid.isdigit():
            existing_nums.append(int(sid))
        elif sid.upper().startswith("SMI_"):
            suffix = sid[4:]
            if suffix.isdigit():
                existing_nums.append(int(suffix))
    if existing_nums:
        return "SMI_" + str(max(existing_nums) + 1)
    else:
        return "SMI_1"


def unify_drugs(mapping_data, all_drugs_rows):
    """
    Unify each canSMILES to a stable_id, dropping build_statuses.
    This logic is analogous to unify_samples, but simpler:
      - We keep track of canSMILES -> stable_id
      - If new canSMILES is found, try using the FIRST improve_drug_id as stable_id (if free & not numeric)
        else use SMI_{#}
      - Return drug_id_mapping (old improve_drug_id -> stable_id)
    """

    drugs_list = mapping_data["drugs"]

    # Build canSMILES -> stable_id from existing data
    canSMILES_to_stable_id = {}
    for drug in drugs_list:
        canSMILES_to_stable_id[drug["canSMILES"]] = drug["stable_id"]

    # Track used stable_ids so we don't accidentally duplicate
    used_stable_ids = set(d["stable_id"] for d in drugs_list)

    # Group new data by canSMILES
    new_drugs_grouped = defaultdict(list)  # canSMILES -> list of improve_drug_ids
    for row in all_drugs_rows:
        imp_id = row.get("improve_drug_id", "").strip()
        csmiles = row.get("canSMILES", "").strip()
        if imp_id and csmiles:
            new_drugs_grouped[csmiles].append(imp_id)

    drug_id_mapping = {}
    next_candidate = generate_new_stable_id(drugs_list)

    for csmiles, drug_ids in new_drugs_grouped.items():
        if csmiles in canSMILES_to_stable_id:
            # Already known canSMILES => use existing stable_id
            stable_id = canSMILES_to_stable_id[csmiles]
            for did in drug_ids:
                drug_id_mapping[did] = stable_id
        else:
            # New canSMILES => new stable_id
            # Try using the FIRST improve_drug_id if it's free and not numeric
            first_id = drug_ids[0]
            if (first_id not in used_stable_ids) and (not first_id.isdigit()):
                stable_id = first_id
            else:
                stable_id = next_candidate
                # Generate a new candidate for future
                tmp_drug_list = drugs_list + [{"stable_id": stable_id}]
                next_candidate = generate_new_stable_id(tmp_drug_list)

            used_stable_ids.add(stable_id)
            canSMILES_to_stable_id[csmiles] = stable_id

            # Map all the improve_drug_ids to that stable_id
            for did in drug_ids:
                drug_id_mapping[did] = stable_id

            # Add to the JSON, minus any build_statuses
            new_drug = {
                "stable_id": stable_id,
                "canSMILES": csmiles
            }
            drugs_list.append(new_drug)
            print(f"Created new stable_id '{stable_id}' for canSMILES '{csmiles}'.")

    return drug_id_mapping


def rewrite_drugs_file(file_path, drug_id_mapping, datasets=None):
    """
    Rewrites {dataset}_drugs.tsv or {dataset}_drug_descriptors.tsv:
      - Replace improve_drug_id with stable_id
      - Only rewrite if file is in the specified datasets (if provided)
    """
    filename = os.path.basename(file_path)
    if datasets and not any(ds in filename for ds in datasets):
        return

    print(f"Rewriting drug file: {file_path}")
    file_path, was_gz = decompress_gz_if_needed(file_path)
    if not os.path.exists(file_path):
        recompress_if_needed(file_path, file_path, was_gz)
        print(f"File not found or empty after decompression: {file_path}")
        return

    with open(file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = list(reader)
        if not rows:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"Empty file: {file_path}")
            return

        header = rows[0]
        if "improve_drug_id" not in header:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"'improve_drug_id' column not found in {file_path}")
            return
        try:
            idx_id = header.index("improve_drug_id")
        except ValueError:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"'improve_drug_id' column index error in {file_path}")
            return

    tmp_path = file_path + ".tmp"
    with open(file_path, 'r', newline='', encoding='utf-8') as fin, \
         open(tmp_path, 'w', newline='', encoding='utf-8') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        hdr = next(reader)
        writer.writerow(hdr)

        for row in reader:
            if len(row) <= idx_id:
                writer.writerow(row)
                continue
            original_id = row[idx_id].strip()
            if original_id in drug_id_mapping:
                new_id = drug_id_mapping[original_id]
                print(f"Replacing improve_drug_id '{original_id}' with stable_id '{new_id}' in {file_path}")
                row[idx_id] = new_id
            writer.writerow(row)

    os.replace(tmp_path, file_path)
    recompress_if_needed(file_path, file_path, was_gz)


def main():
    parser = argparse.ArgumentParser(description="""
Use canSMILES overlaps to assign stable drug IDs (without build_statuses).
On the first build, create improve_drug_mapping.json but do not rewrite.
Subsequent builds unify via canSMILES, update stable IDs, and rewrite files.
""")
    parser.add_argument('--build_date', default=None,
                        help='Build date in YYYY-MM-DD. Default=now.')
    # parser.add_argument('--git_commit', default='abcdef12345',
    #                     help='Git commit hash.')
    # parser.add_argument('--description', default='New build with updated drugs.',
    #                     help='Build description.')
    parser.add_argument('--version', required=True,
                        help='Build version. Must be unique per build.')
    parser.add_argument('--datasets', default='broad_sanger,ccle,ctrpv2,fimm,gcsi,gdscv1,gdscv2,nci60,prism,beataml,mpnst,mpnstpdx',
                        help='Comma-separated list of datasets.')
    parser.add_argument('--local_dir', default='data',
                        help='Directory containing TSV files.')
    parser.add_argument('--other_files', default='drug_descriptors,experiments',
                        help='Comma-separated list of other file types to rewrite.')
    args = parser.parse_args()

    # Set build_date
    build_date = args.build_date or datetime.utcnow().strftime("%Y-%m-%d")

    # Load or initialize improve_drug_mapping.json
    mapping_file = "improve_drug_mapping.json"
    mapping_data, had_prior = load_mapping(mapping_file)

    # Insert current build metadata
    current_build_metadata = get_current_build_metadata(build_date, args.version)

    # Ensure build uniqueness
    if any(
        b["build_date"] == build_date and b["version"] == args.version
        for b in mapping_data["metadata"]["builds"]
    ):
        print(f"Build {build_date}, version {args.version} already exists in {mapping_file}. Exiting.")
        return

    mapping_data["metadata"]["builds"].append(current_build_metadata)

    ds_list = [d.strip() for d in args.datasets.split(',')]
    other_file_types = [f.strip() for f in args.other_files.split(',')]
    local_dir = args.local_dir

    # Gather all drugs from each dataset
    all_drugs_rows = []
    for ds in ds_list:
        path_tsv = os.path.join(local_dir, f"{ds}_drugs.tsv")
        rows = read_drugs_tsv(path_tsv)
        all_drugs_rows.extend(rows)

    # First build check
    is_first_build = (not had_prior) and (not mapping_data["drugs"])
    if is_first_build:
        print("First build detected. Generating improve_drug_mapping.json without rewriting files.")

        # We'll unify once, but skip rewriting
        drug_id_mapping = unify_drugs(mapping_data, all_drugs_rows)
        save_mapping(mapping_data)
        print(f"Created {mapping_file} with {len(mapping_data['drugs'])} drugs. No rewrite performed.")
        return

    # Otherwise, subsequent build
    print("Subsequent build detected. Unifying new drugs and rewriting files.")
    drug_id_mapping = unify_drugs(mapping_data, all_drugs_rows)

    # Sort final "drugs" by numeric portion if SMI_#, or lexically
    def parse_smi_num(sid):
        # E.g., stable_id = "SMI_42" => 42
        if sid.upper().startswith("SMI_"):
            suffix = sid[4:]
            if suffix.isdigit():
                return int(suffix)
        # If purely numeric, handle that
        if sid.isdigit():
            return int(sid)
        return 999999999

    mapping_data["drugs"].sort(key=lambda d: parse_smi_num(d["stable_id"]))
    save_mapping(mapping_data, mapping_file)
    print(f"Updated {mapping_file} with {len(mapping_data['drugs'])} drugs.")

    # Rewrite the main drug TSV
    for ds in ds_list:
        path_tsv = os.path.join(local_dir, f"{ds}_drugs.tsv")
        if os.path.exists(path_tsv) or os.path.exists(path_tsv + '.gz'):
            rewrite_drugs_file(path_tsv, drug_id_mapping, datasets=ds_list)

    # Rewrite other file types (e.g., drug_descriptors)
    for ds in ds_list:
        for ftype in other_file_types:
            other_path = os.path.join(local_dir, f"{ds}_{ftype}.tsv")
            if os.path.exists(other_path) or os.path.exists(other_path + '.gz'):
                rewrite_drugs_file(other_path, drug_id_mapping, datasets=ds_list)

    print("All files have been rewritten with stable IDs (no build_statuses).")
    print(f"Stable IDs have been updated in {mapping_file}.")


if __name__ == "__main__":
    main()
