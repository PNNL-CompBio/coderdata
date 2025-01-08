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

def save_mapping(mapping_data, mapping_file):
    """Saves mapping data to disk as JSON."""
    with open(mapping_file, 'w') as f:
        json.dump(mapping_data, f, indent=2)

def get_current_build_metadata(build_date, version):
    """Returns dict describing current build metadata."""
    return {
        "build_date": build_date,
        "version": version
    }

#### (De)Compression Functions
def decompress_gz_if_needed(file_path):
    """
    If file_path ends with .gz, decompress into a temp file and return that path,
    plus a bool indicating it was gz compressed.
    """
    if file_path.endswith('.gz'):
        decompressed_path = file_path[:-3]
        try:
            with gzip.open(file_path, 'rb') as f_in, open(decompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            return decompressed_path, True  # Return decompressed path, was_gz=True
        except Exception as e:
            print(f"Error decompressing {file_path}: {e}")
            return None, False
    return file_path, False

def recompress_if_needed(decompressed_path, was_gz, original_path=None):
    """
    If was_gz is True, re-compress decompressed_path back into original_path
    and remove decompressed_path. Otherwise do nothing.
    """
    if was_gz and original_path:
        try:
            with open(decompressed_path, 'rb') as f_in, gzip.open(original_path, 'wb', compresslevel=5) as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(decompressed_path)
        except Exception as e:
            print(f"Error recompressing {decompressed_path} to {original_path}: {e}")

##### Read {dataset}_drugs.tsv files
def read_drugs_tsv(file_path, dataset):
    """
    Reads a {dataset}_drugs.tsv or .tsv.gz and returns a list of row dicts.
    Adds the dataset name to each row.
    """
    if not os.path.exists(file_path):
        gz_path = file_path + '.gz'
        if os.path.exists(gz_path):
            file_path = gz_path
        else:
            print(f"File not found: {file_path} or {gz_path}")
            return []
    file_path, was_gz = decompress_gz_if_needed(file_path)
    if file_path is None or not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        recompress_if_needed(file_path, was_gz)
        print(f"Empty or missing file after decompression: {file_path}")
        return []
    rows = []
    print(f"Reading drugs file: {file_path}")
    with open(file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            row['dataset'] = dataset  # Add dataset info
            rows.append(row)
    recompress_if_needed(file_path, was_gz, original_path=file_path + '.gz' if was_gz else None)
    return rows

##### Create or pick stable IDs for drugs
def generate_new_stable_id(assigned_stable_ids):
    """
    Generate a new stable ID in the form SMI_{n} by examining existing stable_ids that 
    match 'SMI_{#}', then picking the next integer.
    """
    existing_nums = []
    for sid in assigned_stable_ids:
        if sid.upper().startswith("SMI_"):
            suffix = sid[4:]
            if suffix.isdigit():
                existing_nums.append(int(suffix))
    if existing_nums:
        return "SMI_" + str(max(existing_nums) + 1)
    else:
        return "SMI_1"

def unify_drugs(mapping_data, all_drugs_rows, dataset_priority):
    """
    Unify each canSMILES to a stable_id, prioritizing datasets.
    """
    drugs_list = mapping_data["drugs"]
    # Build canSMILES -> stable_id from existing data
    canSMILES_to_stable_id = {drug["canSMILES"]: drug["stable_id"] for drug in drugs_list}
    # Build stable_id -> canSMILES mapping from existing data
    stable_id_to_canSMILES = {drug["stable_id"]: drug["canSMILES"] for drug in drugs_list}
    # Track used stable_ids to avoid duplication
    assigned_stable_ids = set(drug["stable_id"] for drug in drugs_list)
    # Group new data by canSMILES
    new_drugs_grouped = defaultdict(set)  # canSMILES -> set of (improve_drug_id, dataset)
    imp_id_to_datasets = defaultdict(set)  # improve_drug_id -> set of datasets
    for row in all_drugs_rows:
        imp_id = row.get("improve_drug_id", "").strip()
        csmiles = row.get("canSMILES", "").strip()
        dataset = row.get("dataset", "").strip()
        if imp_id and csmiles and dataset:
            new_drugs_grouped[csmiles].add((imp_id, dataset))
            imp_id_to_datasets[imp_id].add(dataset)
    # Map improve_drug_id to stable_id for replacements
    drug_id_mapping = {}
    # Process each canSMILES
    for csmiles, imp_id_dataset_pairs in new_drugs_grouped.items():
        if csmiles in canSMILES_to_stable_id:
            # Already known canSMILES
            stable_id = canSMILES_to_stable_id[csmiles]
        else:
            # New canSMILES
            # Collect (priority, improve_drug_id) pairs
            imp_id_info = []
            for imp_id, dataset in imp_id_dataset_pairs:
                priority = dataset_priority.get(dataset, float('inf'))
                imp_id_info.append((priority, imp_id, dataset))
            # Sort by priority (lower number is higher priority)
            imp_id_info.sort()
            # Try to assign an improve_drug_id as stable_id based on priority
            stable_id = None
            for _, imp_id, _ in imp_id_info:
                if imp_id not in assigned_stable_ids:
                    stable_id = imp_id
                    break
            if not stable_id:
                # Assign a new unique stable_id
                stable_id = generate_new_stable_id(assigned_stable_ids)
            # Update mappings
            canSMILES_to_stable_id[csmiles] = stable_id
            stable_id_to_canSMILES[stable_id] = csmiles
            assigned_stable_ids.add(stable_id)
            # Add to mapping_data
            drugs_list.append({
                "stable_id": stable_id,
                "canSMILES": csmiles
            })
            print(f"Assigned stable_id '{stable_id}' to new canSMILES '{csmiles}'.")
        # Map improve_drug_ids to stable_id
        for imp_id, dataset in imp_id_dataset_pairs:
            if imp_id != stable_id:
                drug_id_mapping[(dataset, imp_id)] = stable_id
            else:
                # Even if imp_id == stable_id, add to mapping
                drug_id_mapping[(dataset, imp_id)] = stable_id
    return drug_id_mapping

def rewrite_drugs_file(file_path, drug_id_mapping, datasets=None):
    """
    Rewrites drug files by replacing improve_drug_id with stable_id based on drug_id_mapping.
    """
    filename = os.path.basename(file_path)
    dataset = filename.split('_')[0]
    if datasets and dataset not in datasets:
        return
    print(f"Processing drug file: {file_path}")
    file_path, was_gz = decompress_gz_if_needed(file_path)
    if file_path is None or not os.path.exists(file_path):
        recompress_if_needed(file_path, was_gz, original_path=file_path + '.gz' if was_gz else None)
        print(f"File not found or empty after decompression: {file_path}")
        return
    changes_made = False
    tmp_path = file_path + ".tmp"
    with open(file_path, 'r', newline='', encoding='utf-8') as fin, \
         open(tmp_path, 'w', newline='', encoding='utf-8') as fout:
        reader = csv.DictReader(fin, delimiter='\t')
        fieldnames = reader.fieldnames
        if "improve_drug_id" not in fieldnames:
            recompress_if_needed(file_path, was_gz, original_path=file_path + '.gz' if was_gz else None)
            print(f"'improve_drug_id' column not found in {file_path}")
            return
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            original_id = row.get("improve_drug_id", "").strip()
            key = (dataset, original_id)
            if key in drug_id_mapping:
                new_id = drug_id_mapping[key]
                if original_id != new_id:
                    print(f"Replacing improve_drug_id '{original_id}' with stable_id '{new_id}' in {file_path}")
                    row["improve_drug_id"] = new_id
                    changes_made = True
            writer.writerow(row)
    if changes_made:
        os.replace(tmp_path, file_path)
        print(f"Changes were made to {file_path}. File has been updated.")
    else:
        os.remove(tmp_path)
        print(f"No changes made to {file_path}. File remains unchanged.")
    recompress_if_needed(file_path, was_gz, original_path=file_path + '.gz' if was_gz else None)

def parse_smi_num(sid):
    """Parse numeric part of stable_id for sorting."""
    # E.g., stable_id = "SMI_42" => 42
    if sid.upper().startswith("SMI_"):
        suffix = sid[4:]
        if suffix.isdigit():
            return int(suffix)
    # If purely numeric, handle that
    if sid.isdigit():
        return int(sid)
    return None

def main():
    parser = argparse.ArgumentParser(description="""
Use canSMILES overlaps to assign stable drug IDs (without build_statuses).
Prioritize datasets to retain improve_drug_id values from higher-priority datasets.
Unify via canSMILES, update stable IDs, and rewrite files.
""")
    parser.add_argument('--build_date', default=None,
                        help='Build date in YYYY-MM-DD. Default=now.')
    parser.add_argument('--version', required=True,
                        help='Build version. Must be unique per build.')
    parser.add_argument('--datasets', default='gdscv1,broad_sanger,ccle,ctrpv2,fimm,gcsi,gdscv2,nci60,prism,beataml,mpnst,mpnstpdx',
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
    # Prepare dataset and file type lists
    ds_list = [d.strip() for d in args.datasets.split(',')]
    other_file_types = [x.strip() for x in args.other_files.split(',')]
    local_dir = args.local_dir
    # Define dataset priorities based on the order they appear in ds_list
    # Lower index means higher priority (starting from 1)
    dataset_priority = {ds: idx + 1 for idx, ds in enumerate(ds_list)}
    print("Dataset priorities:", dataset_priority)
    # Gather all drugs from each dataset
    all_drugs_rows = []
    for ds in ds_list:
        path_tsv = os.path.join(local_dir, f"{ds}_drugs.tsv")
        rows = read_drugs_tsv(path_tsv, ds)
        all_drugs_rows.extend(rows)
    # Unify drugs and get mappings
    print("Unifying drugs and updating mappings.")
    drug_id_mapping = unify_drugs(mapping_data, all_drugs_rows, dataset_priority)
    # Sort final "drugs" by numeric portion if SMI_#, or lexically
    mapping_data["drugs"].sort(key=lambda d: parse_smi_num(d["stable_id"]))
    save_mapping(mapping_data, mapping_file)
    print(f"Updated {mapping_file} with {len(mapping_data['drugs'])} drugs.")
    # Proceed to rewrite files
    print("Rewriting files with updated stable IDs.")
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