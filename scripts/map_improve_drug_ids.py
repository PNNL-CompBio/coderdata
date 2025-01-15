#!/usr/bin/env python3
import json
import os
import csv
import argparse
from datetime import datetime
import gzip
import shutil

##### Load, make structure, save functions for improve_drug_mapping.json
def load_mapping(mapping_file='improve_drug_mapping.json'):
    """
    Loads an existing improve_drug_mapping.json if available.
    Otherwise, returns an empty base structure with metadata and drugs.
    """
    if os.path.exists(mapping_file) and os.path.getsize(mapping_file) > 0:
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
            print(f"Loaded mapping data from {mapping_file}.")
            return mapping_data, True
    else:
        print(f"No existing mapping file found or file is empty: {mapping_file}. Initializing new mapping data.")
        return {
            "metadata": {
                "builds": []
            },
            "drugs": []
        }, False

def save_mapping(mapping_data, mapping_file='improve_drug_mapping.json'):
    """Saves mapping data to disk as JSON."""
    with open(mapping_file, 'w') as f:
        json.dump(mapping_data, f, indent=2)
    print(f"Saved mapping data to {mapping_file}.")

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

##### Helper to parse numeric part from SMI IDs
def parse_smi_num(sid):
    """Parse numeric part of stable_id for sorting (e.g. SMI_123 -> 123)."""
    if sid and sid.upper().startswith("SMI_"):
        suffix = sid[4:]
        if suffix.isdigit():
            return int(suffix)
    elif sid and sid.isdigit():
        return int(sid)
    return None


##### Resolve conflicts so only one canSMILES gets an improve_drug_id
def resolve_id_conflicts(all_drugs_rows, dataset_priority):
    """
    Pre-processing step:
      - If an imp_id (SMI_xxx) is used by more than 1 distinct canSMILES,
        keep it ONLY for the highest-priority usage.
      - All other usages become "losing IDs": __loser_id__ is set to the old ID,
        and removed from improve_drug_id so unify_drugs will assign them a new ID.
    Returns: a modified list of rows (copies) so we don't mutate the original in place.
    """
    rows_copy = []
    for row in all_drugs_rows:
        rows_copy.append(dict(row))  # shallow copy

    # imp_id_map: dict of imp_id -> set of (canSMILES, dataset)
    imp_id_map = {}
    for r in rows_copy:
        row_imp_id = r.get('improve_drug_id')
        imp_id_str = row_imp_id.strip() if isinstance(row_imp_id, str) else ''
        csm = r.get('canSMILES', '').strip()
        ds = r.get('dataset', '').strip()

        if imp_id_str and csm and ds:
            if imp_id_str not in imp_id_map:
                imp_id_map[imp_id_str] = set()
            imp_id_map[imp_id_str].add((csm, ds))

    # For each imp_id with more than 1 distinct canSMILES, pick a single winner
    for imp_id, pairs in imp_id_map.items():
        # Only if it starts with "SMI_", so we handle stable IDs
        if not imp_id.upper().startswith("SMI_"):
            continue

        # Distinct canSMILES used by this ID
        canSmi_set = {p[0] for p in pairs}
        if len(canSmi_set) > 1:
            # Conflict: multiple canSMILES share this SMI_xxx
            # Sort them by dataset priority => find the "winner"
            tmp = []
            for (csm, ds) in pairs:
                tmp.append((dataset_priority.get(ds, float('inf')), csm, ds))
            tmp.sort(key=lambda x: x[0])  # by priority ascending
            # The winner is the first entry (lowest priority number)
            _, winner_csm, winner_ds = tmp[0]

            # All other duplicates are forced to get a new ID
            # but we store the old ID in __loser_id__ so we can rewrite later
            for r in rows_copy:
                row_imp_id = r.get('improve_drug_id')
                row_imp_id_str = row_imp_id.strip() if isinstance(row_imp_id, str) else ''
                r_csm = r.get('canSMILES', '').strip()
                r_ds = r.get('dataset', '').strip()
                if row_imp_id_str == imp_id and r_csm != winner_csm:
                    # Mark the old ID in a new field so unify_drugs can produce a mapping
                    r['__loser_id__'] = imp_id
                    # Remove from improve_drug_id so unify_drugs forces a new stable ID
                    r['improve_drug_id'] = None

    return rows_copy

# Main unification function
def unify_drugs(mapping_data, all_drugs_rows, dataset_priority):
    """
    Unify each canSMILES to a stable_id, prioritizing datasets.
    Only reassign an ID if there's a genuine conflict or if it's impossible
    to reuse the original ID.

    Now, if a row has __loser_id__ (the old ID from a conflict),
    we build a mapping from (dataset, old_id) -> new_stable_id
    so rewriting can happen in the original file.
    """
    drugs_list = mapping_data["drugs"]

    # Build canSMILES -> stable_id from existing data
    canSMILES_to_stable_id = {drug["canSMILES"]: drug["stable_id"] for drug in drugs_list}
    
    # Track assigned stable_ids (from existing data)
    assigned_stable_ids = set(drug["stable_id"] for drug in drugs_list)

    # Determine the highest numeric SMI_xxx from previous file
    current_max_stable_id = 0
    for sid in assigned_stable_ids:
        if sid and sid.upper().startswith("SMI_"):
            suffix = sid[4:]
            if suffix.isdigit():
                num = int(suffix)
                if num > current_max_stable_id:
                    current_max_stable_id = num

    # Determine the highest numericSMI_xxx from new data
    max_file_id = 0
    smi_rows = [
        row for row in all_drugs_rows
        if row.get("improve_drug_id") and row["improve_drug_id"].upper().startswith("SMI_")
    ]
    for item in smi_rows:
        suffix = item["improve_drug_id"][4:]
        if suffix.isdigit():
            num = int(suffix)
            if num > max_file_id:
                max_file_id = num

    next_available_id = max(current_max_stable_id, max_file_id)

    # Build imp_id_to_canSMILES to check uniqueness
    imp_id_to_canSMILES = {}
    for row in all_drugs_rows:
        imp_id = row.get("improve_drug_id", "")
        csmiles = row.get("canSMILES", "")
        if imp_id and csmiles:
            if imp_id not in imp_id_to_canSMILES:
                imp_id_to_canSMILES[imp_id] = set()
            imp_id_to_canSMILES[imp_id].add(csmiles)

    # Group new data by canSMILES => set of (imp_id, dataset, loser_id)
    new_drugs_grouped = {}
    for row in all_drugs_rows:
        imp_id = row.get("improve_drug_id", "")
        csmiles = row.get("canSMILES", "")
        ds = row.get("dataset", "")
        loser_id = row.get("__loser_id__")  # might be None
        if csmiles and ds:
            if csmiles not in new_drugs_grouped:
                new_drugs_grouped[csmiles] = []
            new_drugs_grouped[csmiles].append((imp_id, ds, loser_id))

    # Map (dataset, old_improve_drug_id) -> stable_id for rewriting
    drug_id_mapping = {}

    # Assign stable IDs
    for csmiles, rowlist in new_drugs_grouped.items():
        # If we already know this canSMILES, reuse its stable_id
        if csmiles in canSMILES_to_stable_id:
            stable_id = canSMILES_to_stable_id[csmiles]
        else:
            # We need to assign a stable_id for this new canSMILES
            # Sort candidate imp_ids by dataset priority
            imp_id_info = []
            for (imp_id, ds, loser_id) in rowlist:
                priority = dataset_priority.get(ds, float('inf'))
                imp_id_info.append((priority, imp_id, ds, loser_id))
            imp_id_info.sort(key=lambda x: x[0])  # ascending priority

            stable_id = None
            for (_, imp_id, ds, loser_id) in imp_id_info:
                if not imp_id:
                    # Can't resue duplicates
                    continue
                # Conditions to reuse imp_id as stable_id:
                # It starts with "SMI_"
                # Its not already assigned
                #Its unique to this canSMILES in entire input
                if (imp_id.upper().startswith("SMI_") 
                    and imp_id not in assigned_stable_ids 
                    and len(imp_id_to_canSMILES.get(imp_id, set())) == 1):
                    stable_id = imp_id
                    assigned_stable_ids.add(stable_id)
                    break

            # If none reusable, generate a new ID
            if not stable_id:
                next_available_id += 1
                stable_id = f"SMI_{next_available_id}"
                assigned_stable_ids.add(stable_id)

            # Update internal data
            canSMILES_to_stable_id[csmiles] = stable_id
            drugs_list.append({
                "stable_id": stable_id,
                "canSMILES": csmiles
            })
            print(f"Assigned stable_id '{stable_id}' to canSMILES '{csmiles}'.")

        # For each row referencing this canSMILES, build rewriting map
        for (imp_id, ds, loser_id) in rowlist:
            # If a row had an actual old ID, and it differs from stable_id, map it
            if imp_id and imp_id != stable_id:
                drug_id_mapping[(ds, imp_id)] = stable_id
                print(f"Mapping improve_drug_id '{imp_id}' (dataset '{ds}') to stable_id '{stable_id}'.")
            # If the row was a conflict loser, we also map (ds, loser_id) => stable_id
            if loser_id and loser_id != stable_id:
                drug_id_mapping[(ds, loser_id)] = stable_id
                print(f"Loser ID '{loser_id}' (dataset '{ds}') => stable_id '{stable_id}'.")

    return drug_id_mapping

def rewrite_drugs_file(file_path, drug_id_mapping, datasets=None):
    """
    Rewrites drug files by replacing improve_drug_id with stable_id based on drug_id_mapping.
    """
    filename = os.path.basename(file_path)
    dataset = filename.split('_')[0]
    if datasets and dataset not in datasets:
        return
    # Check file existence
    if not os.path.exists(file_path):
        gz_path = file_path + '.gz'
        if os.path.exists(gz_path):
            file_path = gz_path
        else:
            print(f"File not found: {file_path} or {gz_path}")
            return
    print(f"Processing drug file: {file_path}")
    file_path_decompressed, was_gz = decompress_gz_if_needed(file_path)
    if (file_path_decompressed is None or 
        not os.path.exists(file_path_decompressed) or 
        os.path.getsize(file_path_decompressed) == 0):
        recompress_if_needed(file_path_decompressed, was_gz, original_path=file_path if was_gz else None)
        print(f"File not found or empty after decompression: {file_path_decompressed}")
        return

    changes_made = False
    tmp_path = file_path_decompressed + ".tmp"
    with open(file_path_decompressed, 'r', newline='', encoding='utf-8') as fin, \
         open(tmp_path, 'w', newline='', encoding='utf-8') as fout:
        reader = csv.DictReader(fin, delimiter='\t')
        fieldnames = reader.fieldnames
        if "improve_drug_id" not in fieldnames:
            recompress_if_needed(file_path_decompressed, was_gz, original_path=file_path if was_gz else None)
            print(f"'improve_drug_id' column not found in {file_path_decompressed}")
            return
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            original_id = row.get("improve_drug_id", "").strip()
            key = (dataset, original_id)
            if key in drug_id_mapping:
                new_id = drug_id_mapping[key]
                if original_id != new_id:
                    print(f"Replacing improve_drug_id '{original_id}' with stable_id '{new_id}' in {file_path_decompressed}")
                    row["improve_drug_id"] = new_id
                    changes_made = True
            writer.writerow(row)

    if changes_made:
        os.replace(tmp_path, file_path_decompressed)
        print(f"Changes were made to {file_path_decompressed}. File has been updated.")
    else:
        os.remove(tmp_path)
        print(f"No changes made to {file_path_decompressed}. File remains unchanged.")

    recompress_if_needed(file_path_decompressed, was_gz, original_path=file_path if was_gz else None)

def main():
    parser = argparse.ArgumentParser(description="""
Use canSMILES overlaps to assign stable drug IDs (no build_statuses).
Prioritize datasets to retain improve_drug_id values from higher-priority datasets.
Unify via canSMILES, update stable IDs, and rewrite files.

We resolve conflicts so that if a single SMI_xxx is used for multiple canSMILES,
the highest-priority dataset keeps that SMI_xxx, and all other usages get new IDs.
We also ensure the old IDs in losing rows get mapped to the new stable ID
so their TSV files are properly rewritten.
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
    parser.add_argument('--input_files', nargs='+',
                        help='List of input files to process. If specified, only these files will be processed.')
    args = parser.parse_args()

    # Set build_date
    build_date = args.build_date or datetime.utcnow().strftime("%Y-%m-%d")
    # Mapping file path
    mapping_file = 'improve_drug_mapping.json'

    # Load or initialize improve_drug_mapping.json
    mapping_data, had_prior = load_mapping(mapping_file)
    # Ensure mapping_data includes existing "drugs" list
    if not mapping_data.get('drugs'):
        mapping_data['drugs'] = []

    # Print existing stable_ids
    existing_stable_ids = [drug["stable_id"] for drug in mapping_data["drugs"]]
    print(f"Existing stable_ids loaded: {existing_stable_ids}")

    # Insert current build metadata, check uniqueness
    current_build_metadata = get_current_build_metadata(build_date, args.version)
    if any(
        b["build_date"] == build_date and b["version"] == args.version
        for b in mapping_data["metadata"]["builds"]
    ):
        print(f"Build {build_date}, version {args.version} already exists in {mapping_file}. Exiting.")
        return
    mapping_data["metadata"]["builds"].append(current_build_metadata)

    input_files = args.input_files
    if input_files:
        # Prepare dataset priorities (in the order they appear)
        dataset_priority = {}
        all_drugs_rows = []
        for file_path in input_files:
            filename = os.path.basename(file_path)
            dataset = filename.split('_')[0]
            if dataset not in dataset_priority:
                dataset_priority[dataset] = len(dataset_priority) + 1
            rows = read_drugs_tsv(file_path, dataset)
            all_drugs_rows.extend(rows)

        # 1) Resolve ID-level conflicts so highest-priority usage keeps the old ID
        all_drugs_rows = resolve_id_conflicts(all_drugs_rows, dataset_priority)

        print("Unifying drugs and updating mappings.")
        print(f"Dataset priorities: {dataset_priority}")
        drug_id_mapping = unify_drugs(mapping_data, all_drugs_rows, dataset_priority)

        # Sort final "drugs" by numeric portion if SMI_#, or lexically
        mapping_data["drugs"].sort(key=lambda d: parse_smi_num(d["stable_id"]) or 999999999)
        save_mapping(mapping_data, mapping_file)
        print(f"Updated {mapping_file} with {len(mapping_data['drugs'])} drugs.")

        # Rewrite files
        print("Rewriting files with updated stable IDs.")
        for file_path in input_files:
            rewrite_drugs_file(file_path, drug_id_mapping)
    else:
        # Prepare dataset and file type lists
        ds_list = [d.strip() for d in args.datasets.split(',')]
        other_file_types = [x.strip() for x in args.other_files.split(',')]
        local_dir = args.local_dir

        # Define dataset priorities based on the order they appear in ds_list
        dataset_priority = {ds: idx + 1 for idx, ds in enumerate(ds_list)}
        print("Dataset priorities:", dataset_priority)

        # Gather all drugs from each dataset
        all_drugs_rows = []
        for ds in ds_list:
            path_tsv = os.path.join(local_dir, f"{ds}_drugs.tsv")
            if os.path.exists(path_tsv) or os.path.exists(path_tsv + '.gz'):
                rows = read_drugs_tsv(path_tsv, ds)
                all_drugs_rows.extend(rows)
            else:
                print(f"No drug file found for dataset {ds}: {path_tsv} or {path_tsv + '.gz'}")

        # 1) Resolve ID conflicts so highest-priority usage keeps the old ID
        all_drugs_rows = resolve_id_conflicts(all_drugs_rows, dataset_priority)

        # 2) Unify drugs and get mappings
        print("Unifying drugs and updating mappings.")
        drug_id_mapping = unify_drugs(mapping_data, all_drugs_rows, dataset_priority)

        # 3) Sort final "drugs" by numeric portion  If not a SMI_number, reassign to 999999999999.
        mapping_data["drugs"].sort(key=lambda d: parse_smi_num(d["stable_id"]) or 999999999999)
        save_mapping(mapping_data, mapping_file)
        print(f"Updated {mapping_file} with {len(mapping_data['drugs'])} drugs.")

        # 4) Rewrite the main drug TSV
        print("Rewriting files with updated stable IDs.")
        for ds in ds_list:
            path_tsv = os.path.join(local_dir, f"{ds}_drugs.tsv")
            if os.path.exists(path_tsv) or os.path.exists(path_tsv + '.gz'):
                rewrite_drugs_file(path_tsv, drug_id_mapping, datasets=ds_list)
            else:
                print(f"No drug file found for dataset {ds}: {path_tsv} or {path_tsv + '.gz'}")

        # 5) Rewrite other file types (drug_descriptors, experiments, etc.)
        for ds in ds_list:
            for ftype in other_file_types:
                other_path = os.path.join(local_dir, f"{ds}_{ftype}.tsv")
                if os.path.exists(other_path) or os.path.exists(other_path + '.gz'):
                    rewrite_drugs_file(other_path, drug_id_mapping, datasets=ds_list)
                else:
                    print(f"No {ftype} file found for dataset {ds}: {other_path} or {other_path + '.gz'}")

        print("All files have been rewritten with stable IDs.")
        print(f"Stable IDs have been updated in {mapping_file}.")

if __name__ == "__main__":
    main()
