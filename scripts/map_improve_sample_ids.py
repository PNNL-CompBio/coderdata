import json
import os
import csv
import argparse
from datetime import datetime
import gzip
import shutil
from collections import defaultdict

##### Load, make structure, save functions for improve_sample_mapping.json
def load_mapping(mapping_file='improve_sample_mapping.json'):
    """
    Loads an existing improve_sample_mapping.json if available.
    Otherwise returns an empty base structure with metadata and samples.
    """
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            return json.load(f), True
    else:
        return {
            "metadata": {
                "builds": []
            },
            "samples": []
        }, False

def save_mapping(mapping_data):
    """Saves mapping data to disk as JSON."""
    mapping_file='tmp/improve_sample_mapping.json'
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

##### Read {dataset}_samples.csv files

def read_samples_csv(file_path, dataset):
    """
    Reads a {dataset}_samples.csv or .csv.gz and returns a list of row dicts,
    adding the dataset name to each row.
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
    print(f"Reading samples file: {file_path}")
    with open(file_path,'r',newline='',encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            row['dataset'] = dataset  # Add the dataset name to each row
            rows.append(row)
    recompress_if_needed(file_path, file_path, was_gz)
    return rows

##### Create New stable ID (improve_sample_id)

def generate_new_stable_id(assigned_stable_ids):
    """Generate new stable ID by finding the max existing stable_id and adding 1."""
    existing_ids = [int(sid) for sid in assigned_stable_ids if sid.isdigit()]
    if existing_ids:
        return str(max(existing_ids) + 1)
    else:
        return "1"

##### Assign stable ID based on improve_sample_id, quadruplet, and previous improve_sample_mapping.json file

def unify_samples(mapping_data, all_samples_rows, current_build_metadata):
    """
    Assign stable IDs to samples based on quadruplet overlaps.
    Parameters:
        mapping_data: existing mapping data loaded from improve_sample_mapping.json
        all_samples_rows: list of row dicts from the new build's samples CSVs
        current_build_metadata: dict containing build_date, version, etc.
    Returns:
        quadruplet_to_stable_id: dict mapping quadruplet tuples to stable_id
        sample_id_mapping: dict mapping new build's improve_sample_id (with dataset) to stable_id
    """
    samples_list = mapping_data["samples"]
    # List of builds sorted by date and version
    builds = sorted(mapping_data["metadata"]["builds"], key=lambda x: (x["build_date"], x["version"]))
    # Build quadruplet_to_stable_id mapping from existing samples
    quadruplet_to_stable_id = {}
    for sample in samples_list:
        stable_id = sample["stable_id"]
        for quadruplet in sample.get("quadruples", []):
            quad = (quadruplet["other_id"], quadruplet["other_id_source"], quadruplet["model_type"], quadruplet.get("other_names", "").strip())
            quadruplet_to_stable_id[quad] = stable_id
    # Initialize assigned stable_ids with existing stable_ids
    assigned_stable_ids = set(sample["stable_id"] for sample in samples_list)
    # Group new samples by (dataset, improve_sample_id)
    new_samples_grouped = defaultdict(set)  # (dataset, improve_sample_id) -> set of quadruplet tuples
    for row in all_samples_rows:
        improve_sample_id = row.get("improve_sample_id", "").strip()
        dataset = row.get("dataset", "").strip()
        if not improve_sample_id:
            continue  # skip rows without improve_sample_id
        quadruplet = (
            row.get("other_id", "").strip(),
            row.get("other_id_source", "").strip(),
            row.get("model_type", "").strip(),
            row.get("other_names", "").strip()
        )
        new_samples_grouped[(dataset, improve_sample_id)].add(quadruplet)
    # Mapping from new build's (dataset, improve_sample_id) to stable_id
    sample_id_mapping = {}
    # Set to track which stable_ids are present in current build
    present_stable_ids = set()
    # Process each improve_sample_id
    for (dataset, improve_id), quadruplets in new_samples_grouped.items():
        # Check if any quadruplet overlaps with existing samples
        matched_stable_id = None
        for quadruplet in quadruplets:
            if quadruplet in quadruplet_to_stable_id:
                matched_stable_id = quadruplet_to_stable_id[quadruplet]
                break  # assume one match is enough
        if matched_stable_id:
            stable_id = matched_stable_id
            sample_id_mapping[(dataset, improve_id)] = stable_id
            present_stable_ids.add(stable_id)
            # Update the sample's quadruplet list with any new quadruplets
            for sample in samples_list:
                if sample["stable_id"] == stable_id:
                    existing_quadruplets = set(
                        (q["other_id"], q["other_id_source"], q["model_type"], q.get("other_names", "").strip())
                        for q in sample.get("quadruples", [])
                    )
                    new_unique_quadruplets = quadruplets - existing_quadruplets
                    for quadruplet in new_unique_quadruplets:
                        sample["quadruples"].append({
                            "other_id": quadruplet[0],
                            "other_id_source": quadruplet[1],
                            "model_type": quadruplet[2],
                            "other_names": quadruplet[3]
                        })
                        quadruplet_to_stable_id[quadruplet] = stable_id  # update quadruplet mapping
                    break
        else:
            # No overlapping quadruplet found
            # Try to assign improve_sample_id as stable_id if available
            if improve_id not in assigned_stable_ids:
                stable_id = improve_id
            else:
                # Assign a new unique stable_id
                stable_id = generate_new_stable_id(assigned_stable_ids)
            sample_id_mapping[(dataset, improve_id)] = stable_id
            present_stable_ids.add(stable_id)
            assigned_stable_ids.add(stable_id)
            # Create a new sample record
            new_sample = {
                "stable_id": stable_id,
                "quadruples": [
                    {"other_id": quadruplet[0], "other_id_source": quadruplet[1], "model_type": quadruplet[2], "other_names": quadruplet[3]}
                    for quadruplet in quadruplets
                ],
                "build_statuses": []
            }
            # Mark "absent" for all previous builds
            for build in builds:
                if not (build["build_date"] == current_build_metadata["build_date"] and
                        build["version"] == current_build_metadata["version"]):
                    new_sample["build_statuses"].append({
                        "build_date": build["build_date"],
                        "version": build["version"],
                        "status": "absent"
                    })
            # Mark "present" for current build
            new_sample["build_statuses"].append({
                "build_date": current_build_metadata["build_date"],
                "version": current_build_metadata["version"],
                "status": "present"
            })
            samples_list.append(new_sample)
            print(f"Initialized new sample '{stable_id}' with status 'present' for build '{current_build_metadata['version']}'.")
            # Update quadruplet_to_stable_id
            for quadruplet in quadruplets:
                quadruplet_to_stable_id[quadruplet] = stable_id
    # Now, iterate over all existing samples to update build_statuses
    for sample in samples_list:
        stable_id = sample["stable_id"]
        # Determine status for current build
        if stable_id in present_stable_ids:
            status = "present"
        else:
            status = "absent"
        # Append the status for current build if not already present
        # Check if a status for the current build (build_date + version) already exists
        already_present = any(
            (bs.get("build_date") == current_build_metadata["build_date"] and
             bs.get("version") == current_build_metadata["version"])
            for bs in sample["build_statuses"]
        )
        if not already_present:
            sample["build_statuses"].append({
                "build_date": current_build_metadata["build_date"],
                "version": current_build_metadata["version"],
                "status": status
            })
            # Uncomment the following line if you want to see the update messages
            # print(f"Appended status '{status}' for sample '{stable_id}' in build '{current_build_metadata['version']}'.")
    return quadruplet_to_stable_id, sample_id_mapping

#### Rewrite samples files based on newly assigned (Stable) IDs

def rewrite_samples_file(file_path, sample_id_mapping, datasets=None, dataset=None):
    """
    Rewrites {dataset}_samples.csv or {dataset}_samples.csv.gz:
      - Replace improve_sample_id with stable_id based on sample_id_mapping
      - If the file is part of the specified datasets
    """
    if datasets and not any(ds in os.path.basename(file_path) for ds in datasets):
        return
    if dataset is None:
        print(f"Dataset not specified for {file_path}. Skipping.")
        return
    print(f"Rewriting samples file: {file_path}")
    file_path, was_gz = decompress_gz_if_needed(file_path)
    if not os.path.exists(file_path):
        recompress_if_needed(file_path, file_path, was_gz)
        print(f"File not found or empty after decompression: {file_path}")
        return
    with open(file_path,'r',newline='',encoding='utf-8') as f:
        reader = csv.reader(f)
        rows = list(reader)
        if not rows:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"Empty file: {file_path}")
            return
        header = rows[0]
        if "improve_sample_id" not in header:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"'improve_sample_id' column not found in {file_path}")
            return
        try:
            idx_id = header.index("improve_sample_id")
        except ValueError:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"'improve_sample_id' column index error in {file_path}")
            return
    tmp = file_path + ".tmp"
    with open(file_path,'r',newline='',encoding='utf-8') as fin, \
         open(tmp,'w',newline='',encoding='utf-8') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        hdr = next(reader)
        writer.writerow(hdr)
        for row in reader:
            if len(row) <= idx_id:
                writer.writerow(row)
                continue
            original_id = row[idx_id].strip()
            mapping_key = (dataset, original_id)
            if mapping_key in sample_id_mapping:
                new_id = sample_id_mapping[mapping_key]
                if new_id != original_id:
                    print(f"Replacing improve_sample_id '{original_id}' with stable_id '{new_id}' in {file_path}")
                row[idx_id] = new_id
            writer.writerow(row)
    os.replace(tmp, file_path)
    recompress_if_needed(file_path, file_path, was_gz)

#### Rewrite all other files based on Stable IDs.

def rewrite_other_file(file_path, sample_id_mapping, datasets=None, dataset=None):
    """
    Rewrites other files (e.g., transcriptomics, proteomics):
      - Replace improve_sample_id with stable_id based on sample_id_mapping
      - If the file is part of the specified datasets
    """
    if datasets and not any(ds in os.path.basename(file_path) for ds in datasets):
        return
    if dataset is None:
        print(f"Dataset not specified for {file_path}. Skipping.")
        return
    print(f"Rewriting other file: {file_path}")
    file_path, was_gz = decompress_gz_if_needed(file_path)
    if not os.path.exists(file_path):
        recompress_if_needed(file_path, file_path, was_gz)
        print(f"File not found or empty after decompression: {file_path}")
        return
    with open(file_path,'r',newline='',encoding='utf-8') as f:
        reader = csv.reader(f)
        rows = list(reader)
        if not rows:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"Empty file: {file_path}")
            return
        header = rows[0]
        if "improve_sample_id" not in header:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"'improve_sample_id' column not found in {file_path}")
            return
        try:
            idx_id = header.index("improve_sample_id")
        except ValueError:
            recompress_if_needed(file_path, file_path, was_gz)
            print(f"'improve_sample_id' column index error in {file_path}")
            return
    tmp = file_path + ".tmp"
    with open(file_path,'r',newline='',encoding='utf-8') as fin, \
         open(tmp,'w',newline='',encoding='utf-8') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        hdr = next(reader)
        writer.writerow(hdr)
        for row in reader:
            if len(row) <= idx_id:
                writer.writerow(row)
                continue
            original_id = row[idx_id].strip()
            mapping_key = (dataset, original_id)
            if mapping_key in sample_id_mapping:
                new_id = sample_id_mapping[mapping_key]
                if new_id != original_id:
                    print(f"Replacing improve_sample_id '{original_id}' with stable_id '{new_id}' in {file_path}")
                row[idx_id] = new_id
            writer.writerow(row)
    os.replace(tmp, file_path)
    recompress_if_needed(file_path, file_path, was_gz)

#### Call everything in Main

def main():
    parser = argparse.ArgumentParser(description="""
Use quadruplet overlaps to assign stable IDs across builds.
Generate improve_sample_mapping.json.
Match samples via quadruplet overlaps, assign stable IDs,
update improve_sample_mapping.json, and rewrite files by replacing improve_sample_id with stable_id.
""")
    parser.add_argument('--build_date', default=None,
                        help='Build date in YYYY-MM-DD. Default=now.')
    parser.add_argument('--version', required=True,
                        help='Build version. Must be unique per build.')
    parser.add_argument('--datasets', default='ccle,ctrpv2,fimm,gcsi,gdscv1,gdscv2,nci60,prism,hcmi,beataml,cptac,mpnst,mpnstpdx',
                        help='Comma-separated list of datasets, e.g., beataml,ccle')
    parser.add_argument('--local_dir', default='data',
                        help='Directory containing all CSV/TSV files.')
    parser.add_argument('--other_files', default='transcriptomics,proteomics,mutations,copy_number,experiments',
                        help='Comma-separated list of other file types to rewrite, e.g., transcriptomics,proteomics')
    args = parser.parse_args()
    # Set build_date
    build_date = args.build_date or datetime.utcnow().strftime("%Y-%m-%d")
    # Load or initialize improve_sample_mapping.json
    mapping_file = "improve_sample_mapping.json"
    mapping_data, had_prior = load_mapping(mapping_file)
    # Insert current build metadata
    current_build_metadata = get_current_build_metadata(build_date, args.version)
    # Ensure that each build has a unique combination of build_date and version
    if not any(b["build_date"] == build_date and b["version"] == args.version for b in mapping_data["metadata"]["builds"]):
        mapping_data["metadata"]["builds"].append(current_build_metadata)
    else:
        print(f"Build with date {build_date} and version {args.version} already exists in improve_sample_mapping.json.")
        return
    # Prepare dataset and file type lists
    ds_list = [d.strip() for d in args.datasets.split(',')]
    other_file_types = [x.strip() for x in args.other_files.split(',')]
    local_dir = args.local_dir
    # Gather all sample rows from each dataset
    all_samples_rows = []
    for ds in ds_list:
        samp_path = os.path.join(local_dir, f"{ds}_samples.csv")
        rows = read_samples_csv(samp_path, ds)
        all_samples_rows.extend(rows)
    # Process samples to assign stable IDs
    print("Processing samples to assign stable IDs using quadruplets.")
    # Unify samples and get mappings
    quadruplet_to_stable_id, sample_id_mapping = unify_samples(mapping_data, all_samples_rows, current_build_metadata)
    # Sort samples by stable_id ascending
    mapping_data["samples"] = sorted(
        mapping_data["samples"],
        key=lambda x: int(x["stable_id"]) if x["stable_id"].isdigit() else x["stable_id"]
    )
    # Save improve_sample_mapping.json
    save_mapping(mapping_data)
    print(f"improve_sample_mapping.json updated with {len(mapping_data['samples'])} samples.")
    # Determine if any improve_sample_ids were modified
    ids_modified = any(
        (mapping_key[1] != sample_id_mapping[mapping_key])
        for mapping_key in sample_id_mapping
    )
    if not ids_modified:
        if not had_prior:
            print("First build detected. No IDs modified. No file rewriting needed.")
            return
        else:
            print("No improve_sample_id modifications detected. No file rewriting needed.")
            return
    # Proceed to rewrite files
    print("Rewriting files with updated stable IDs.")
    # Rewrite samples files
    for ds in ds_list:
        samp_path = os.path.join(local_dir, f"{ds}_samples.csv")
        rewrite_samples_file(samp_path, sample_id_mapping, datasets=ds_list, dataset=ds)
        samp_gz_path = samp_path + '.gz'
        if os.path.exists(samp_gz_path):
            rewrite_samples_file(samp_gz_path, sample_id_mapping, datasets=ds_list, dataset=ds)
    # Rewrite other files
    for ds in ds_list:
        for ftype in other_file_types:
            # Possible file extensions
            for ext in ['.csv', '.tsv']:
                other_file = os.path.join(local_dir, f"{ds}_{ftype}{ext}")
                if os.path.exists(other_file):
                    rewrite_other_file(other_file, sample_id_mapping, datasets=ds_list, dataset=ds)
                gz_other_file = other_file + '.gz'
                if os.path.exists(gz_other_file):
                    rewrite_other_file(gz_other_file, sample_id_mapping, datasets=ds_list, dataset=ds)
    print("All files have been rewritten with stable IDs.")
    print("Stable IDs have been updated in improve_sample_mapping.json.")

if __name__ == "__main__":
    main()