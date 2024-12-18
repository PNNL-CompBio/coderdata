import json
import os
import csv
import argparse
from datetime import datetime
from collections import defaultdict
import gzip
import shutil

def load_mapping(mapping_file='mapping.json'):
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
        
def save_mapping(data, mapping_file='mapping.json'):
    with open(mapping_file, 'w') as f:
        json.dump(data, f, indent=2)

def get_current_build_metadata(build_date, git_commit, description, version):
    return {
        "build_date": build_date,
        "git_commit": git_commit,
        "description": description,
        "version": version
    }

def read_current_samples(sample_file):
    # Check if file exists as-is
    if not os.path.exists(sample_file):
        # If not found, check if a gzipped version exists
        gz_file = sample_file + '.gz'
        if os.path.exists(gz_file):
            sample_file = gz_file
        else:
            # Neither plain nor gz file found, return empty
            return []

    is_gz = sample_file.endswith('.gz')
    file_to_read = sample_file
    if is_gz:
        decompressed_path = sample_file[:-3]
        with gzip.open(sample_file, 'rb') as f_in, open(decompressed_path, 'wb', buffering=1024*1024) as f_out:
            shutil.copyfileobj(f_in, f_out)
        file_to_read = decompressed_path

    current_samples = []
    print(f"reading {file_to_read}")
    with open(file_to_read, 'r', newline='', encoding='utf-8', buffering=2**20) as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            cancer_type = row.get("cancer_type", "").strip()
            common_name = row.get("common_name", "").strip()
            model_type = row.get("model_type", "").strip()
            temporary_id = row.get("improve_sample_id", "").strip()
            current_samples.append({
                "temporary_improve_id": temporary_id,
                "cancer_type": cancer_type,
                "common_name": common_name,
                "model_type": model_type
            })

    if is_gz:
        os.remove(file_to_read)

    return current_samples

def generate_new_id(samples):
    # Generate a new stable ID.
    return str(len(samples) + 1)

def update_sample_status(mapping_data, all_current_samples, build_date):
    """
    Update sample statuses and unify ephemeral IDs. Returns a dictionary mapping any new ephemeral IDs
    to stable IDs for known triples.
    """
    def triple_key(c_type, c_name, m_type):
        return (c_type, c_name, m_type)

    samples = mapping_data["samples"]
    existing_triples = {
        (s["cancer_type"], s.get("common_name",""), s.get("model_type","")): s
        for s in samples
    }

    # Build a lookup of improve_sample_id -> triple for conflict checking
    id_to_triple = {}
    for s in samples:
        t = (s["cancer_type"], s.get("common_name",""), s.get("model_type",""))
        id_to_triple[s["improve_sample_id"]] = t

    current_triples = [(cs["cancer_type"], cs["common_name"], cs["model_type"]) for cs in all_current_samples]
    current_triples_set = set(current_triples)

    # Sort builds once
    build_list = mapping_data["metadata"]["builds"]
    build_list_sorted = sorted(build_list, key=lambda b: b["build_date"])

    # Dictionary to store ephemeral_id -> stable_id mappings discovered for known triples
    ephemeral_id_to_stable_map = {}

    # Assign improve_sample_ids to current samples
    for cs in all_current_samples:
        t = (cs["cancer_type"], cs["common_name"], cs["model_type"])
        temporary_id = cs["temporary_improve_id"]
        if t in existing_triples:
            # Known triple, just mark present
            smpl = existing_triples[t]
            smpl["build_statuses"].append({"build_date": build_date, "status": "present"})
            
            # If there's a temporary_id that differs from the known stable_id, record a mapping
            stable_id = smpl["improve_sample_id"]
            if temporary_id and temporary_id != stable_id:
                # This means another dataset is using a different ephemeral ID for the same triple
                ephemeral_id_to_stable_map[temporary_id] = stable_id

        else:
            # New triple
            if temporary_id and temporary_id in id_to_triple:
                # temporary_id already in use
                existing_triple_for_id = id_to_triple[temporary_id]
                if existing_triple_for_id != t:
                    # Conflict: this ID is used for a different triple, must generate a new stable ID
                    new_id = generate_new_id(samples)
                else:
                    # Same triple scenario shouldn't happen as triple is new, but let's keep safe
                    new_id = temporary_id
            elif temporary_id and temporary_id not in id_to_triple:
                # temporary_id is free to use
                new_id = temporary_id
            else:
                # No temporary_id or empty, generate a new stable ID
                new_id = generate_new_id(samples)

            new_sample = {
                "improve_sample_id": new_id,
                "common_name": t[1],
                "cancer_type": t[0],
                "model_type": t[2],
                "build_statuses": []
            }
            samples.append(new_sample)
            existing_triples[t] = new_sample
            id_to_triple[new_id] = t

            # Add absent for all previous builds
            for bd in build_list_sorted:
                if bd["build_date"] < build_date:
                    new_sample["build_statuses"].append({
                        "build_date": bd["build_date"],
                        "status": "absent"
                    })
            # Now add present for current build
            new_sample["build_statuses"].append({
                "build_date": build_date,
                "status": "present"
            })

    # Mark samples that appear in previous builds but not this one as absent
    for smpl in samples:
        t = (smpl["cancer_type"], smpl.get("common_name",""), smpl.get("model_type",""))
        if t not in current_triples_set:
            # If last known status was present, now absent
            bs = smpl["build_statuses"]
            if bs and bs[-1]["status"] == "present":
                bs.append({
                    "build_date": build_date,
                    "status": "absent"
                })

    return ephemeral_id_to_stable_map

def build_triple_to_improve_id_map(mapping_data):
    triple_to_improve_id = {}
    for s in mapping_data["samples"]:
        t = (s["cancer_type"], s.get("common_name",""), s.get("model_type",""))
        triple_to_improve_id[t] = s["improve_sample_id"]
    return triple_to_improve_id

def rewrite_all_files(local_dir='local', temporary_samples=None, mapping_data=None, datasets=None, extra_ephemeral_map=None):
    triple_to_improve_id = build_triple_to_improve_id_map(mapping_data)

    # Map temporary_id to stable_id from all_current_samples
    temporary_id_to_stable_id = {}
    for es in temporary_samples:
        t = (es["cancer_type"], es["common_name"], es["model_type"])
        stable_id = triple_to_improve_id.get(t)
        if stable_id and es["temporary_improve_id"]:
            temporary_id_to_stable_id[es["temporary_improve_id"]] = stable_id

    # Merge in extra_ephemeral_map discovered from known triples
    if extra_ephemeral_map:
        for eph_id, st_id in extra_ephemeral_map.items():
            temporary_id_to_stable_id[eph_id] = st_id

    print(f"temporary_id_to_stable_id:{temporary_id_to_stable_id}")
    print(f"triple_to_improve_id:{triple_to_improve_id}")

    valid_extensions = ['.csv', '.tsv', '.csv.gz', '.tsv.gz']

    def detect_delimiter(filename):
        if '.tsv' in filename:
            return '\t'
        return ','

    def maybe_decompress(file_path):
        if file_path.endswith('.gz'):
            decompressed_path = file_path[:-3]
            try:
                with gzip.open(file_path, 'rb') as f_in, open(decompressed_path, 'wb', buffering=1024*1024) as f_out:
                    shutil.copyfileobj(f_in, f_out)
                return decompressed_path, True
            except gzip.BadGzipFile:
                print(f"Warning: {file_path} is not a valid gzip file (CRC error or corruption). Skipping decompression.")
                return file_path, False
            except Exception as e:
                print(f"Error decompressing {file_path}: {e}")
                return file_path, False
        return file_path, False

    def maybe_recompress(original_path, decompressed_path, was_gz):
        if was_gz:
            with open(decompressed_path, 'rb') as f_in, gzip.open(original_path, 'wb', compresslevel=5) as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(decompressed_path)
            
    for root, dirs, files in os.walk(local_dir):
        for name in files:
            filepath = os.path.join(root, name)

            # Check extension
            if not any(filepath.endswith(ext) for ext in valid_extensions):
                continue

            # Filter by datasets: Only process if it matches one of the dataset prefixes
            if datasets and not any(ds in name for ds in datasets):
                # Skip files that don't match known datasets
                continue
            print(f"File to rewrite: {filepath}")
            decompressed_path, was_gz = maybe_decompress(filepath)
            delimiter = detect_delimiter(decompressed_path)

            # Read just the header
            with open(decompressed_path, 'r', newline='', encoding='utf-8', buffering=2**20) as f:
                reader = csv.reader(f, delimiter=delimiter)
                try:
                    header = next(reader)
                except StopIteration:
                    maybe_recompress(filepath, decompressed_path, was_gz)
                    continue

            if "improve_sample_id" not in header:
                maybe_recompress(filepath, decompressed_path, was_gz)
                continue

            idx = header.index("improve_sample_id")

            temp_out_path = decompressed_path + ".tmp"
            with open(decompressed_path, 'r', newline='', encoding='utf-8', buffering=2**20) as f_in, \
                 open(temp_out_path, 'w', newline='', encoding='utf-8', buffering=2**20) as f_out:
                reader = csv.reader(f_in, delimiter=delimiter)
                writer = csv.writer(f_out, delimiter=delimiter)
                writer.writerow(header)  # Write header
                next(reader, None)  # skip header row in input

                for row in reader:
                    if len(row) <= idx:
                        writer.writerow(row)
                        continue
                    old_id = row[idx].strip()
                    if old_id in temporary_id_to_stable_id:
                        row[idx] = temporary_id_to_stable_id[old_id]
                    writer.writerow(row)

            os.replace(temp_out_path, decompressed_path)
            maybe_recompress(filepath, decompressed_path, was_gz)

def main():
    parser = argparse.ArgumentParser(description="Update coderdata sample mapping and rewrite files.")
    parser.add_argument('--build_date', type=str, help='Build date in YYYY-MM-DD format. If not provided, uses current UTC date.')
    parser.add_argument('--git_commit', type=str, default='abcdef12345', help='Git commit hash.')
    parser.add_argument('--description', type=str, default='New build with updated samples.', help='Build description.')
    parser.add_argument('--version', type=str, default='0.1.41', help='Version string.')
    parser.add_argument('--datasets', type=str, default='fimm,mpnst,ccle,hcmi,beataml,cptac', help='Comma-separated list of datasets.')
    parser.add_argument('--local_dir', type=str, default='local', help='Local directory containing files to rewrite.')

    args = parser.parse_args()

    build_date = args.build_date if args.build_date else datetime.utcnow().strftime("%Y-%m-%d")
    git_commit = args.git_commit
    description = args.description
    version = args.version
    datasets = [d.strip() for d in args.datasets.split(',')]
    local_dir = args.local_dir

    mapping_data, previous_mapping_found = load_mapping()

    current_build_meta = get_current_build_metadata(build_date, git_commit, description, version)
    if not any(b["build_date"] == build_date and b["version"] == version for b in mapping_data["metadata"]["builds"]):
        mapping_data["metadata"]["builds"].append(current_build_meta)

    all_current_samples = []
    for ds in datasets:
        print(f"Dataset in progress: {ds}")
        sample_file = f"{local_dir}/{ds}_samples.csv"
        cs = read_current_samples(sample_file)
        all_current_samples.extend(cs)
        
    # Deduplicate all_current_samples so that each (cancer_type, common_name, model_type) appears only once
    unique_samples = {}
    for cs in all_current_samples:
        triple = (cs["cancer_type"], cs["common_name"], cs["model_type"])
        if triple not in unique_samples:
            unique_samples[triple] = cs
    all_current_samples = list(unique_samples.values())

    print(f"All samples: {all_current_samples}")
    print(f"Updating sample status")
    # Update sample status and get ephemeral_id_to_stable_map for known triples
    extra_ephemeral_map = update_sample_status(mapping_data, all_current_samples, build_date)
    print(f"Saving mapping")
    save_mapping(mapping_data)
    print(f"Rewriting all files")    
    rewrite_all_files(local_dir=local_dir, temporary_samples=all_current_samples, mapping_data=mapping_data, datasets=datasets, extra_ephemeral_map=extra_ephemeral_map)
    print("Mapping updated and all files rewritten with stable improve_sample_ids.")

if __name__ == "__main__":
    main()
