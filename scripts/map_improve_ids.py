import json
import os
import csv
# NOTE: Will add argparse later
from datetime import datetime
from collections import defaultdict
import gzip
import shutil
from glob import glob

def load_mapping(mapping_file='mapping.json'):
    # Just returning a minimal structure for now
    # handle JSON schema 
    if os.path.exists(mapping_file):
        with open(mapping_file, 'r') as f:
            return json.load(f)
    else:
        return {
            "metadata": {
                "builds": []
            },
            "samples": []
        }

def save_mapping(data, mapping_file='mapping.json'):
    # Maybe add backups later
    with open(mapping_file, 'w') as f:
        json.dump(data, f, indent=2)

def get_current_build_metadata(build_date, git_commit, description, version):
    # maybe validate or something
    return {
        "build_date": build_date,
        "git_commit": git_commit,
        "description": description,
        "version": version
    }

def read_current_samples(sample_file):
    # Just handle CSV/TSV and gz. Maybe unify later.
    if not os.path.exists(sample_file) or os.path.getsize(sample_file) == 0:
        return []

    _, ext = os.path.splitext(sample_file)
    delimiter = ',' if '.csv' in ext else '\t'

    is_gz = sample_file.endswith('.gz')
    if is_gz:
        # decompress, maybe store in tmp
        decompressed_path = sample_file[:-3]
        with gzip.open(sample_file, 'rb') as f_in, open(decompressed_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        file_to_read = decompressed_path
    else:
        file_to_read = sample_file

    current_samples = []
    with open(file_to_read, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            # maybe check fields some day
            temporary_id = row.get("improve_sample_id", "").strip()
            current_samples.append({
                "temporary_improve_id": temporary_id,
                "other_id": row.get("other_id", "").strip(),
                "common_name": row.get("common_name", "").strip(),
                "model_type": row.get("model_type", "").strip()
            })

    if is_gz:
        os.remove(file_to_read)

    return current_samples

def update_sample_status(mapping_data, all_current_samples, build_date):
    # might need optimization if huge
    def triple_key(s):
        return (s.get("other_id", ""), s.get("common_name", ""), s.get("model_type", ""))

    existing_triples = {
        (smpl["other_id"], smpl.get("common_name",""), smpl.get("model_type","")): smpl
        for smpl in mapping_data["samples"]
    }

    current_triples_set = {
        (cs["other_id"], cs["common_name"], cs["model_type"]) for cs in all_current_samples
    }

    build_list = mapping_data["metadata"]["builds"]
    build_list_sorted = sorted(build_list, key=lambda b: b["build_date"])

    for cs in all_current_samples:
        t = (cs["other_id"], cs["common_name"], cs["model_type"])
        if t in existing_triples:
            smpl = existing_triples[t]
            smpl["build_statuses"].append({"build_date": build_date, "status": "present"})
        else:
            # in progress
            new_id = str(len(mapping_data["samples"]) + 1)
            new_sample = {
                "improve_sample_id": new_id,
                "common_name": t[1],
                "other_id": t[0],
                "model_type": t[2],
                "build_statuses": []
            }
            mapping_data["samples"].append(new_sample)
            existing_triples[t] = new_sample

            # mark old builds as absent
            for bd in build_list_sorted:
                if bd["build_date"] < build_date:
                    new_sample["build_statuses"].append({
                        "build_date": bd["build_date"],
                        "status": "absent"
                    })
            new_sample["build_statuses"].append({
                "build_date": build_date,
                "status": "present"
            })

    # mark absent now
    for smpl in mapping_data["samples"]:
        t = (smpl["other_id"], smpl.get("common_name",""), smpl.get("model_type",""))
        if t not in current_triples_set:
            if smpl.get("build_statuses") and smpl["build_statuses"][-1]["status"] == "present":
                smpl["build_statuses"].append({
                    "build_date": build_date,
                    "status": "absent"
                })

def build_triple_to_improve_id_map(mapping_data):
    # Maybe store this elsewhere
    triple_to_improve_id = {}
    for s in mapping_data["samples"]:
        t = (s["other_id"], s.get("common_name",""), s.get("model_type",""))
        triple_to_improve_id[t] = s["improve_sample_id"]
    return triple_to_improve_id

def rewrite_all_files(local_dir='local', temporary_samples=None, mapping_data=None):
    # handle files later more carefully
    triple_to_improve_id = build_triple_to_improve_id_map(mapping_data)

    temporary_id_to_stable_id = {}
    for es in temporary_samples:
        t = (es["other_id"], es["common_name"], es["model_type"])
        stable_id = triple_to_improve_id.get(t)
        if stable_id:
            temporary_id_to_stable_id[es["temporary_improve_id"]] = stable_id

    def detect_delimiter(filename):
        # might guess from file content someday
        if '.tsv' in filename:
            return '\t'
        return ','

    def maybe_decompress(file_path):
        # handle memory-only in future
        if file_path.endswith('.gz'):
            decompressed_path = file_path[:-3]
            with gzip.open(file_path, 'rb') as f_in, open(decompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            return decompressed_path, True
        return file_path, False

    def maybe_recompress(original_path, decompressed_path, was_gz):
        if was_gz:
            with open(decompressed_path, 'rb') as f_in, gzip.open(original_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(decompressed_path)

    for root, dirs, files in os.walk(local_dir):
        for name in files:
            filepath = os.path.join(root, name)
            if any(filepath.endswith(ext) for ext in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']):
                decompressed_path, was_gz = maybe_decompress(filepath)
                delimiter = detect_delimiter(decompressed_path)

                with open(decompressed_path, 'r', newline='', encoding='utf-8') as f:
                    reader = csv.reader(f, delimiter=delimiter)
                    rows = list(reader)

                if not rows:
                    maybe_recompress(filepath, decompressed_path, was_gz)
                    continue

                header = rows[0]
                if "improve_sample_id" not in header:
                    maybe_recompress(filepath, decompressed_path, was_gz)
                    continue

                idx = header.index("improve_sample_id")

                for i in range(1, len(rows)):
                    old_id = rows[i][idx].strip()
                    if old_id in temporary_id_to_stable_id:
                        rows[i][idx] = temporary_id_to_stable_id[old_id]

                with open(decompressed_path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f, delimiter=delimiter)
                    writer.writerows(rows)

                maybe_recompress(filepath, decompressed_path, was_gz)

def main():
    # NOTE: add argparse here 
    # For now just hardcode
    build_date = datetime.utcnow().strftime("%Y-%m-%d")
    git_commit = 'abcdef12345'
    description = 'New build with updated samples.'
    version = '0.1.41'
    datasets = ['mpnst','ccle','hcmi','beataml','cptac']
    local_dir = 'local'

    mapping_data = load_mapping()
    current_build_meta = get_current_build_metadata(build_date, git_commit, description, version)
    if not any(b["build_date"] == build_date and b["version"] == version for b in mapping_data["metadata"]["builds"]):
        mapping_data["metadata"]["builds"].append(current_build_meta)

    # just loop over datasets
    all_current_samples = []
    for ds in datasets:
        # change to tmp when run in docker
        sample_file = f"local/{ds}_samples.csv"
        cs = read_current_samples(sample_file)
        all_current_samples.extend(cs)

    update_sample_status(mapping_data, all_current_samples, build_date)
    save_mapping(mapping_data)
    rewrite_all_files(local_dir=local_dir, temporary_samples=all_current_samples, mapping_data=mapping_data)

    print("Mapping updated. Incomplete and still sloppy.")

if __name__ == "__main__":
    main()
