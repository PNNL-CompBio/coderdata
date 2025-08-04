#!/usr/bin/env python3
import synapseclient
import pandas as pd
import numpy as np
import argparse
import os
import re
import subprocess

# Helper functions
def _clean_geo_id(s):
    """
    Normalise GEO sample IDs so they match Synapse naming.
      • 11.2  → 11_2
      • **_Tumor → *_Parental
      • *_orgP2 → *_Organoid_P2
      • *_xenoorgP4 → *_XenoOrganoid_P4
    """
    s = s.strip()
    s = re.sub(r"(?<=\d)\.(?=\d)", "_", s)          # dots between digits
    s = s.replace("_tumor", "_Parental")            # tumour alias
    # lower-case 'orgP' / 'xenoorgP' fix
    s = re.sub(r"_(org)P(\d+)",      r"_Organoid_P\2",      s, flags=re.IGNORECASE)
    s = re.sub(r"_(xenoorg)P(\d+)",  r"_XenoOrganoid_P\2",  s, flags=re.IGNORECASE)
    return s


def _parse_model_type(sample_id):
    """Derive model_type from Sample ID."""
    low = sample_id.lower()
    if "_xenoorganoid" in low:
        return "xenograft derived organoid"
    if "_organoid" in low:
        return "patient derived organoid"
    if "_xenograft" in low:
        return "patient derived xenograft"
    if "_parental" in low:
        return "tumor"
    return "unknown"

#Generate Samples Data
def get_bladder_pdo_samples(synLoginObject, maxval):
    
    
    #Part 1: Get Data from Synapse
    
    # download from Synapse..
    samples_syn = synLoginObject.get('syn64765486')
    # and read the file
    samples_df = pd.read_csv(samples_syn.path, sep="\t")

    samples = samples_df[['Sample ID', 'Patient ID', 'Cancer Type Detailed', 'Sample Class']]
    samples = samples.rename({"Sample ID" : 'other_id', 'Patient ID' : 'common_name', 'Cancer Type Detailed': 'cancer_type', 'Sample Class' : 'model_type'}, axis=1)

    samples.loc[:,['species']] = 'Homo sapiens(Human)'
    samples.loc[:,['other_id_source']] = 'Synapse'
    samples.loc[:,['other_names'] ]= ''
    samples.loc[:,['cancer_type']]=samples['cancer_type'].str.lower()
    samples["model_type"] = samples["other_id"].apply(_parse_model_type)

    #Part 2: Get Data from Geo
    subprocess.call (["Rscript", "--vanilla", "obtainGSMidLink.R"])
    GEO_ids_link = "./gsmlinkDf.csv"

    geo_map  = pd.read_csv(GEO_ids_link)
    geo_ids  = geo_map["sampleid"].dropna().map(_clean_geo_id).unique()
    missing  = sorted(set(geo_ids) - set(samples["other_id"]))

    if missing:                    
        print(f"Adding {len(missing)} GEO samples not in Synapse sheet")

    rows = []                      
    for oid in missing:            
        common = oid.split("_")[0]
        ctype  = (
            samples.loc[samples["common_name"] == common, "cancer_type"]
            .iloc[0]
            if (samples["common_name"] == common).any()
            else "bladder urothelial carcinoma"
        )
        rows.append(
            {
                "other_id":        oid,
                "common_name":     common,
                "cancer_type":     ctype,
                "model_type":      _parse_model_type(oid),
                "species":         "Homo sapiens(Human)",
                "other_id_source": "GEO",
                "other_names":     "",
            }
        )
    if rows:                       
        samples = pd.concat([samples, pd.DataFrame(rows)], ignore_index=True)

    samples = samples.sort_values("other_id").reset_index(drop=True)

    samples['improve_sample_id'] = range(maxval+1, maxval+1+samples.shape[0])

    return samples


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="This script handles downloading, processing and formatting of sample files for the Sarcoma PDO project into a single samplesheet")
    parser.add_argument('-t', '--token', type=str, help='Synapse Token')
    parser.add_argument("-p", '--prevSamples', nargs="?", type=str, default ="", const  = "", help = "Use this to provide previous sample file, will run sample file generation")
    args = parser.parse_args()
   
    print("Logging into Synapse")
    PAT = args.token
    synObject = synapseclient.login(authToken=PAT)

    if (args.prevSamples):
        prev_max_improve_id = max(pd.read_csv(args.prevSamples).improve_sample_id)
    else: 
        prev_max_improve_id = 0

    bladder_pdo_samples = get_bladder_pdo_samples(synObject, prev_max_improve_id)
    bladder_pdo_samples.to_csv("/tmp/bladderpdo_samples.csv", index=False)