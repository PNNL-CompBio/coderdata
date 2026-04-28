# cNF — Cutaneous Neurofibroma Organoid Drug Screen

A coderdata dataset of cutaneous neurofibroma (cNF) patient-derived organoids
with drug-response measurements and matched multi-omic profiling.

## Overview

Cutaneous neurofibromas are benign nerve-sheath tumors that develop in patients
with neurofibromatosis type 1 (NF1). This dataset contains drug-screen
viability measurements for 238 small-molecule compounds across 23 cNF tumor
specimens from 10 patients, paired with bulk RNA-seq and global LC-MS/MS
proteomics for many of the same specimens.

| Property | Value |
|---|---|
| Cancer type | Cutaneous Neurofibroma |
| Model type | patient derived organoid |
| Patients | 10 |
| Specimens | 23 |
| Drugs | 238 |
| Modalities included | Transcriptomics (TPM), Proteomics (log ratio) |
| Modalities excluded | Phospho-proteomics (no schema slot) |
| Cohorts | 2 (encoded in `other_id_source`) |

## Multi-tumor patient structure

Each patient (`NF0017`, `NF0018`, ...) contributed up to three independent
tumor specimens, named `NFxxxx_T1`, `NFxxxx_T2`, `NFxxxx_T3`. This dataset
treats each specimen as an independent sample (its own `improve_sample_id`).
The patient ID is preserved in the `other_names` column so downstream models
can group specimens by patient when needed (e.g. patient-level holdout).

```
improve_sample_id  other_id    common_name  other_names  other_id_source
2401               NF0019_T1   NF0019_T1    NF0019       cNF_Cohort_1
2402               NF0019_T2   NF0019_T2    NF0019       cNF_Cohort_1
2403               NF0019_T3   NF0019_T3    NF0019       cNF_Cohort_1
```

## Cohort structure

Specimens were collected and processed in two cohorts. Cohort assignment is
encoded in the `other_id_source` column of the Sample table:

| Cohort         | Patients                                       |
|----------------|------------------------------------------------|
| cNF_Cohort_1   | NF0017, NF0018, NF0019, NF0020, NF0021, NF0022, NF0023 |
| cNF_Cohort_2   | NF0025, NF0027, NF0035                         |

Note: per-modality cohort membership can differ from patient-level cohort
(e.g. one specimen had RNA processed in Cohort 1 but proteomics in Cohort 2).
For coderdata, sample-level cohort follows the drug-screen batch.

## Drug-response metrics

The cNF screen mixes two designs:

* **Multi-dose drugs**: a subset of drugs were screened at multiple
  concentrations. These are run through `coderbuild/utils/fit_curve.py` and
  contribute the standard schema metrics (`fit_auc`, `fit_ic50`, `fit_einf`,
  `fit_hs`, `fit_r2`).
* **Single-dose drugs**: remaining drugs were screened at a single 1 μM
  concentration. These are recorded with `dose_response_metric = uM_viability`
  and `dose_response_value = Viability_percentage / 100` (range 0–1, where
  1.0 = full viability and 0.0 = full kill).

`uM_viability` is a new value in the `ResponseMetric` enum and requires a
schema PR alongside this dataset. See `schema_patch.md`.

## Data sources

| Modality | Synapse ID | Notes |
|---|---|---|
| Drug-screen file index | `syn51301431` (table) | dataType='drug screen' |
| Drug-screen files (parent) | `syn69947322` | One file per cohort |
| RNA-seq (gene-level TPM) | set via `CNF_RNA_TPM_SYN_ID` env var | |
| Global proteomics | `syn70078416` | log-ratio (correctedAbundance) |
| Phospho-proteomics | `syn70078415` | **excluded** from coderdata |

## Build pipeline

The build follows coderdata conventions: four shell scripts wrapping Python
implementations, packaged in a Docker image.

```
coderbuild/cnf/
├── build_samples.py      # mint improve_sample_ids, write cnf_samples.csv
├── build_samples.sh
├── build_omics.py        # write cnf_transcriptomics.csv + cnf_proteomics.csv
├── build_omics.sh
├── build_drugs.py        # call pubchem_retrieval + descriptor utility
├── build_drugs.sh
├── build_exp.py          # split single/multi-dose, run fit_curve, write cnf_experiments.tsv
├── build_exp.sh
├── requirements.txt
└── README.md             (this file)

coderbuild/docker/
└── Dockerfile.cnf
```

### Required environment variables

| Variable | Purpose |
|---|---|
| `SYNAPSE_AUTH_TOKEN` | Synapse auth (required for all scripts) |
| `CNF_RNA_TPM_SYN_ID` | Synapse ID of gene-level TPM file (required for omics build) |
| `CNF_GLOBAL_PROT_SYN_ID` | Optional override for global proteomics (default `syn70078416`) |

### Local build

```bash
# From repository root
docker build -f coderbuild/docker/Dockerfile.cnf -t coderdata-cnf .

docker run --rm \
    -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN \
    -e CNF_RNA_TPM_SYN_ID=syn... \
    -v $(pwd)/build_outputs:/coderbuild/cnf/out \
    coderdata-cnf bash build_samples.sh /coderbuild/cnf/out/prev_samples.csv

# repeat for build_omics.sh, build_drugs.sh, build_exp.sh in order
```

In the full coderdata build, this dataset is invoked by `build_dataset.py` /
`build_all.py` after at least one prior dataset has populated the running
`samples.csv` / `drugs.tsv` / `genes.csv` files.

## Schema dependency

This dataset uses `dose_response_metric = uM_viability`, which is **not**
currently in the `ResponseMetric` enum. Either submit the schema patch with
this dataset's PR or land the schema patch first:

```yaml
# schema/coderdata.yaml, ResponseMetric enum
uM_viability:
  description: >-
    Single-dose viability fraction at 1 μM. Range 0–1, where 1.0 = no
    effect (full viability) and 0.0 = full kill. Used for organoid drug
    screens where dose-response curves were not collected for that drug.
```

There is also a pre-existing inconsistency where the `dose_response_metric`
slot has `range: CurveMetric` while the enum is named `ResponseMetric`. The
schema patch should resolve this — flagged in `schema_patch.md`.

## What's intentionally not included

* **Phospho-proteomics** — useful in our internal modeling but no current
  schema slot. Re-evaluate when coderdata adds a `Phosphoproteomics` class.
* **Mutations** — not generated for this cohort.
* **Copy number** — not generated for this cohort.
* **Drug-class annotations** — handled by the standard
  `pubchem_retrieval.py` / descriptor pipeline rather than a custom file.

## Citation

To be added when the cNF organoid drug-screening manuscript is published.
The drug-screen analysis pipeline lives at
[PNNL-CompBio/cNFDrugScreening](https://github.com/PNNL-CompBio/cNFDrugScreening).