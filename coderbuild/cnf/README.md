# Cutaneous Neurofibroma (cnf) Organoid Drug Screen

A coderdata dataset of cutaneous neurofibroma (cNF) patient-derived organoids
with drug-response measurements and matched multi-omic profiling.

## Overview

Cutaneous neurofibromas are benign nerve-sheath tumors that develop in patients
with neurofibromatosis type 1 (NF1). This dataset contains drug-screen
viability measurements for 238 small-molecule compounds across 23 cNF tumor
specimens from 10 patients, paired with bulk RNA-seq, global LC-MS/MS
proteomics, and batch-corrected phosphoproteomics for many of the same
specimens.

| Property | Value |
|---|---|
| Cancer type | Cutaneous Neurofibroma |
| Model type | patient derived organoid |
| Patients | 10 |
| Specimens | 23 (untreated organoids; additional tissue and treated organoid samples also present) |
| Drugs | 238 |
| Modalities included | Transcriptomics (TPM), Proteomics (log ratio), Phosphoproteomics (log ratio) |

## Multi-tumor patient structure

Each patient (`NF0017`, `NF0018`, ...) contributed up to three independent
tumor specimens, named `NFxxxx_T1`, `NFxxxx_T2`, `NFxxxx_T3`. This dataset
treats each specimen as an independent sample (its own `improve_sample_id`).
The patient ID is preserved in the `other_names` column so downstream models
can group specimens by patient when needed (e.g. patient-level holdout).

```
improve_sample_id  other_id    common_name  other_names  other_id_source
2401               NF0019_T1   NF0019_T1    NF0019       NF1_cNF_Project
2402               NF0019_T2   NF0019_T2    NF0019       NF1_cNF_Project
2403               NF0019_T3   NF0019_T3    NF0019       NF1_cNF_Project
```

The sample table also includes matched normal skin (`NFxxxx_skin`), primary
tumor tissue (`NFxxxx_Tx_tissue`), and drug-treated organoids
(`NFxxxx_Tx_<treatment_label>`). Each gets its own `improve_sample_id` and
`model_type`. Organoids (untreated) receive the lowest IDs, sorted by
specimen ID, to make drug-modelable samples easy to filter.


## Drug-response metrics

The cNF screen mixes two designs:

- **Multi-dose drugs**: a subset of drugs were screened at multiple
  concentrations. These are run through `coderbuild/utils/fit_curve.py` and
  contribute the standard schema metrics (`fit_auc`, `fit_ic50`, `fit_einf`,
  `fit_hs`, `fit_r2`).

- **Single-dose drugs**: remaining drugs were screened at a single 1 μM
  concentration. These are recorded with `dose_response_metric = uM_viability`
  and `dose_response_value = Viability_percentage / 100` (range 0–1, where
  1.0 = full viability and 0.0 = full kill).

`uM_viability` is a value in the `ResponseMetric` enum. See the Schema
dependency section below.

## SPECIMEN_DUAL_MAPPINGS edge case

One viability file covers the treated organoid `NF0021_T1_Onalespid_1uM` but
its drug-response measurements should also be attributed to the untreated
`NF0021_T1`. The `SPECIMEN_DUAL_MAPPINGS` dict in `04-experiments-cnf.py`
handles this by duplicating the rows for every additional specimen listed:

```python
SPECIMEN_DUAL_MAPPINGS = {
    "NF0021_T1_Onalespid_1uM": ["NF0021_T1"],
}
```

Add new entries here if additional treated organoids need the same dual
attribution. Both specimens must already exist in `cnf_samples.csv` or the
duplicate rows will be dropped during the sample-map join.

## Protocol optimization samples

Six RNA columns in the cohort TPM matrices correspond to protocol optimization
runs rather than real specimens. They are excluded in `02-omics-cnf.py` before
the wide-to-long melt. Column stems matching any of these substrings are
dropped:

```python
DROP_SAMPLE_SUBSTRINGS = [
    "cNF_organoid_DIA_G_02_11Feb25",
    "cNF_organoid_DIA_G_05_11Feb25",
    "cNF_organoid_DIA_G_06_11Feb25",
    "cNF_organoid_DIA_P_02_29Jan25",
    "cNF_organoid_DIA_P_05_11Feb25",
    "cNF_organoid_DIA_P_06_11Feb25",
]
```

These are matched on the column stem (basename without extension), so both
plain paths and filename columns are handled correctly. If new protocol
optimization runs appear in future cohort data, add their column stems to
this list.

## Phosphosites dependency and known unmatched sites

`02-omics-cnf.py` requires `phosphosites.csv` as a positional argument.
The phosphosites reference must be built before the cNF omics step runs.
`build_dataset.py` and `build_all.py` enforce this sequencing.

The phosphoproteomics output file (`cnf_phosphoproteomics.csv`) is produced
by an inner join between the raw Synapse phospho data and `phosphosites.csv`.
Sites that don't appear in the reference are silently dropped. As of the
current release, **12 sites** in the raw cNF phospho data (syn70078415) are
permanently unmatched because their gene symbols (`BAP18`, `C11orf96`,
`C14orf93`, `C1orf21`, `C2orf49`, etc). The remaining 2,974 of 2,986 unique sites
(99.6%) are present in the output.

## Data sources

| Modality | Synapse ID | Notes |
|---|---|---|
| Drug-screen file index | `syn51301431` (table) | `dataType='drug screen'` |
| Drug-screen file parents | `syn51301414`, `syn51301420`, `syn51301426` | One folder per cohort |
| RNA-seq TPM (Cohort 1) | `syn66352931` | `salmon.merged.gene_tpm.tsv` |
| RNA-seq TPM (Cohort 2) | `syn70765053` | `salmon.merged.gene_tpm.tsv` |
| Global proteomics | `syn74815895` | Batch-corrected; `correctedAbundance` column |
| Phospho-proteomics | `syn70078415` | Batch-corrected; `correctedAbundance` column; mapped via `phosphosites.csv` |
| Sample discovery (RNA) | `syn71333780` | All-conditions RNA file for sample enumeration only |
| Normal Skin folder | `syn74284682` | One child per patient for skin sample discovery |

## Build pipeline

The build follows coderdata conventions: four shell scripts wrapping Python
implementations, plus shared utilities, packaged in a Docker image.

```
coderbuild/cnf/
├── 01-samples-cnf.py     # mint improve_sample_ids, write cnf_samples.csv
├── build_samples.sh
├── 02-omics-cnf.py       # write cnf_transcriptomics.csv, cnf_proteomics.csv,
├── build_omics.sh        #   cnf_phosphoproteomics.csv
├── 03-drugs-cnf.py       # call pubchem_retrieval + descriptor utility
├── build_drugs.sh
├── 04-experiments-cnf.py # split single/multi-dose, run fit_curve, write cnf_experiments.tsv
├── build_exp.sh
├── cnf_utils.py          # shared specimen canonicalization utilities
├── requirements.txt
└── README.md             (this file)

coderbuild/docker/
└── Dockerfile.cnf
```

### Local build

```bash
# From repository root
python coderbuild/build_dataset.py  --build --validate --dataset cnf
```

In the full coderdata build, this dataset is orchestrated by `build_dataset.py`
/ `build_all.py`, which ensures `genes.csv`, `samples.csv`, `drugs.tsv`, and
`phosphosites.csv` are all present before starting the cNF containers.


## What's intentionally not included

- **Mutations** - not yet generated for this dataset.
- **Copy number** - not yet generated for this dataset.

