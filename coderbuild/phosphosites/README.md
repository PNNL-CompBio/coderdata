# Phosphosites Reference Builder

Produces `phosphosites.csv`, the global human phosphosite reference used by
dataset omics pipelines (currently: cNF) to map raw phospho measurements to
stable integer identifiers.

## Sources

Three sources are merged in priority order. A site from an earlier source is
never overwritten by a later one (deduplication on `other_id`).

| Priority | Source | Sites (approx.) |
|---|---|---|
| 1 | Ochoa et al. (2020) Nat Biotechnol 38:365-373, Supp. Table S3 | ~116k |
| 2 | UniProt PTM REST API | ~12k additional |
| 3 | Synapse supplement `syn70078415` (cNF raw phospho data) | ~188 additional |

## Site notation

```
<GENE>-<Residue><Position><mod_code>
e.g.  AAAS-S495s  (s=phosphoserine, t=phosphothreonine, y=phosphotyrosine)
```

`other_id` in `phosphosites.csv` stores these strings and is the join key for
raw phospho data.

## Output columns

`phosphosite_id`, `entrez_id`, `gene_symbol`, `residue`, `position`,
`modification`, `other_id`

## Files

```
coderbuild/phosphosites/
├── 00-buildPhosphositeFile.py   # Builder script
└── build_phosphosites.sh        # Docker entry point
coderbuild/docker/
└── Dockerfile.phosphosites
```

## Builder — `00-buildPhosphositeFile.py`

**Ochoa (primary):** Downloads `41587_2019_344_MOESM4_ESM.xlsx` from the
Springer static CDN (no login required; requires `openpyxl`). Reads the
`annotated_phosphoproteome` sheet (`uniprot`, `position`, `residue` columns).
UniProt accessions are resolved to gene symbols via a separate UniProt REST
call (TSV, all reviewed human proteins, paginated via `Link: rel="next"`
header).

**UniProt PTM (secondary):** Downloads `Modified residue` features with
phosphoserine/phosphothreonine/phosphotyrosine descriptions from the UniProt
REST API (JSON, same Link-header pagination). Returns gene symbols directly.

**Synapse supplement (`--synapse_supplements`):** Downloads a raw phospho
file from Synapse, auto-detects the site column (`site`, `Site`,
`phosphosite`, or `feature_id`), and appends any `GENE-ResPositionmod` sites
not already in the reference. Requires `SYNAPSE_AUTH_TOKEN`. Failures are
logged as warnings and do not abort the build.

### Key arguments

| Argument | Default | Description |
|---|---|---|
| `genes` (positional) | — | Path to `genes.csv` |
| `--out` | `/tmp/phosphosites.csv` | Output path |
| `--prev` | None | Previous CSV to preserve `phosphosite_id` values across builds |
| `--synapse_supplements` | None | Synapse IDs of raw phospho files (space-separated) |
| `--synapse_site_col` | `site` | Site column name in Synapse supplement files |
| `--supplement` | None | Local CSV/TSV file(s) with site strings |
| `--site_col` | `site` | Site column name in local supplement files |

## Build sequencing

Phosphosites must run **after** `genes.csv` is ready. `build_dataset.py` and
`build_all.py` wait on the genes future before submitting the phosphosites
container.

## ID stability

Pass `--prev <prev_phosphosites.csv>` to preserve IDs across rebuilds. New
sites receive IDs starting from `max(prev.phosphosite_id) + 1` — same pattern
as `improve_sample_id` and `improve_drug_id`.

## Local run (Only for debugging purposes, this is already build into build_all.py script)

```bash
pip install pandas synapseclient openpyxl

export SYNAPSE_AUTH_TOKEN=<token>
python coderbuild/phosphosites/00-buildPhosphositeFile.py \
    local/genes.csv \
    --out local/phosphosites.csv \
    --synapse_supplements syn70078415 \
    --synapse_site_col site
```
