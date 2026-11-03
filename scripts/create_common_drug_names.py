#!/usr/bin/env python
"""Build a common drug-name table from CoderData drug files.

The script uses local files supplied with ``--data-dir`` or downloads CoderData
files through the package, then collects the ``improve_drug_id``/``chem_name``
pairs from each drug file. It keeps every ``improve_drug_id`` unchanged,
checks that each ID has a consistent InChIKey or canonical SMILES, and selects
one source ``chem_name`` for display. Selection favors concise names over
registry or catalog identifiers (such as ChEMBL IDs, NSC numbers, and
InChIKeys), deprioritizes names marked obsolete, and then uses standardized
name annotations, source-file frequency, and alphabetical order to break ties.
The output name must be an exact ``chem_name`` value associated with that same
ID in the input files.
"""

import argparse
from collections import Counter, defaultdict
from pathlib import Path
import re
import tempfile
from typing import Iterable, Optional

import pandas as pd
import coderdata


_DRUG_FILE_PATTERN = re.compile(r"_drugs\.(?:tsv|csv)(?:\.gz)?$", re.IGNORECASE)
# These patterns identify names that are useful as aliases, but poor display
# names (for example database accessions, catalog numbers, or InChIKeys).
_IDENTIFIER_NAME_PATTERNS = (
    re.compile(r"[a-z]{14}-[a-z]{10}-[a-z]", re.IGNORECASE),
    re.compile(r"dtxsid\d+", re.IGNORECASE),
    re.compile(r"refchem[:\s-]*\w+", re.IGNORECASE),
    re.compile(r"nci60[_-]\d+", re.IGNORECASE),
    re.compile(r"schembl\d+", re.IGNORECASE),
    re.compile(r"chembl\d+", re.IGNORECASE),
    re.compile(r"bdbm\d+", re.IGNORECASE),
    re.compile(r"akos\d+", re.IGNORECASE),
    re.compile(r"nsc[-\s]?\d+", re.IGNORECASE),
    re.compile(r"(?:pubchem\s*)?(?:cid\s*)?\d+", re.IGNORECASE),
    re.compile(r"unii[-\s]?[a-z0-9]+", re.IGNORECASE),
    re.compile(r"\d{2,7}-\d{2}-\d"),
    re.compile(r"(?:hms|mls|tox21|ncgc|brd-k|sr-)\w+", re.IGNORECASE),
    re.compile(r"[a-z]{1,5}\d{3,}[a-z0-9-]*", re.IGNORECASE),
)
# These annotations mark a standardized drug-name synonym in source files.
_STANDARD_NAME_ANNOTATION = re.compile(
    r"(?:\[(?:inn|usan|jan|who-dd|mart\.?)\]"
    r"|\((?:inn|usan|jan/usan)\))\s*$",
    re.IGNORECASE,
)


def _is_identifier_name(name: str) -> bool:
    return any(pattern.fullmatch(name.strip()) for pattern in _IDENTIFIER_NAME_PATTERNS)


def _normalize_chem_name(name: str) -> str:
    # Normalize case, strip trailing source annotations, and standardize code
    # separators for ranking. The chosen output is later taken from the raw
    # source value so normalization never rewrites a source chem_name.
    normalized = name.strip().casefold()
    while True:
        cleaned = re.sub(r"\s*(?:\[[^\]]*\]|\([^)]*\))\s*$", "", normalized)
        if cleaned == normalized:
            normalized = normalized.strip(" ,;-")
            return re.sub(
                r"\b([a-z]{1,5})[ -]+(\d{2,})\b",
                r"\1-\2",
                normalized,
            )
        normalized = cleaned


def _has_standard_name_annotation(name: str) -> bool:
    return _STANDARD_NAME_ANNOTATION.search(name.strip()) is not None


def _has_obsolete_name_annotation(name: str) -> bool:
    return re.search(r"\bobsolete\b", name, re.IGNORECASE) is not None


def _is_common_name_candidate(name: str) -> bool:
    # Short alphabetic names and drug codes are preferred to long chemical
    # descriptions, while identifier patterns are scored separately.
    return (
        len(name.strip()) <= 50
        and len(name.strip().split()) <= 5
        and re.fullmatch(r"[a-zA-Z][a-zA-Z0-9 -]*", name.strip()) is not None
    )


def _common_names_from_files(drug_files: Iterable[Path]) -> pd.DataFrame:
    """Choose one source chem_name per ID, rejecting structural ID conflicts."""
    mapping_file_counts = Counter()
    standard_name_file_counts = Counter()
    obsolete_name_file_counts = Counter()
    variant_file_counts = Counter()
    source_pairs = set()
    identity_values = {
        column: defaultdict(set) for column in ("InChIKey", "canSMILES")
    }

    for drug_file in drug_files:
        # Accept the repository's tab- and comma-delimited files, including
        # compressed files, and read structure columns only when present.
        separator = (
            "," if drug_file.name.lower().endswith((".csv", ".csv.gz")) else "\t"
        )
        available_columns = pd.read_csv(drug_file, sep=separator, nrows=0).columns
        identity_columns = [
            column for column in identity_values if column in available_columns
        ]
        drugs = pd.read_csv(
            drug_file,
            sep=separator,
            compression="infer",
            usecols=["improve_drug_id", "chem_name", *identity_columns],
            dtype=str,
        )
        drugs = drugs.dropna(subset=["improve_drug_id", "chem_name"])

        for column in identity_columns:
            values = drugs[["improve_drug_id", column]].dropna()
            for drug_id, value in values.itertuples(index=False, name=None):
                value = value.strip()
                if value:
                    identity_values[column][drug_id].add(
                        value.upper() if column == "InChIKey" else value
                    )

        # Count each normalized alias at most once per file, so a dataset with
        # repeated synonym rows does not outweigh other dataset files.
        pairs = set()
        standard_name_pairs = set()
        obsolete_name_pairs = set()
        variants = set()
        named_drugs = drugs.loc[
            (drugs["improve_drug_id"] != "")
            & (drugs["chem_name"].str.strip() != ""),
            ["improve_drug_id", "chem_name"],
        ]
        for drug_id, chem_name in named_drugs.itertuples(index=False, name=None):
            normalized_name = _normalize_chem_name(chem_name)
            if not normalized_name:
                continue
            pair = (drug_id, normalized_name)
            pairs.add(pair)
            # Keep raw pairs for the final exact source-mapping validation.
            source_pairs.add((drug_id, chem_name))
            variants.add((drug_id, normalized_name, chem_name))
            if _has_standard_name_annotation(chem_name):
                standard_name_pairs.add(pair)
            if _has_obsolete_name_annotation(chem_name):
                obsolete_name_pairs.add(pair)
        mapping_file_counts.update(pairs)
        standard_name_file_counts.update(standard_name_pairs)
        obsolete_name_file_counts.update(obsolete_name_pairs)
        variant_file_counts.update(variants)

    # IDs must resolve to one chemical structure across all downloaded files.
    conflicts = [
        (drug_id, column, sorted(values))
        for column, values_by_id in identity_values.items()
        for drug_id, values in values_by_id.items()
        if len(values) > 1
    ]
    if conflicts:
        details = "; ".join(
            f"{drug_id} has conflicting {column} values: {values}"
            for drug_id, column, values in conflicts[:10]
        )
        raise ValueError(
            "CoderData improve_drug_id values do not consistently identify "
            f"one chemical structure. {details}"
        )

    names_by_id = defaultdict(list)
    for (drug_id, chem_name), file_count in mapping_file_counts.items():
        names_by_id[drug_id].append(
            (
                chem_name,
                file_count,
                standard_name_file_counts[(drug_id, chem_name)],
                obsolete_name_file_counts[(drug_id, chem_name)],
            )
        )

    variants_by_name = defaultdict(list)
    for (drug_id, normalized_name, source_name), file_count in variant_file_counts.items():
        variants_by_name[(drug_id, normalized_name)].append((source_name, file_count))

    rows = []
    for drug_id, candidates in names_by_id.items():
        # Rank aliases by: non-identifier, concise/common-looking, not marked
        # obsolete, standardized-name support, file frequency, then alphabetic
        # order. Choose the raw source spelling within the winning alias.
        selected_name = min(
            candidates,
            key=lambda candidate: (
                _is_identifier_name(candidate[0]),
                not _is_common_name_candidate(candidate[0]),
                bool(candidate[3]),
                -candidate[2],
                -candidate[1],
                candidate[0].casefold(),
                candidate[0],
            ),
        )[0]
        source_name = min(
            variants_by_name[(drug_id, selected_name)],
            key=lambda variant: (
                _is_identifier_name(variant[0]),
                _has_obsolete_name_annotation(variant[0]),
                _has_standard_name_annotation(variant[0]),
                not _is_common_name_candidate(variant[0]),
                -variant[1],
                variant[0].casefold(),
                variant[0],
            ),
        )[0]
        rows.append((drug_id, source_name))

    # This guards against changing IDs or synthesizing names during selection.
    if not set(rows).issubset(source_pairs):
        raise ValueError(
            "Generated mappings include improve_drug_id/chem_name pairs not "
            "present in the downloaded source drug files."
        )

    return pd.DataFrame(
        sorted(rows),
        columns=["improve_drug_id", "common_name"],
    )


def create_common_drug_names(
    output_path: Path, data_dir: Optional[Path] = None
) -> int:
    """Build the mapping from local files, or download files if none are given."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def find_drug_files(directory: Path):
        if not directory.is_dir():
            raise NotADirectoryError(f"CoderData directory does not exist: {directory}")
        files = sorted(
            path
            for path in directory.iterdir()
            if path.is_file() and _DRUG_FILE_PATTERN.search(path.name)
        )
        if not files:
            raise FileNotFoundError(f"No CoderData drug files found in {directory}.")
        return files

    if data_dir is not None:
        # Reuse an existing CoderData download and avoid network or extraction.
        mapping = _common_names_from_files(find_drug_files(Path(data_dir)))
    else:
        with tempfile.TemporaryDirectory(prefix="coderdata-drugs-") as download_dir:
            # Use the package's public downloader. It downloads all dataset files;
            # only drug files are selected and read, and temporary files are
            # removed automatically after the mapping has been written.
            coderdata.download(name="all", local_path=download_dir)
            mapping = _common_names_from_files(find_drug_files(Path(download_dir)))
    mapping.to_csv(output_path, index=False, encoding="utf-8-sig")

    return len(mapping)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Download CoderData files and create a common drug-name CSV."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("common_drug_names.csv"),
        help="destination CSV path (default: ./common_drug_names.csv)",
    )
    parser.add_argument(
        "-d",
        "--data-dir",
        type=Path,
        help=(
            "directory containing already-downloaded CoderData drug files "
            "(*_drugs.tsv[.gz] or *_drugs.csv[.gz]); skips downloading."
        ),
    )
    args = parser.parse_args(argv)
    row_count = create_common_drug_names(args.output, data_dir=args.data_dir)
    print(f"Wrote {row_count} drug-name mappings to {args.output}")


if __name__ == "__main__":
    main()
