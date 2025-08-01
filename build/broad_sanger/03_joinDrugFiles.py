#!/usr/bin/env python3
import polars as pl
import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Merge drug TSVs, dedupe, and write single output.")
    parser.add_argument("inputs", nargs="+", help="Input TSV file(s) to merge (must share same header).")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path.")
    args = parser.parse_args()

    if len(args.inputs) < 1:
        print("Need at least one input file.", file=sys.stderr)
        sys.exit(1)

    dfs = []
    for p in args.inputs:
        if not Path(p).exists():
            print(f"Input file '{p}' does not exist; skipping.", file=sys.stderr)
            continue
        try:
            df = pl.read_csv(p, sep="\t", ignore_errors=True)
            dfs.append(df)
        except Exception as e:
            print(f"Failed to read '{p}': {e}", file=sys.stderr)

    if not dfs:
        print("No input dataframes could be read. Exiting.", file=sys.stderr)
        sys.exit(1)

    combined = pl.concat(dfs, how="vertical").unique()

    # Write out as TSV
    combined.write_csv(args.output, separator="\t")
    print(f"Wrote {combined.height} unique rows to {args.output}")

if __name__ == "__main__":
    main()
