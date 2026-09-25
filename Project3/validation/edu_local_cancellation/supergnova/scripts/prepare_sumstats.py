#!/usr/bin/env python3
"""Prepare complete, strand-unambiguous summary statistics for SUPERGNOVA."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from pathlib import Path

import pandas as pd


RSID = re.compile(r"^rs[0-9]+$")
VALID_ALLELES = {"A", "C", "G", "T"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--sample-size", required=True, type=float)
    parser.add_argument("--chunksize", type=int, default=500_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    counters = {
        "input_rows": 0,
        "missing_required": 0,
        "nonfinite_z": 0,
        "non_rsid": 0,
        "non_snv": 0,
        "strand_ambiguous": 0,
        "exact_duplicate_rows_collapsed": 0,
        "conflicting_duplicate_rows_removed": 0,
        "output_rows": 0,
    }
    kept_chunks: list[pd.DataFrame] = []
    required = ["SNP", "A1", "A2", "Z"]

    for chunk in pd.read_csv(
        args.input,
        sep=r"\s+",
        usecols=required,
        dtype={"SNP": "string", "A1": "string", "A2": "string", "Z": "float64"},
        chunksize=args.chunksize,
    ):
        counters["input_rows"] += len(chunk)
        missing = chunk[required].isna().any(axis=1)
        counters["missing_required"] += int(missing.sum())
        chunk = chunk.loc[~missing].copy()

        finite_z = chunk["Z"].map(math.isfinite)
        counters["nonfinite_z"] += int((~finite_z).sum())
        chunk = chunk.loc[finite_z].copy()

        chunk["A1"] = chunk["A1"].str.upper()
        chunk["A2"] = chunk["A2"].str.upper()

        is_rsid = chunk["SNP"].str.match(RSID)
        counters["non_rsid"] += int((~is_rsid).sum())
        chunk = chunk.loc[is_rsid].copy()

        is_snv = chunk["A1"].isin(VALID_ALLELES) & chunk["A2"].isin(VALID_ALLELES)
        counters["non_snv"] += int((~is_snv).sum())
        chunk = chunk.loc[is_snv].copy()

        ambiguous = (
            ((chunk["A1"] == "A") & (chunk["A2"] == "T"))
            | ((chunk["A1"] == "T") & (chunk["A2"] == "A"))
            | ((chunk["A1"] == "C") & (chunk["A2"] == "G"))
            | ((chunk["A1"] == "G") & (chunk["A2"] == "C"))
        )
        counters["strand_ambiguous"] += int(ambiguous.sum())
        kept_chunks.append(chunk.loc[~ambiguous, required])

    data = pd.concat(kept_chunks, ignore_index=True)
    exact_duplicate = data.duplicated(subset=["SNP", "A1", "A2", "Z"], keep="first")
    counters["exact_duplicate_rows_collapsed"] = int(exact_duplicate.sum())
    data = data.loc[~exact_duplicate].copy()
    conflicting_duplicate = data["SNP"].duplicated(keep=False)
    counters["conflicting_duplicate_rows_removed"] = int(conflicting_duplicate.sum())
    data = data.loc[~conflicting_duplicate].copy()
    data.insert(3, "N", args.sample_size)
    data = data[["SNP", "A1", "A2", "N", "Z"]]
    counters["output_rows"] = len(data)

    with gzip.open(args.output, "wt") as handle:
        data.to_csv(handle, sep="\t", index=False, float_format="%.12g")

    args.report.write_text(json.dumps(counters, indent=2) + "\n")


if __name__ == "__main__":
    main()
