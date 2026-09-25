#!/usr/bin/env python3
"""Diagnose trait/reference SNP and allele overlap for one chromosome."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def classify(ref_a1: str, ref_a2: str, a1: str, a2: str) -> str:
    if (a1, a2) == (ref_a1, ref_a2) or (a1.translate(COMPLEMENT), a2.translate(COMPLEMENT)) == (ref_a1, ref_a2):
        return "same"
    if (a2, a1) == (ref_a1, ref_a2) or (a2.translate(COMPLEMENT), a1.translate(COMPLEMENT)) == (ref_a1, ref_a2):
        return "reversed"
    return "mismatch"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", required=True, type=Path)
    parser.add_argument("--sumstats", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    bim = pd.read_csv(
        args.bim,
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        usecols=[1, 4, 5],
        dtype="string",
    )
    reference = dict(zip(bim["SNP"], zip(bim["A1"], bim["A2"])))
    counts = {"reference_variants": len(reference), "shared_ids": 0, "same": 0, "reversed": 0, "mismatch": 0}
    for chunk in pd.read_csv(
        args.sumstats,
        sep=r"\s+",
        usecols=["SNP", "A1", "A2"],
        dtype="string",
        chunksize=500_000,
    ):
        chunk = chunk.loc[chunk["SNP"].isin(reference)]
        counts["shared_ids"] += len(chunk)
        for snp, a1, a2 in chunk.itertuples(index=False, name=None):
            counts[classify(*reference[snp], a1, a2)] += 1
    args.output.write_text(json.dumps(counts, indent=2) + "\n")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
