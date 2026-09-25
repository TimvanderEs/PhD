#!/usr/bin/env python3
"""Convert published ASN LDetect BED intervals to SUPERGNOVA format."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    blocks = pd.read_csv(args.input, sep=r"\s+")
    blocks.columns = [column.strip().lower() for column in blocks.columns]
    blocks["CHR"] = blocks["chr"].str.replace("chr", "", regex=False).astype(int)
    blocks["START"] = blocks["start"].astype(int)
    # The published BED stop is half-open; SUPERGNOVA uses inclusive comparisons.
    blocks["END"] = blocks["stop"].astype(int) - 1
    output = blocks[["CHR", "START", "END"]].sort_values(["CHR", "START"]).reset_index(drop=True)
    if len(output) != 1445:
        raise ValueError(f"Expected 1445 Asian LDetect blocks, found {len(output)}")
    if output["CHR"].min() != 1 or output["CHR"].max() != 22:
        raise ValueError("Expected autosomal chromosomes 1-22")
    if (output["START"] > output["END"]).any():
        raise ValueError("At least one block has START > END")
    for chromosome, group in output.groupby("CHR"):
        if (group["START"].iloc[1:].to_numpy() <= group["END"].iloc[:-1].to_numpy()).any():
            raise ValueError(f"Overlapping intervals remain on chromosome {chromosome}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
