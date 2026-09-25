#!/usr/bin/env python3
"""Convert a LAVA locus file to the headered region format SUPERGNOVA reads."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--locfile", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    loci = pd.read_csv(args.locfile, sep=r"\s+")
    required = ["CHR", "START", "STOP"]
    missing = set(required) - set(loci.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    out = loci[required].rename(columns={"STOP": "END"}).sort_values(["CHR", "START"]).reset_index(drop=True)
    if len(out) != 3064:
        raise ValueError(f"Expected 3064 LAVA regions, found {len(out)}")
    if (out["START"] > out["END"]).any():
        raise ValueError("At least one region has START > END")
    for chromosome, group in out.groupby("CHR"):
        if (group["START"].iloc[1:].to_numpy() <= group["END"].iloc[:-1].to_numpy()).any():
            raise ValueError(f"Overlapping intervals found on chromosome {chromosome}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
