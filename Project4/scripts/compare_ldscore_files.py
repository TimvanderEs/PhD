#!/usr/bin/env python3
"""Compare two LD-score files by SNP and report numerical equivalence."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def read(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().split()
        snp_index = header.index("SNP")
        l2_index = header.index("L2")
        for line in handle:
            fields = line.split()
            values[fields[snp_index]] = float(fields[l2_index])
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    original = read(args.original)
    candidate = read(args.candidate)
    if set(original) != set(candidate):
        raise ValueError("LD-score SNP sets differ")
    differences = [abs(original[marker] - candidate[marker]) for marker in original]
    row = {
        "original": str(args.original),
        "candidate": str(args.candidate),
        "rows": len(differences),
        "max_absolute_difference": max(differences),
        "mean_absolute_difference": sum(differences) / len(differences),
        "exact_equal_at_stored_precision": int(all(value == 0 for value in differences)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
