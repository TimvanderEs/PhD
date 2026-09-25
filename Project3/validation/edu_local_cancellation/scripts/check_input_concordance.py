#!/usr/bin/env python3
"""Compare LDSC and LAVA versions of the same summary statistics."""

from __future__ import annotations

import argparse
import gzip
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Effect:
    a1: str
    a2: str
    z: float


def read_effects(path: str) -> dict[str, Effect]:
    effects: dict[str, Effect] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        columns = {name: index for index, name in enumerate(header)}
        required = {"SNP", "A1", "A2", "Z"}
        missing = required.difference(columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            try:
                snp = fields[columns["SNP"]]
                a1 = fields[columns["A1"]].upper()
                a2 = fields[columns["A2"]].upper()
                z = float(fields[columns["Z"]])
            except (IndexError, ValueError):
                continue
            if snp and a1 and a2 and math.isfinite(z):
                effects[snp] = Effect(a1, a2, z)
    return effects


def compare(reference: dict[str, Effect], path: str) -> dict[str, float | int]:
    n_shared = 0
    n_same = 0
    n_swapped = 0
    n_mismatch = 0
    n_same_sign = 0
    sum_x = sum_y = sum_x2 = sum_y2 = sum_xy = 0.0
    sum_squared_difference = 0.0

    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        columns = {name: index for index, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            try:
                snp = fields[columns["SNP"]]
                current = Effect(
                    fields[columns["A1"]].upper(),
                    fields[columns["A2"]].upper(),
                    float(fields[columns["Z"]]),
                )
            except (IndexError, ValueError):
                continue
            original = reference.get(snp)
            if original is None or not math.isfinite(current.z):
                continue
            if (current.a1, current.a2) == (original.a1, original.a2):
                aligned_z = current.z
                n_same += 1
            elif (current.a1, current.a2) == (original.a2, original.a1):
                aligned_z = -current.z
                n_swapped += 1
            else:
                n_mismatch += 1
                continue
            n_shared += 1
            x = original.z
            y = aligned_z
            sum_x += x
            sum_y += y
            sum_x2 += x * x
            sum_y2 += y * y
            sum_xy += x * y
            sum_squared_difference += (y - x) ** 2
            n_same_sign += int((x >= 0) == (y >= 0))

    numerator = n_shared * sum_xy - sum_x * sum_y
    denominator = math.sqrt(
        (n_shared * sum_x2 - sum_x * sum_x)
        * (n_shared * sum_y2 - sum_y * sum_y)
    )
    correlation = numerator / denominator if denominator > 0 else math.nan
    slope = sum_xy / sum_x2 if sum_x2 > 0 else math.nan
    rmse = math.sqrt(sum_squared_difference / n_shared) if n_shared else math.nan
    return {
        "reference_variants": len(reference),
        "shared_aligned_variants": n_shared,
        "same_allele_order": n_same,
        "swapped_allele_order": n_swapped,
        "allele_mismatches": n_mismatch,
        "z_correlation": correlation,
        "aligned_z_slope": slope,
        "aligned_z_rmse": rmse,
        "sign_concordance": n_same_sign / n_shared if n_shared else math.nan,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--ldsc", required=True)
    parser.add_argument("--lava", required=True)
    args = parser.parse_args()

    result = compare(read_effects(args.ldsc), args.lava)
    print("comparison\tmetric\tvalue")
    for metric, value in result.items():
        print(f"{args.label}\t{metric}\t{value}")


if __name__ == "__main__":
    main()
