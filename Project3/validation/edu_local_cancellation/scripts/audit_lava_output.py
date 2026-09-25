#!/usr/bin/env python3
"""Check whether a LAVA bivariate file contains every eligible trait pair."""

from __future__ import annotations

import argparse
import csv
import itertools


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--univ", required=True)
    parser.add_argument("--bivar", required=True)
    parser.add_argument("--threshold", type=float, required=True)
    parser.add_argument("--phenotype-a", required=True)
    parser.add_argument("--phenotype-b", required=True)
    args = parser.parse_args()

    eligible: dict[str, list[str]] = {}
    with open(args.univ, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row["p"]) < args.threshold:
                eligible.setdefault(row["locus"], []).append(row["phen"])

    expected: set[tuple[str, str, str]] = set()
    for locus, phenotypes in eligible.items():
        for first, second in itertools.combinations(phenotypes, 2):
            expected.add((locus, first, second))

    actual: set[tuple[str, str, str]] = set()
    with open(args.bivar, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            actual.add((row["locus"], row["phen1"], row["phen2"]))

    target = {
        item
        for item in actual
        if {item[1], item[2]} == {args.phenotype_a, args.phenotype_b}
    }

    print("metric\tvalue")
    print(f"expected_bivariate_rows\t{len(expected)}")
    print(f"observed_bivariate_rows\t{len(actual)}")
    print(f"missing_expected_rows\t{len(expected - actual)}")
    print(f"unexpected_rows\t{len(actual - expected)}")
    print(f"target_pair_rows\t{len(target)}")


if __name__ == "__main__":
    main()
