#!/usr/bin/env python3
"""Center and unit-normalize PLINK eigenvectors for cov-LDSC."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eigenvec", required=True, type=Path)
    parser.add_argument("--components", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    identifiers: list[tuple[str, str]] = []
    rows: list[list[float]] = []
    with args.eigenvec.open() as handle:
        header = handle.readline().split()
        if len(header) < args.components + 2:
            raise ValueError("Eigenvector file has fewer components than requested")
        for line_number, line in enumerate(handle, start=2):
            fields = line.split()
            if len(fields) < args.components + 2:
                raise ValueError(f"Malformed eigenvector line {line_number}")
            identifiers.append((fields[0], fields[1]))
            rows.append([float(value) for value in fields[2 : args.components + 2]])
    if not rows or len({iid for _, iid in identifiers}) != len(identifiers):
        raise ValueError("Eigenvector file is empty or has duplicate IIDs")

    n = len(rows)
    means = [sum(row[j] for row in rows) / n for j in range(args.components)]
    norms = []
    for j in range(args.components):
        norm = math.sqrt(sum((row[j] - means[j]) ** 2 for row in rows))
        if norm == 0:
            raise ValueError(f"PC{j + 1} has zero norm")
        norms.append(norm)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for (fid, iid), row in zip(identifiers, rows):
            values = [(row[j] - means[j]) / norms[j] for j in range(args.components)]
            handle.write(f"{fid}\t{iid}\t" + "\t".join(f"{value:.12g}" for value in values) + "\n")

    checks = []
    for j in range(args.components):
        values = [(row[j] - means[j]) / norms[j] for row in rows]
        checks.append(
            {
                "component": j + 1,
                "mean": sum(values) / n,
                "l2_norm": math.sqrt(sum(value * value for value in values)),
            }
        )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(
            {
                "source": str(args.eigenvec),
                "output": str(args.output),
                "samples": n,
                "components": args.components,
                "checks": checks,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
