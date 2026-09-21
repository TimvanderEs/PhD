#!/usr/bin/env python3
"""Summarize shared-variant Z-score concordance across HELIOS GWAS freezes."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


def read_sumstats(path: Path) -> dict[str, tuple[float, float]]:
    values: dict[str, tuple[float, float]] = {}
    with gzip.open(path, "rt") as stream:
        header = stream.readline().split()
        snp_i, n_i, z_i = (header.index(name) for name in ("SNP", "N", "Z"))
        for line in stream:
            fields = line.split()
            try:
                values[fields[snp_i]] = (float(fields[n_i]), float(fields[z_i]))
            except (IndexError, ValueError):
                continue
    return values


def moments(xs: list[float], ys: list[float]) -> dict[str, float]:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return {
        "pearson_r": sxy / math.sqrt(sxx * syy),
        "slope_y_on_x": sxy / sxx,
        "rmse": math.sqrt(sum((x - y) ** 2 for x, y in zip(xs, ys)) / n),
        "sign_concordance": sum((x >= 0) == (y >= 0) for x, y in zip(xs, ys)) / n,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ten-k", type=Path, required=True)
    parser.add_argument("--twelve-k", type=Path, required=True)
    parser.add_argument("--twenty-two-k", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = {
        "trans10k": read_sumstats(args.ten_k),
        "extension12k": read_sumstats(args.twelve_k),
        "combined22k": read_sumstats(args.twenty_two_k),
    }
    rows: list[dict[str, object]] = []
    pairs = (("trans10k", "extension12k"), ("trans10k", "combined22k"), ("extension12k", "combined22k"))
    for first, second in pairs:
        shared = sorted(set(data[first]).intersection(data[second]))
        xs = [data[first][snp][1] for snp in shared]
        ys = [data[second][snp][1] for snp in shared]
        rows.append({"qc": args.label, "comparison": f"{first}_vs_{second}", "n_shared": len(shared), **moments(xs, ys)})

    shared_three = sorted(set.intersection(*(set(values) for values in data.values())))
    observed = [data["combined22k"][snp][1] for snp in shared_three]
    for orientation in (1, -1):
        predicted = []
        for snp in shared_three:
            n1, z1 = data["trans10k"][snp]
            n2, z2 = data["extension12k"][snp]
            denominator = math.sqrt(n1 + n2)
            predicted.append((math.sqrt(n1) * z1 + orientation * math.sqrt(n2) * z2) / denominator)
        rows.append(
            {
                "qc": args.label,
                "comparison": f"combined22k_vs_Nweighted_components_orientation_{orientation:+d}",
                "n_shared": len(shared_three),
                **moments(predicted, observed),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
