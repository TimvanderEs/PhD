#!/usr/bin/env python3
"""Quantify reference-size LD-score stability on three sentinel chromosomes."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics
from pathlib import Path


def read_scores(path: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().split()
        snp_index = header.index("SNP")
        l2_index = header.index("L2")
        for line in handle:
            fields = line.split()
            scores[fields[snp_index]] = float(fields[l2_index])
    return scores


def pearson(a: list[float], b: list[float]) -> float:
    mean_a, mean_b = statistics.fmean(a), statistics.fmean(b)
    numerator = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    denom = math.sqrt(sum((x - mean_a) ** 2 for x in a) * sum((y - mean_b) ** 2 for y in b))
    return numerator / denom


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    full_root = args.root / "panels" / "chinese_standard"
    raw_root = args.root / "panels" / "reference_size_raw"
    sentinel_root = args.root / "panels" / "reference_size_sentinel"
    rows: list[dict[str, object]] = []
    for chromosome in (1, 6, 22):
        full = read_scores(full_root / f"{chromosome}.l2.ldscore.gz")
        for size in (500, 2500, 10000):
            for replicate in range(1, 6):
                if replicate == 1:
                    path = raw_root / f"chinese_n{size}" / f"{chromosome}.l2.ldscore.gz"
                else:
                    path = sentinel_root / f"chinese_n{size}_r{replicate}" / f"{chromosome}.l2.ldscore.gz"
                subset = read_scores(path)
                common = sorted(set(full).intersection(subset))
                observed = [subset[snp] for snp in common]
                expected = [full[snp] for snp in common]
                differences = [x - y for x, y in zip(observed, expected)]
                rows.append(
                    {
                        "chromosome": chromosome,
                        "sample_size": size,
                        "replicate": replicate,
                        "common_snps": len(common),
                        "pearson_with_full": pearson(observed, expected),
                        "mean_ldscore_difference": statistics.fmean(differences),
                        "rmse_vs_full": math.sqrt(statistics.fmean(value * value for value in differences)),
                    }
                )

    summary: list[dict[str, object]] = []
    for size in (500, 2500, 10000):
        subset = [row for row in rows if row["sample_size"] == size]
        summary.append(
            {
                "sample_size": size,
                "chromosome_replicate_comparisons": len(subset),
                "mean_pearson_with_full": statistics.fmean(float(row["pearson_with_full"]) for row in subset),
                "sd_pearson_with_full": statistics.stdev(float(row["pearson_with_full"]) for row in subset),
                "mean_ldscore_difference": statistics.fmean(float(row["mean_ldscore_difference"]) for row in subset),
                "mean_rmse_vs_full": statistics.fmean(float(row["rmse_vs_full"]) for row in subset),
            }
        )

    for filename, output in (("reference_size_ldscore_sentinel.tsv", rows), ("reference_size_ldscore_summary.tsv", summary)):
        path = args.root / "results" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(output)


if __name__ == "__main__":
    main()
