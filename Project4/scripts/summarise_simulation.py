#!/usr/bin/env python3
"""Summarize reference-panel accuracy in held-out genotype simulations."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


PANELS = ("1kgeas", "chinese_full_overlap", "chinese_disjoint_n10000")


def parse_h2(path: Path) -> tuple[float, float, float, float]:
    text = path.read_text(errors="replace")
    h2 = re.search(r"Total Observed scale h2:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", text)
    intercept = re.search(r"Intercept:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", text)
    if h2 is None or intercept is None:
        raise ValueError(f"Could not parse {path}")
    return (
        float(h2.group(1)),
        float(h2.group(2)),
        float(intercept.group(1)),
        float(intercept.group(2)),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--target-h2", type=float, required=True)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    for replicate in range(1, args.replicates + 1):
        rep = args.root / "simulation" / f"replicate_{replicate}"
        phenotype_qc = json.loads((rep / "phenotype_qc.json").read_text())
        for panel in PANELS:
            h2, se, intercept, intercept_se = parse_h2(rep / f"h2_{panel}.log")
            rows.append(
                {
                    "replicate": replicate,
                    "panel": panel,
                    "target_h2": args.target_h2,
                    "realized_h2": phenotype_qc["realized_h2"],
                    "estimated_h2": h2,
                    "h2_se": se,
                    "bias": h2 - args.target_h2,
                    "absolute_error": abs(h2 - args.target_h2),
                    "squared_error": (h2 - args.target_h2) ** 2,
                    "ci95_covers_target": int(h2 - 1.96 * se <= args.target_h2 <= h2 + 1.96 * se),
                    "intercept": intercept,
                    "intercept_se": intercept_se,
                }
            )

    summary: list[dict[str, object]] = []
    for panel in PANELS:
        subset = [row for row in rows if row["panel"] == panel]
        n = len(subset)
        summary.append(
            {
                "panel": panel,
                "replicates": n,
                "mean_estimated_h2": sum(float(row["estimated_h2"]) for row in subset) / n,
                "mean_bias": sum(float(row["bias"]) for row in subset) / n,
                "mean_absolute_error": sum(float(row["absolute_error"]) for row in subset) / n,
                "rmse": math.sqrt(sum(float(row["squared_error"]) for row in subset) / n),
                "mean_se": sum(float(row["h2_se"]) for row in subset) / n,
                "coverage_95": sum(int(row["ci95_covers_target"]) for row in subset) / n,
                "mean_intercept": sum(float(row["intercept"]) for row in subset) / n,
            }
        )

    for path, output_rows in (
        (args.root / "results" / "simulation_replicates.tsv", rows),
        (args.root / "results" / "simulation_accuracy_summary.tsv", summary),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(output_rows)


if __name__ == "__main__":
    main()
