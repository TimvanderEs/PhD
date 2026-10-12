#!/usr/bin/env python3
"""Summarize real-GWAS h2 estimates across reference-panel sample sizes."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


PANELS = ("1kgeas", "chinese_n500", "chinese_n2500", "chinese_n10000", "chinese_full")
TRAITS = ("c1_sbp_avg", "c11_whr", "dbi13_height", "c8_bmi", "c37_dlab_ldl_mmol_l_log10", "c20_t2d")


def parse_log(path: Path) -> tuple[float, float, float, float, float]:
    text = path.read_text(errors="replace")
    h2 = re.search(r"Total Observed scale h2:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", text)
    intercept = re.search(r"Intercept:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", text)
    mean_chi2 = re.search(r"Mean Chi\^2:\s+([-+0-9.eE]+)", text)
    if h2 is None or intercept is None or mean_chi2 is None:
        raise ValueError(f"Could not parse {path}")
    return float(h2.group(1)), float(h2.group(2)), float(intercept.group(1)), float(intercept.group(2)), float(mean_chi2.group(1))


def read_delete(path: Path) -> list[float]:
    values = [float(line.split()[0]) for line in path.read_text().splitlines() if line.split()]
    if len(values) < 2:
        raise ValueError(f"Too few delete values: {path}")
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    estimates: dict[tuple[str, str], tuple[float, list[float]]] = {}
    rows: list[dict[str, object]] = []
    for panel in PANELS:
        for trait in TRAITS:
            prefix = args.root / "results" / "reference_size_h2" / panel / trait
            h2, se, intercept, intercept_se, mean_chi2 = parse_log(prefix.with_suffix(".log"))
            delete = read_delete(prefix.with_suffix(".delete"))
            estimates[(panel, trait)] = (h2, delete)
            rows.append(
                {
                    "panel": panel,
                    "trait": trait,
                    "h2": h2,
                    "h2_se": se,
                    "h2_z": h2 / se,
                    "intercept": intercept,
                    "intercept_se": intercept_se,
                    "mean_chi2": mean_chi2,
                    "blocks": len(delete),
                }
            )

    paired: list[dict[str, object]] = []
    for panel in PANELS[:-1]:
        for trait in TRAITS:
            h2_a, delete_a = estimates[(panel, trait)]
            h2_b, delete_b = estimates[("chinese_full", trait)]
            if len(delete_a) != len(delete_b):
                raise ValueError("Unequal jackknife block counts")
            blocks = len(delete_a)
            # The logs round full estimates to four decimals, while delete values
            # retain full precision.  Use the reported full-estimate contrast and
            # derive only its covariance-aware SE from paired delete values.
            difference = h2_a - h2_b
            delete_differences = [a - b for a, b in zip(delete_a, delete_b)]
            mean_delete_difference = sum(delete_differences) / blocks
            variance = (blocks - 1) / blocks * sum(
                (value - mean_delete_difference) ** 2 for value in delete_differences
            )
            se = math.sqrt(variance)
            z = difference / se if se else math.nan
            paired.append(
                {
                    "panel": panel,
                    "comparator": "chinese_full",
                    "trait": trait,
                    "h2_difference": difference,
                    "paired_se": se,
                    "z": z,
                    "p": math.erfc(abs(z) / math.sqrt(2)) if math.isfinite(z) else math.nan,
                    "blocks": blocks,
                }
            )

    panel_summary: list[dict[str, object]] = []
    for panel in PANELS:
        subset = [row for row in rows if row["panel"] == panel]
        panel_summary.append(
            {
                "panel": panel,
                "traits": len(subset),
                "mean_h2": sum(float(row["h2"]) for row in subset) / len(subset),
                "mean_h2_se": sum(float(row["h2_se"]) for row in subset) / len(subset),
                "mean_h2_z": sum(float(row["h2_z"]) for row in subset) / len(subset),
                "mean_intercept": sum(float(row["intercept"]) for row in subset) / len(subset),
            }
        )

    for filename, output in (
        ("reference_size_h2.tsv", rows),
        ("reference_size_h2_panel_summary.tsv", panel_summary),
        ("reference_size_paired.tsv", paired),
    ):
        path = args.root / "results" / filename
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(output)


if __name__ == "__main__":
    main()
