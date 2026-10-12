#!/usr/bin/env python3
"""Summarize controlled LDSC h2 results and paired block-jackknife contrasts."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


PANELS = (
    "1kgeas_matched",
    "chinese_standard",
    "chinese_cov20pc",
    "pooled_standard",
    "pooled_cov20pc",
)
TRAITS = (
    "c1_sbp_avg",
    "c11_whr",
    "dbi13_height",
    "c8_bmi",
    "c37_dlab_ldl_mmol_l_log10",
    "c20_t2d",
)
CONTRASTS = (
    ("chinese_standard", "1kgeas_matched", "Chinese Singapore versus 1000 Genomes EAS reference"),
    ("chinese_cov20pc", "1kgeas_matched", "PC-adjusted Chinese Singapore versus 1000 Genomes EAS reference"),
    ("chinese_cov20pc", "chinese_standard", "PC adjustment within Chinese reference"),
    ("pooled_cov20pc", "pooled_standard", "PC adjustment within pooled reference"),
    ("pooled_standard", "chinese_standard", "Pooled versus Chinese reference without PC adjustment"),
    ("pooled_cov20pc", "chinese_cov20pc", "Pooled versus Chinese reference with PC adjustment"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def parse_log(path: Path) -> dict[str, float]:
    text = path.read_text(errors="replace")
    patterns = {
        "h2": r"Total Observed scale h2:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)",
        "intercept": r"Intercept:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)",
        "mean_chi2": r"Mean Chi\^2:\s+([-+0-9.eE]+)",
        "lambda_gc": r"Lambda GC:\s+([-+0-9.eE]+)",
    }
    h2_match = re.search(patterns["h2"], text)
    intercept_match = re.search(patterns["intercept"], text)
    if h2_match is None or intercept_match is None:
        raise ValueError(f"Could not parse LDSC log: {path}")
    result = {
        "h2": float(h2_match.group(1)),
        "h2_se": float(h2_match.group(2)),
        "intercept": float(intercept_match.group(1)),
        "intercept_se": float(intercept_match.group(2)),
    }
    for key in ("mean_chi2", "lambda_gc"):
        match = re.search(patterns[key], text)
        result[key] = float(match.group(1)) if match else math.nan
    return result


def read_delete_values(path: Path) -> list[float]:
    values = []
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if fields:
                values.append(float(fields[0]))
    if len(values) < 2:
        raise ValueError(f"Too few delete values: {path}")
    return values


def bh_adjust(rows: list[dict], p_key: str, output_key: str) -> None:
    order = sorted(range(len(rows)), key=lambda i: rows[i][p_key])
    adjusted = [1.0] * len(rows)
    running = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        row_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, rows[row_index][p_key] * len(rows) / rank)
        adjusted[row_index] = min(1.0, running)
    for row, value in zip(rows, adjusted):
        row[output_key] = value


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    estimates: dict[tuple[str, str], dict] = {}
    summary_rows = []
    for panel in PANELS:
        for trait in TRAITS:
            prefix = args.root / "results" / "h2" / panel / trait
            log_path = prefix.with_suffix(".log")
            delete_path = prefix.with_suffix(".delete")
            estimate = parse_log(log_path)
            delete_values = read_delete_values(delete_path)
            estimate["delete_values"] = delete_values
            estimates[(panel, trait)] = estimate
            summary_rows.append(
                {
                    "panel": panel,
                    "trait": trait,
                    "h2": estimate["h2"],
                    "h2_se": estimate["h2_se"],
                    "h2_z": estimate["h2"] / estimate["h2_se"],
                    "intercept": estimate["intercept"],
                    "intercept_se": estimate["intercept_se"],
                    "mean_chi2": estimate["mean_chi2"],
                    "lambda_gc": estimate["lambda_gc"],
                    "jackknife_blocks": len(delete_values),
                    "log": str(log_path),
                }
            )

    paired_rows = []
    for panel_a, panel_b, label in CONTRASTS:
        contrast_rows = []
        for trait in TRAITS:
            a = estimates[(panel_a, trait)]
            b = estimates[(panel_b, trait)]
            if len(a["delete_values"]) != len(b["delete_values"]):
                raise ValueError(f"Unequal jackknife block count for {trait}: {panel_a} vs {panel_b}")
            blocks = len(a["delete_values"])
            # LDSC prints the full estimate to four decimals but writes delete-one-
            # block estimates at full precision.  Bias-correcting the rounded full
            # estimate would multiply its rounding error by the block count.  Use
            # the reported full-estimate contrast as the point estimate and the
            # paired delete-value differences only for its jackknife variance.
            difference = a["h2"] - b["h2"]
            delete_differences = [x - y for x, y in zip(a["delete_values"], b["delete_values"])]
            mean_delete_difference = sum(delete_differences) / blocks
            variance = (blocks - 1) / blocks * sum(
                (value - mean_delete_difference) ** 2 for value in delete_differences
            )
            se = math.sqrt(variance)
            z = difference / se if se else math.nan
            p = math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else math.nan
            contrast_rows.append(
                {
                    "contrast": label,
                    "panel_a": panel_a,
                    "panel_b": panel_b,
                    "trait": trait,
                    "h2_difference": difference,
                    "paired_se": se,
                    "z": z,
                    "p": p,
                    "blocks": blocks,
                }
            )
        bh_adjust(contrast_rows, "p", "p_bh_within_contrast")
        paired_rows.extend(contrast_rows)

    write_tsv(
        args.root / "results" / "controlled_h2.tsv",
        summary_rows,
        [
            "panel",
            "trait",
            "h2",
            "h2_se",
            "h2_z",
            "intercept",
            "intercept_se",
            "mean_chi2",
            "lambda_gc",
            "jackknife_blocks",
            "log",
        ],
    )
    panel_summary_rows = []
    for panel in PANELS:
        subset = [row for row in summary_rows if row["panel"] == panel]
        panel_summary_rows.append(
            {
                "panel": panel,
                "traits": len(subset),
                "mean_h2": sum(float(row["h2"]) for row in subset) / len(subset),
                "mean_h2_se": sum(float(row["h2_se"]) for row in subset) / len(subset),
                "mean_h2_z": sum(float(row["h2_z"]) for row in subset) / len(subset),
                "mean_intercept": sum(float(row["intercept"]) for row in subset) / len(subset),
            }
        )
    write_tsv(
        args.root / "results" / "controlled_h2_panel_summary.tsv",
        panel_summary_rows,
        ["panel", "traits", "mean_h2", "mean_h2_se", "mean_h2_z", "mean_intercept"],
    )
    write_tsv(
        args.root / "results" / "paired_h2_differences.tsv",
        paired_rows,
        [
            "contrast",
            "panel_a",
            "panel_b",
            "trait",
            "h2_difference",
            "paired_se",
            "z",
            "p",
            "p_bh_within_contrast",
            "blocks",
        ],
    )


if __name__ == "__main__":
    main()
