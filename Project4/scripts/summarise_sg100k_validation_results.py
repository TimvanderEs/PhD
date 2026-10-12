#!/usr/bin/env python3
"""Create compact comparison tables from the SG100K validation outputs."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def value(row: dict[str, str], key: str) -> float | None:
    try:
        number = float(row[key])
    except (KeyError, TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def mean(values: list[float]) -> str:
    return f"{statistics.fmean(values):.6g}" if values else "NA"


def median(values: list[float]) -> str:
    return f"{statistics.median(values):.6g}" if values else "NA"


def pearson(first: list[float], second: list[float]) -> str:
    if len(first) < 2 or len(first) != len(second):
        return "NA"
    first_mean = statistics.fmean(first)
    second_mean = statistics.fmean(second)
    numerator = sum((x - first_mean) * (y - second_mean) for x, y in zip(first, second))
    denominator = math.sqrt(
        sum((x - first_mean) ** 2 for x in first)
        * sum((y - second_mean) ** 2 for y in second)
    )
    return f"{numerator / denominator:.6g}" if denominator else "NA"


def h2_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["gwas_ancestry"], row["ld_panel"])].append(row)

    output: list[dict[str, object]] = []
    for (ancestry, panel), group in sorted(grouped.items()):
        h2 = [number for row in group if (number := value(row, "h2")) is not None]
        se = [number for row in group if (number := value(row, "h2_se")) is not None]
        z = [number for row in group if (number := value(row, "h2_z")) is not None]
        chi2 = [number for row in group if (number := value(row, "mean_chi2")) is not None]
        intercept = [number for row in group if (number := value(row, "intercept")) is not None]
        output.append(
            {
                "gwas_ancestry": ancestry,
                "ld_panel": panel,
                "completed_traits": sum(row["status"] == "complete" for row in group),
                "total_traits": len(group),
                "mean_h2": mean(h2),
                "median_h2": median(h2),
                "mean_h2_se": mean(se),
                "mean_h2_z": mean(z),
                "traits_h2_z_ge_1_96": sum(abs(number) >= 1.96 for number in z),
                "mean_chi2": mean(chi2),
                "mean_intercept": mean(intercept),
            }
        )
    return output


def chinese_h2_comparison(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    index = {(row["trait_id"], row["ld_panel"]): row for row in rows if row["gwas_ancestry"] == "Chinese"}
    output: list[dict[str, object]] = []
    traits = sorted({key[0] for key in index})
    for trait in traits:
        singapore = index.get((trait, "Singapore_Chinese_1cM"))
        public = index.get((trait, "1KG_EAS_1cM"))
        if not singapore or not public:
            continue
        sg_h2 = value(singapore, "h2")
        public_h2 = value(public, "h2")
        difference = sg_h2 - public_h2 if sg_h2 is not None and public_h2 is not None else None
        percent = (
            100 * difference / public_h2
            if difference is not None and public_h2 not in (None, 0)
            else None
        )
        output.append(
            {
                "trait_id": trait,
                "singapore_chinese_h2": singapore["h2"],
                "singapore_chinese_se": singapore["h2_se"],
                "onekg_eas_h2": public["h2"],
                "onekg_eas_se": public["h2_se"],
                "h2_difference_singapore_minus_1kg": f"{difference:.6g}" if difference is not None else "NA",
                "relative_difference_percent": f"{percent:.6g}" if percent is not None else "NA",
                "singapore_intercept": singapore["intercept"],
                "onekg_eas_intercept": public["intercept"],
                "note": "Paired point-estimate comparison; no difference P value because cross-panel estimate covariance is unavailable.",
            }
        )
    return output


def rg_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["gwas_ancestry"], row["ld_panel"])].append(row)

    output: list[dict[str, object]] = []
    for (ancestry, panel), group in sorted(grouped.items()):
        estimates = [number for row in group if (number := value(row, "rg")) is not None]
        ses = [number for row in group if (number := value(row, "rg_se")) is not None]
        output.append(
            {
                "gwas_ancestry": ancestry,
                "ld_panel": panel,
                "completed_pairs": sum(row["status"] == "complete" for row in group),
                "total_pairs": len(group),
                "mean_rg": mean(estimates),
                "median_rg": median(estimates),
                "mean_rg_se": mean(ses),
                "nominal_p_lt_0_05": sum((number := value(row, "p")) is not None and number < 0.05 for row in group),
                "bh_fdr_lt_0_05": sum(row.get("significant_bh_0_05") == "yes" for row in group),
            }
        )
    return output


def chinese_rg_comparison(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    index = {
        (row["trait_1"], row["trait_2"], row["ld_panel"]): row
        for row in rows
        if row["gwas_ancestry"] == "Chinese"
    }
    pairs = sorted({(key[0], key[1]) for key in index})
    singapore: list[float] = []
    public: list[float] = []
    absolute_differences: list[float] = []
    direction_agreement = 0
    for first, second in pairs:
        sg = index.get((first, second, "Singapore_Chinese_1cM"))
        kg = index.get((first, second, "1KG_EAS_1cM"))
        if not sg or not kg:
            continue
        sg_rg = value(sg, "rg")
        kg_rg = value(kg, "rg")
        if sg_rg is None or kg_rg is None:
            continue
        singapore.append(sg_rg)
        public.append(kg_rg)
        absolute_differences.append(abs(sg_rg - kg_rg))
        direction_agreement += (sg_rg >= 0) == (kg_rg >= 0)

    return [
        {
            "gwas_ancestry": "Chinese",
            "panel_1": "Singapore_Chinese_1cM",
            "panel_2": "1KG_EAS_1cM",
            "complete_pairwise_comparisons": len(singapore),
            "pearson_rg": pearson(singapore, public),
            "mean_absolute_rg_difference": mean(absolute_differences),
            "direction_agreement": f"{direction_agreement}/{len(singapore)}",
        }
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tables", type=Path, help="Directory containing ancestry_matched_h2.tsv and ancestry_matched_rg.tsv")
    args = parser.parse_args()
    h2_rows = read_tsv(args.tables / "ancestry_matched_h2.tsv")
    rg_rows = read_tsv(args.tables / "ancestry_matched_rg.tsv")

    h2_grouped = h2_summary(h2_rows)
    write_tsv(args.tables / "h2_summary_by_analysis.tsv", list(h2_grouped[0]), h2_grouped)
    h2_comparison = chinese_h2_comparison(h2_rows)
    write_tsv(args.tables / "chinese_h2_panel_comparison.tsv", list(h2_comparison[0]), h2_comparison)
    rg_grouped = rg_summary(rg_rows)
    write_tsv(args.tables / "rg_summary_by_analysis.tsv", list(rg_grouped[0]), rg_grouped)
    rg_comparison = chinese_rg_comparison(rg_rows)
    write_tsv(args.tables / "chinese_rg_panel_comparison.tsv", list(rg_comparison[0]), rg_comparison)


if __name__ == "__main__":
    main()
