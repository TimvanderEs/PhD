#!/usr/bin/env python3
"""Recalculate the aggregate Project 4 LDSC summaries from recovered tables."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


PANEL_LABELS = {
    "1KG_EAS_1cM_old": "1000 Genomes EAS 1 cM",
    "HELIOS_Chinese_1cM_new": "Singapore Chinese 1 cM",
    "HELIOS_Indian_1cM_new": "Singapore Indian 1 cM",
    "HELIOS_Malay_1cM_new": "Singapore Malay 1 cM",
    "Singapore_20cM_covLDSC_old": "Pooled Singapore cov-LDSC 20 cM",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def pearson(x: list[float], y: list[float]) -> float:
    mx, my = mean(x), mean(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator = math.sqrt(
        sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)
    )
    return numerator / denominator


def summarise_h2(rows: list[dict[str, str]]) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["panel"]].append(row)

    print("HERITABILITY SUMMARY")
    print("panel\tn_traits\tmean_h2\tmean_h2_se\tmean_h2_over_se\tmean_intercept")
    for panel in PANEL_LABELS:
        panel_rows = grouped.get(panel, [])
        if not panel_rows:
            continue
        h2 = [float(row["h2"]) for row in panel_rows]
        se = [float(row["h2_se"]) for row in panel_rows]
        intercept = [float(row["intercept"]) for row in panel_rows]
        z = [estimate / error for estimate, error in zip(h2, se)]
        print(
            f"{PANEL_LABELS[panel]}\t{len(panel_rows)}\t{mean(h2):.5f}\t"
            f"{mean(se):.5f}\t{mean(z):.4f}\t{mean(intercept):.5f}"
        )


def summarise_rg(rows: list[dict[str, str]]) -> None:
    by_panel: dict[str, dict[str, float]] = defaultdict(dict)
    nominal_by_pair: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["status"] != "complete":
            continue
        panel = row["ancestry_ld"]
        pair_id = row["pair_id"]
        by_panel[panel][pair_id] = float(row["rg"])
        if row["significant_p05"] == "yes":
            nominal_by_pair[pair_id].add(panel)

    print("\nGENETIC-CORRELATION CONCORDANCE")
    print("panel_1\tpanel_2\tn_pairs\tpearson_r\tsame_direction_n")
    comparisons = [("Chinese", "Indian"), ("Chinese", "Malay"), ("Indian", "Malay")]
    for first, second in comparisons:
        common = sorted(set(by_panel[first]) & set(by_panel[second]))
        x = [by_panel[first][pair] for pair in common]
        y = [by_panel[second][pair] for pair in common]
        same_direction = sum((a >= 0) == (b >= 0) for a, b in zip(x, y))
        print(
            f"{first}\t{second}\t{len(common)}\t{pearson(x, y):.4f}\t{same_direction}"
        )

    complete_panels = set(by_panel)
    common_nominal = sorted(
        pair for pair, panels in nominal_by_pair.items() if panels == complete_panels
    )
    print("\nNOMINALLY SIGNIFICANT IN ALL THREE PANELS")
    for pair in common_nominal:
        print(pair)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h2",
        type=Path,
        default=Path("Project4/results/recovered/helios_h2_panel_comparison_long.tsv"),
    )
    parser.add_argument(
        "--rg",
        type=Path,
        default=Path("Project4/results/recovered/helios_rg_summary.tsv"),
    )
    args = parser.parse_args()
    summarise_h2(read_tsv(args.h2))
    summarise_rg(read_tsv(args.rg))


if __name__ == "__main__":
    main()
