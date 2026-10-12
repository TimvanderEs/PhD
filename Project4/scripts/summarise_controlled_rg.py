#!/usr/bin/env python3
"""Summarise controlled genetic-correlation runs and panel concordance."""

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


def parse_log(path: Path) -> dict[str, float]:
    text = path.read_text(errors="replace")
    rg = re.search(
        r"Genetic Correlation:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)"
        r"\s+Z-score:\s+([-+0-9.eE]+)\s+P:\s+([-+0-9.eE]+)",
        text,
    )
    gencov = re.search(r"Total Observed scale gencov:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", text)
    h2 = re.findall(r"Total Observed scale h2:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", text)
    if rg is None or gencov is None or len(h2) != 2:
        raise ValueError(f"Could not parse {path}")
    return {
        "rg": float(rg.group(1)),
        "rg_se": float(rg.group(2)),
        "z": float(rg.group(3)),
        "p": float(rg.group(4)),
        "gencov": float(gencov.group(1)),
        "gencov_se": float(gencov.group(2)),
        "h2_1": float(h2[0][0]),
        "h2_1_se": float(h2[0][1]),
        "h2_2": float(h2[1][0]),
        "h2_2_se": float(h2[1][1]),
    }


def pearson(x: list[float], y: list[float]) -> float:
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator = math.sqrt(
        sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)
    )
    return numerator / denominator


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    by_panel: dict[str, dict[str, dict[str, float]]] = {}
    for panel in PANELS:
        panel_rows: dict[str, dict[str, float]] = {}
        for path in sorted((args.root / "results" / "rg" / panel).glob("*.log")):
            pair = path.stem
            trait_1, trait_2 = pair.split("__", 1)
            result = parse_log(path)
            panel_rows[pair] = result
            rows.append({"panel": panel, "trait_1": trait_1, "trait_2": trait_2, **result})
        if len(panel_rows) != 15:
            raise ValueError(f"Expected 15 rg logs for {panel}, found {len(panel_rows)}")
        by_panel[panel] = panel_rows

    reference = by_panel["1kgeas_matched"]
    concordance: list[dict[str, object]] = []
    for panel in PANELS[1:]:
        common = sorted(set(reference) & set(by_panel[panel]))
        reference_rg = [reference[pair]["rg"] for pair in common]
        panel_rg = [by_panel[panel][pair]["rg"] for pair in common]
        differences = [value - comparator for value, comparator in zip(panel_rg, reference_rg)]
        concordance.append(
            {
                "panel": panel,
                "comparator": "1kgeas_matched",
                "pairs": len(common),
                "pearson_rg": pearson(panel_rg, reference_rg),
                "mean_rg_difference": sum(differences) / len(differences),
                "mean_absolute_rg_difference": sum(abs(value) for value in differences) / len(differences),
                "max_absolute_rg_difference": max(abs(value) for value in differences),
                "direction_concordant": sum(
                    (value >= 0) == (comparator >= 0)
                    for value, comparator in zip(panel_rg, reference_rg)
                ),
            }
        )

    write_tsv(args.root / "results" / "controlled_rg.tsv", rows)
    write_tsv(args.root / "results" / "controlled_rg_concordance.tsv", concordance)


if __name__ == "__main__":
    main()
