#!/usr/bin/env python3
"""Parse ancestry-matched SG100K LDSC logs into compact TSV tables."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


TRAITS = [
    "c1_sbp_avg",
    "c11_whr",
    "dbi13_height",
    "c8_bmi",
    "c37_dlab_ldl_mmol_l_log10",
    "c20_t2d",
]
ANCESTRIES = ["Chinese", "Indian", "Malay"]


def number(pattern: str, text: str) -> str:
    match = re.search(pattern, text, re.MULTILINE)
    return match.group(1) if match else "NA"


def pair(pattern: str, text: str) -> tuple[str, str]:
    match = re.search(pattern, text, re.MULTILINE)
    return (match.group(1), match.group(2)) if match else ("NA", "NA")


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def confidence_interval(estimate: str, se: str) -> tuple[str, str]:
    if not finite(estimate) or not finite(se):
        return ("NA", "NA")
    value = float(estimate)
    error = float(se)
    return (f"{value - 1.96 * error:.6g}", f"{value + 1.96 * error:.6g}")


def add_bh_fdr(rows: list[dict[str, str]]) -> None:
    """Correct the 15 trait-pair tests within each ancestry and LD panel."""
    groups: dict[tuple[str, str], list[int]] = {}
    for index, row in enumerate(rows):
        if finite(row["p"]):
            key = (row["gwas_ancestry"], row["ld_panel"])
            groups.setdefault(key, []).append(index)

    for indices in groups.values():
        ordered = sorted(indices, key=lambda index: float(rows[index]["p"]))
        adjusted: dict[int, float] = {}
        running = 1.0
        count = len(ordered)
        for rank_index in range(count - 1, -1, -1):
            original_index = ordered[rank_index]
            rank = rank_index + 1
            candidate = min(1.0, float(rows[original_index]["p"]) * count / rank)
            running = min(running, candidate)
            adjusted[original_index] = running
        for original_index, value in adjusted.items():
            rows[original_index]["p_bh"] = f"{value:.6g}"
            rows[original_index]["significant_bh_0_05"] = "yes" if value < 0.05 else "no"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_h2_log(log: Path, gwas_ancestry: str, ld_panel: str, trait: str) -> dict[str, str]:
    text = log.read_text(errors="replace") if log.exists() else ""
    h2, h2_se = pair(r"Total Observed scale h2:\s+([^\s]+)\s+\(([^)]+)\)", text)
    ci_lower, ci_upper = confidence_interval(h2, h2_se)
    intercept, intercept_se = pair(r"Intercept:\s+([^\s]+)\s+\(([^)]+)\)", text)
    ratio, ratio_se = pair(r"Ratio:\s+([^\s]+)\s+\(([^)]+)\)", text)
    return {
        "gwas_ancestry": gwas_ancestry,
        "ld_panel": ld_panel,
        "trait_id": trait,
        "scale": "observed",
        "status": "complete" if "Analysis finished" in text and finite(h2) and finite(h2_se) else "failed",
        "h2": h2,
        "h2_se": h2_se,
        "h2_ci_lower": ci_lower,
        "h2_ci_upper": ci_upper,
        "h2_z": f"{float(h2) / float(h2_se):.6g}" if finite(h2) and finite(h2_se) and float(h2_se) else "NA",
        "lambda_gc": number(r"Lambda GC:\s+([^\s]+)", text),
        "mean_chi2": number(r"Mean Chi\^2:\s+([^\s]+)", text),
        "intercept": intercept,
        "intercept_se": intercept_se,
        "ratio": ratio,
        "ratio_se": ratio_se,
        "n_snps": number(r"After merging with regression SNP LD,\s+(\d+) SNPs remain", text),
        "log": str(log),
    }


def parse_h2(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for ancestry in ANCESTRIES:
        lower = ancestry.lower()
        for trait in TRAITS:
            log = root / "h2" / ancestry / f"{trait}.{lower}_matched.log"
            rows.append(parse_h2_log(log, ancestry, f"Singapore_{ancestry}_1cM", trait))

    for trait in TRAITS:
        log = root / "h2" / "Chinese_1KGEAS" / f"{trait}.chinese_1kgeas.log"
        if log.exists():
            rows.append(parse_h2_log(log, "Chinese", "1KG_EAS_1cM", trait))
    return rows


def parse_rg(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for ancestry in ANCESTRIES:
        lower = ancestry.lower()
        for i, first in enumerate(TRAITS):
            for second in TRAITS[i + 1 :]:
                log = root / "rg" / ancestry / f"{first}__{second}.{lower}_matched.log"
                text = log.read_text(errors="replace") if log.exists() else ""
                values: list[str] | None = None
                for line in text.splitlines():
                    fields = line.split()
                    if len(fields) >= 12 and fields[0].endswith(".sumstats.gz") and fields[1].endswith(".sumstats.gz"):
                        values = fields
                row = {
                    "gwas_ancestry": ancestry,
                    "ld_panel": f"Singapore_{ancestry}_1cM",
                    "trait_1": first,
                    "trait_2": second,
                    "status": "complete" if "Analysis finished" in text and values else "failed",
                    "rg": "NA",
                    "rg_se": "NA",
                    "z": "NA",
                    "p": "NA",
                    "p_bh": "NA",
                    "significant_bh_0_05": "NA",
                    "h2_obs_2": "NA",
                    "h2_obs_se_2": "NA",
                    "h2_int_2": "NA",
                    "h2_int_se_2": "NA",
                    "gcov_int": "NA",
                    "gcov_int_se": "NA",
                    "log": str(log),
                }
                if values:
                    keys = [
                        "rg",
                        "rg_se",
                        "z",
                        "p",
                        "h2_obs_2",
                        "h2_obs_se_2",
                        "h2_int_2",
                        "h2_int_se_2",
                        "gcov_int",
                        "gcov_int_se",
                    ]
                    row.update(dict(zip(keys, values[2:12])))
                rows.append(row)

    comparator_dir = root / "rg" / "Chinese_1KGEAS"
    if comparator_dir.exists():
        for i, first in enumerate(TRAITS):
            for second in TRAITS[i + 1 :]:
                log = comparator_dir / f"{first}__{second}.chinese_1kgeas.log"
                text = log.read_text(errors="replace") if log.exists() else ""
                values: list[str] | None = None
                for line in text.splitlines():
                    fields = line.split()
                    if len(fields) >= 12 and fields[0].endswith(".sumstats.gz") and fields[1].endswith(".sumstats.gz"):
                        values = fields
                row = {
                    "gwas_ancestry": "Chinese",
                    "ld_panel": "1KG_EAS_1cM",
                    "trait_1": first,
                    "trait_2": second,
                    "status": "complete" if "Analysis finished" in text and values else "failed",
                    "rg": "NA",
                    "rg_se": "NA",
                    "z": "NA",
                    "p": "NA",
                    "p_bh": "NA",
                    "significant_bh_0_05": "NA",
                    "h2_obs_2": "NA",
                    "h2_obs_se_2": "NA",
                    "h2_int_2": "NA",
                    "h2_int_se_2": "NA",
                    "gcov_int": "NA",
                    "gcov_int_se": "NA",
                    "log": str(log),
                }
                if values:
                    keys = [
                        "rg", "rg_se", "z", "p", "h2_obs_2", "h2_obs_se_2",
                        "h2_int_2", "h2_int_se_2", "gcov_int", "gcov_int_se",
                    ]
                    row.update(dict(zip(keys, values[2:12])))
                rows.append(row)
    add_bh_fdr(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    tables = args.root / "tables"
    h2_rows = parse_h2(args.root)
    rg_rows = parse_rg(args.root)
    for row in h2_rows + rg_rows:
        log = Path(row["log"])
        if log.is_absolute():
            row["log"] = str(log.relative_to(args.root))
    write_tsv(tables / "ancestry_matched_h2.tsv", list(h2_rows[0]), h2_rows)
    write_tsv(tables / "ancestry_matched_rg.tsv", list(rg_rows[0]), rg_rows)

    failed_h2 = sum(row["status"] != "complete" for row in h2_rows)
    failed_rg = sum(row["status"] != "complete" for row in rg_rows)
    print(f"h2: {len(h2_rows) - failed_h2}/{len(h2_rows)} complete")
    print(f"rg: {len(rg_rows) - failed_rg}/{len(rg_rows)} complete")
    raise SystemExit(1 if failed_h2 or failed_rg else 0)


if __name__ == "__main__":
    main()
