#!/usr/bin/env python3
"""Validate the compact controlled-reference result bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--reader-results", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    results = args.root / "results"
    qc = args.root / "qc"
    expected_rows = {
        "controlled_h2.tsv": 30,
        "controlled_h2_panel_summary.tsv": 5,
        "paired_h2_differences.tsv": 36,
        "controlled_rg.tsv": 75,
        "controlled_rg_concordance.tsv": 4,
        "reference_size_h2.tsv": 30,
        "reference_size_h2_panel_summary.tsv": 5,
        "reference_size_paired.tsv": 24,
        "reference_size_ldscore_summary.tsv": 3,
        "simulation_replicates.tsv": 15,
        "simulation_accuracy_summary.tsv": 3,
    }
    observed_rows: dict[str, int] = {}
    for filename, expected in expected_rows.items():
        path = results / filename
        require(path.is_file(), f"Missing result table: {path}")
        observed_rows[filename] = len(read_tsv(path))
        require(observed_rows[filename] == expected, f"Unexpected row count for {filename}")

    statuses = {
        "controlled_ldsc_status.tsv": 105,
        "reference_size_h2_status.tsv": 30,
        "simulation_ldsc_status.tsv": 15,
    }
    status_rows: dict[str, int] = {}
    for filename, expected in statuses.items():
        rows = read_tsv(qc / filename)
        status_rows[filename] = len(rows)
        require(len(rows) == expected, f"Unexpected status-row count for {filename}")
        require(all(row["status"] == "complete" for row in rows), f"Incomplete run in {filename}")

    intersection = read_tsv(qc / "final_five_panel_intersection.tsv")
    require(len(intersection) == 110, "Final intersection must contain 22 chromosomes x 5 panels")
    total_variants = 0
    for chromosome in range(1, 23):
        rows = [row for row in intersection if int(row["chromosome"]) == chromosome]
        require(len(rows) == 5, f"Missing final panel row for chromosome {chromosome}")
        counts = {int(row["matched_rows"]) for row in rows}
        require(len(counts) == 1, f"Panel row mismatch on chromosome {chromosome}")
        total_variants += counts.pop()
    require(total_variants == 874025, f"Unexpected shared-variant total: {total_variants}")

    paired = read_tsv(results / "paired_h2_differences.tsv")
    require("jackknife_h2_difference" not in paired[0], "Obsolete rounded pseudo-value column retained")
    require(all(float(row["paired_se"]) > 0 for row in paired), "Non-positive paired SE")

    copied_results: dict[str, bool] = {}
    if args.reader_results is not None:
        for filename in expected_rows:
            copied = args.reader_results / filename
            require(copied.is_file(), f"Missing reader-facing copy: {copied}")
            copied_results[filename] = digest(results / filename) == digest(copied)
            require(copied_results[filename], f"Reader-facing copy differs: {filename}")

    completion = qc / "controlled_experiment_complete.timestamp"
    require(completion.is_file() and completion.read_text().strip(), "Missing completion timestamp")
    report = {
        "status": "pass",
        "shared_variants": total_variants,
        "result_rows": observed_rows,
        "status_rows": status_rows,
        "reader_copies_identical": copied_results,
        "completion_timestamp": completion.read_text().strip(),
    }
    if args.report is not None:
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
