#!/usr/bin/env python3
"""Restrict LD-score panels to an identical per-chromosome SNP universe."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel",
        action="append",
        required=True,
        metavar="NAME=DIR",
        help="Panel label and directory; repeat for each panel.",
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def read_panel(path: Path) -> tuple[str, list[str], dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n")
        order: list[str] = []
        rows: dict[str, str] = {}
        for line_number, line in enumerate(handle, start=2):
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed LD-score row {line_number}: {path}")
            marker = fields[1]
            if marker in rows:
                raise ValueError(f"Duplicate SNP {marker}: {path}")
            order.append(marker)
            rows[marker] = line if line.endswith("\n") else line + "\n"
    if not rows:
        raise ValueError(f"No LD-score rows: {path}")
    return header, order, rows


def main() -> None:
    args = parse_args()
    panels: list[tuple[str, Path]] = []
    for value in args.panel:
        if "=" not in value:
            raise ValueError(f"Panel must use NAME=DIR syntax: {value}")
        name, directory = value.split("=", 1)
        panels.append((name, Path(directory)))
    if len({name for name, _ in panels}) != len(panels):
        raise ValueError("Panel names must be unique")

    report_rows: list[dict[str, object]] = []
    for chromosome in range(1, 23):
        loaded = []
        for name, directory in panels:
            loaded.append((name, *read_panel(directory / f"{chromosome}.l2.ldscore.gz")))
        common = set(loaded[0][3])
        for _, _, _, rows in loaded[1:]:
            common.intersection_update(rows)
        if not common:
            raise ValueError(f"No common SNPs on chromosome {chromosome}")
        common_order = [marker for marker in loaded[0][2] if marker in common]
        for name, header, _, rows in loaded:
            output = args.output_root / name
            output.mkdir(parents=True, exist_ok=True)
            with gzip.open(output / f"{chromosome}.l2.ldscore.gz", "wt") as handle:
                handle.write(header + "\n")
                for marker in common_order:
                    handle.write(rows[marker])
            for suffix in ("l2.M", "l2.M_5_50"):
                (output / f"{chromosome}.{suffix}").write_text(f"{len(common_order)}\n")
            report_rows.append(
                {
                    "chromosome": chromosome,
                    "panel": name,
                    "source_rows": len(rows),
                    "matched_rows": len(common_order),
                }
            )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)


if __name__ == "__main__":
    main()
