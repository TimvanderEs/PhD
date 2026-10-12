#!/usr/bin/env python3
"""Align identical LD-score panels to one authoritative per-chromosome SNP order."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


PANELS = ("1kgeas_matched", "chinese_standard", "chinese_cov20pc", "pooled_standard", "pooled_cov20pc")


def read(path: Path) -> tuple[str, list[str], dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n")
        order: list[str] = []
        rows: dict[str, str] = {}
        for line in handle:
            marker = line.split()[1]
            if marker in rows:
                raise ValueError(f"Duplicate marker {marker}: {path}")
            order.append(marker)
            rows[marker] = line if line.endswith("\n") else line + "\n"
    return header, order, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    panel_root = args.root / "panels"
    report: list[dict[str, object]] = []
    for chromosome in range(1, 23):
        authority_path = panel_root / "chinese_standard" / f"{chromosome}.l2.ldscore.gz"
        _, authority_order, authority_rows = read(authority_path)
        authority_set = set(authority_rows)
        for panel in PANELS:
            path = panel_root / panel / f"{chromosome}.l2.ldscore.gz"
            header, order, rows = read(path)
            if set(rows) != authority_set:
                raise ValueError(f"SNP-set mismatch on chr{chromosome}: {panel}")
            changed = order != authority_order
            if changed:
                temporary = path.with_suffix(path.suffix + ".tmp")
                with gzip.open(temporary, "wt") as handle:
                    handle.write(header + "\n")
                    for marker in authority_order:
                        handle.write(rows[marker])
                temporary.replace(path)
            report.append(
                {
                    "chromosome": chromosome,
                    "panel": panel,
                    "rows": len(order),
                    "order_changed": int(changed),
                }
            )
    report_path = args.root / "qc" / "panel_order_alignment.tsv"
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(report)


if __name__ == "__main__":
    main()
