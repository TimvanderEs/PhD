#!/usr/bin/env python3
"""Map the controlled Singapore SNP set back to raw 1000 Genomes EAS rsIDs."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--source-ldscore-root", type=Path, required=True)
    parser.add_argument("--panel-work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--chromosomes", type=int, nargs="+", default=list(range(1, 23)))
    args = parser.parse_args()

    rsid_to_marker: dict[str, str] = {}
    with args.crosswalk.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(f"Malformed crosswalk line {line_number}")
            rsid, marker = fields
            previous = rsid_to_marker.setdefault(rsid, marker)
            if previous != marker:
                raise ValueError(f"Ambiguous mapping for {rsid}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, int]] = []
    seen_markers_global: set[str] = set()
    for chromosome in args.chromosomes:
        if chromosome < 1 or chromosome > 22:
            raise ValueError(f"Invalid autosome: {chromosome}")
        shared_path = args.panel_work_root / f"chr{chromosome}" / "shared.snplist"
        shared = {line.split()[0] for line in shared_path.open() if line.split()}
        mapping: dict[str, str] = {}
        source_path = args.source_ldscore_root / f"{chromosome}.l2.ldscore.gz"
        with gzip.open(source_path, "rt") as handle:
            header = handle.readline()
            if not header:
                raise ValueError(f"Empty source LD-score file: {source_path}")
            for line in handle:
                fields = line.split()
                if len(fields) < 2:
                    continue
                rsid = fields[1]
                marker = rsid_to_marker.get(rsid)
                if marker in shared and marker not in mapping and marker not in seen_markers_global:
                    mapping[marker] = rsid

        missing = sorted(shared.difference(mapping))
        if missing:
            preview = ", ".join(missing[:5])
            raise ValueError(
                f"Chromosome {chromosome}: {len(missing)} shared markers lack a unique source rsID: {preview}"
            )
        seen_markers_global.update(mapping)

        extract_path = args.output_root / f"chr{chromosome}.extract"
        rename_path = args.output_root / f"chr{chromosome}.rsid_to_marker.tsv"
        with extract_path.open("w") as extract, rename_path.open("w") as rename:
            for marker in sorted(mapping):
                rsid = mapping[marker]
                extract.write(rsid + "\n")
                rename.write(f"{rsid}\t{marker}\n")
        rows.append(
            {
                "chromosome": chromosome,
                "shared_markers": len(shared),
                "mapped_rsids": len(mapping),
                "missing": 0,
            }
        )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
