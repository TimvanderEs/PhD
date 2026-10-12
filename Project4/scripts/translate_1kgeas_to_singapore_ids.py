#!/usr/bin/env python3
"""Translate 1KG EAS LD-score rsIDs to the Singapore chr:position:alleles IDs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Allow source rsIDs absent from the target-ID crosswalk; all chromosomes must still retain rows.",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    rsid_to_marker: dict[str, str] = {}
    with args.mapping.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split()
            if len(fields) != 2:
                raise ValueError(f"Malformed mapping line {line_number}")
            rsid, marker = fields
            previous = rsid_to_marker.setdefault(rsid, marker)
            if previous != marker:
                raise ValueError(f"Ambiguous rsID mapping: {rsid}")

    count_rows: list[dict[str, int]] = []
    seen_markers: set[str] = set()
    for chromosome in range(1, 23):
        source = args.source / f"{chromosome}.l2.ldscore.gz"
        output = args.output / source.name
        source_rows = 0
        translated_rows = 0
        missing_rows = 0
        duplicate_rows = 0
        with gzip.open(source, "rt") as input_handle, gzip.open(output, "wt") as output_handle:
            header = input_handle.readline()
            if not header:
                raise ValueError(f"Empty source LD-score file: {source}")
            output_handle.write(header)
            for line in input_handle:
                source_rows += 1
                fields = line.rstrip("\n").split()
                if len(fields) < 4:
                    raise ValueError(f"Malformed LD-score row in {source}")
                marker = rsid_to_marker.get(fields[1])
                if marker is None:
                    missing_rows += 1
                    continue
                if marker in seen_markers:
                    duplicate_rows += 1
                    continue
                seen_markers.add(marker)
                fields[1] = marker
                output_handle.write("\t".join(fields) + "\n")
                translated_rows += 1

        for suffix in ("l2.M", "l2.M_5_50"):
            (args.output / f"{chromosome}.{suffix}").write_text(f"{translated_rows}\n")
        count_rows.append(
            {
                "chromosome": chromosome,
                "source_rows": source_rows,
                "translated_rows": translated_rows,
                "missing_mapping": missing_rows,
                "duplicate_marker": duplicate_rows,
            }
        )

    counts = args.output / "id_translation_counts.tsv"
    with counts.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(count_rows[0]))
        writer.writeheader()
        writer.writerows(count_rows)

    checksum_path = args.output / "CHECKSUMS.sha256"
    with checksum_path.open("w") as handle:
        for path in sorted(args.output.glob("*.l2.*")):
            handle.write(f"{sha256(path)}  {path.name}\n")
        handle.write(f"{sha256(counts)}  {counts.name}\n")
        handle.write(f"{sha256(args.mapping)}  source_mapping/{args.mapping.name}\n")

    source_total = sum(row["source_rows"] for row in count_rows)
    translated_total = sum(row["translated_rows"] for row in count_rows)
    missing_total = sum(row["missing_mapping"] for row in count_rows)
    duplicate_total = sum(row["duplicate_marker"] for row in count_rows)
    empty_chromosomes = [row["chromosome"] for row in count_rows if row["translated_rows"] == 0]
    print(
        f"source={source_total} translated={translated_total} "
        f"missing={missing_total} duplicates={duplicate_total}"
    )
    if translated_total == 0 or empty_chromosomes or duplicate_total:
        raise SystemExit(1)
    if missing_total and not args.allow_missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
