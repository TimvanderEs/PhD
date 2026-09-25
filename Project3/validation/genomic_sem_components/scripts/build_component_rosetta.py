#!/usr/bin/env python3
"""Build the allele-aware HELIOS hg38-to-rsID/hg19 bridge.

This combines the retained SG100K marker set with the existing rsID-to-marker
map and the verified hg19/hg38 coordinate bridge. The output schema is the one
consumed by ``prepare_component_sumstats.py``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rsid-to-marker", required=True, type=Path)
    parser.add_argument("--retained-markers", required=True, type=Path)
    parser.add_argument("--build-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rsid_by_marker: dict[str, str] = {}
    with args.rsid_to_marker.open() as handle:
        for line in handle:
            rsid, marker = line.rstrip("\n").split(maxsplit=1)
            rsid_by_marker[marker] = rsid
    print(f"Loaded {len(rsid_by_marker):,} marker-to-rsID records", file=sys.stderr)

    selected_rsids = set(rsid_by_marker.values())
    coordinates: dict[str, tuple[str, str, str]] = {}
    with args.build_map.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"CHR", "BP", "BP_38", "SNP"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"Build map lacks required columns: {sorted(required)}")
        for row in reader:
            if row["SNP"] in selected_rsids:
                # This recovered legacy file has misleading headers: direct
                # checks against known variants show that BP is the original
                # hg38 coordinate and BP_38 is the lifted hg19 coordinate.
                coordinates[row["SNP"]] = (row["CHR"], row["BP_38"], row["BP"])
    print(f"Recovered coordinates for {len(coordinates):,} selected rsIDs", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    missing_rsid = 0
    missing_coordinate = 0
    coordinate_mismatch = 0
    with args.retained_markers.open() as retained, gzip.open(args.output, "wt", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["RSID", "CHR", "BP_37", "BP_38", "A1", "A2"])
        for line in retained:
            fields = line.rstrip("\n").split()
            if len(fields) != 2:
                continue
            _, marker = fields
            rsid = rsid_by_marker.get(marker)
            if rsid is None:
                missing_rsid += 1
                continue
            marker_fields = marker.split(":")
            if len(marker_fields) != 4:
                missing_rsid += 1
                continue
            marker_chrom = marker_fields[0].removeprefix("chr")
            marker_bp38, a1, a2 = marker_fields[1:]
            coordinate = coordinates.get(rsid)
            if coordinate is None:
                # rsID and alleles remain valid for HapMap3 munging. BP is not
                # consumed by munge_sumstats; retain the marker with unknown
                # hg19 position rather than discarding it silently.
                missing_coordinate += 1
                bp37 = "NA"
            else:
                chrom, bp37, mapped_bp38 = coordinate
                if marker_chrom != chrom or marker_bp38 != mapped_bp38:
                    # A small set of legacy rsIDs points to a different build38
                    # coordinate in the broad bridge. Preserve the allele-aware
                    # marker/rsID match and flag hg19 position as unavailable.
                    coordinate_mismatch += 1
                    bp37 = "NA"
            writer.writerow([rsid, marker_chrom, bp37, marker_bp38, a1.upper(), a2.upper()])
            written += 1

    print(f"Written: {written:,}", file=sys.stderr)
    print(f"Missing rsID: {missing_rsid:,}", file=sys.stderr)
    print(f"Missing coordinate: {missing_coordinate:,}", file=sys.stderr)
    print(f"Coordinate mismatch: {coordinate_mismatch:,}", file=sys.stderr)
    if written < 2_600_000 or (missing_coordinate + coordinate_mismatch) > 5_000:
        raise SystemExit("Rosetta construction failed validation")


if __name__ == "__main__":
    main()
