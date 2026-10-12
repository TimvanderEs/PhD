#!/usr/bin/env python3
"""Rename raw 1KG EAS variants to Singapore IDs and update to GRCh38 positions."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


MARKER = re.compile(r"^chr([0-9]+):([0-9]+):[ACGT]+:[ACGT]+$")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    mapping: dict[str, str] = {}
    with args.mapping.open() as handle:
        for line in handle:
            rsid, marker = line.split()
            mapping[rsid] = marker

    written = 0
    seen: set[str] = set()
    with args.bim.open() as source, args.output.open("w") as target:
        for line_number, line in enumerate(source, start=1):
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"Malformed BIM line {line_number}: {args.bim}")
            marker = mapping.get(fields[1])
            if marker is None:
                raise ValueError(f"Extracted 1KG rsID absent from mapping: {fields[1]}")
            match = MARKER.match(marker)
            if match is None:
                raise ValueError(f"Unexpected Singapore marker format: {marker}")
            if marker in seen:
                raise ValueError(f"Duplicate target marker: {marker}")
            seen.add(marker)
            fields[0] = match.group(1)
            fields[1] = marker
            fields[2] = "0"
            fields[3] = match.group(2)
            target.write("\t".join(fields) + "\n")
            written += 1

    # The precomputed public LD-score panel and the archived raw genotype file
    # are not guaranteed to contain exactly the same release-specific rsIDs.
    # Missing requested rsIDs are handled by a final five-panel intersection;
    # every row actually extracted from the raw BED must still map uniquely.
    print(f"requested={len(mapping)} extracted={written} missing={len(mapping) - written}")


if __name__ == "__main__":
    main()
