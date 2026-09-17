#!/usr/bin/env python3
"""Apply the retained HELIOS 10k post-METAL variant filters."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import TextIO


MARKER_PATTERN = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]+|X|Y):(?P<pos>[0-9]+)(?::.*)?$", re.I
)


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="METAL output, plain text or gzip")
    parser.add_argument("output", type=Path, help="Filtered output, plain text or gzip")
    parser.add_argument("--heterogeneity-p", type=float, default=0.05)
    parser.add_argument("--max-frequency-difference", type=float, default=0.5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    retained = 0
    with open_text(args.input, "r") as source, open_text(args.output, "w") as destination:
        header_line = source.readline()
        if not header_line:
            raise ValueError("Input is empty")

        header = header_line.split()
        required = {"MarkerName", "MinFreq", "MaxFreq", "HetDf", "HetPVal"}
        missing = required.difference(header)
        if missing:
            raise ValueError(f"Missing METAL columns: {', '.join(sorted(missing))}")

        index = {name: position for position, name in enumerate(header)}
        has_coordinates = "CHR" in index and "BP" in index
        output_header = header if has_coordinates else header + ["CHR", "BP"]
        destination.write("\t".join(output_header) + "\n")

        for line_number, line in enumerate(source, start=2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != len(header):
                raise ValueError(
                    f"Line {line_number} has {len(fields)} fields; expected {len(header)}"
                )
            total += 1

            heterogeneity_df = float(fields[index["HetDf"]])
            heterogeneity_p = float(fields[index["HetPVal"]])
            frequency_difference = (
                float(fields[index["MaxFreq"]]) - float(fields[index["MinFreq"]])
            )
            if not (
                heterogeneity_df >= 1
                and heterogeneity_p > args.heterogeneity_p
                and frequency_difference < args.max_frequency_difference
            ):
                continue

            if not has_coordinates:
                marker = fields[index["MarkerName"]]
                match = MARKER_PATTERN.match(marker)
                if match is None:
                    raise ValueError(
                        f"Cannot parse chromosome and position from marker {marker!r} "
                        f"on line {line_number}"
                    )
                fields.extend([match.group("chrom"), match.group("pos")])

            destination.write("\t".join(fields) + "\n")
            retained += 1

    print(f"Read {total} variants; retained {retained}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
