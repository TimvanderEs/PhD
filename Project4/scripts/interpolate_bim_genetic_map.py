#!/usr/bin/env python3
"""Interpolate centimorgan positions from a PLINK-format genetic map."""

from __future__ import annotations

import argparse
import bisect
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", required=True, type=Path)
    parser.add_argument("--map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def interpolate(bp: int, positions: list[int], cms: list[float]) -> tuple[float, str]:
    if bp <= positions[0]:
        return cms[0], "below" if bp < positions[0] else "within"
    if bp >= positions[-1]:
        return cms[-1], "above" if bp > positions[-1] else "within"
    right = bisect.bisect_left(positions, bp)
    if positions[right] == bp:
        return cms[right], "within"
    left = right - 1
    fraction = (bp - positions[left]) / float(positions[right] - positions[left])
    return cms[left] + fraction * (cms[right] - cms[left]), "within"


def main() -> None:
    args = parse_args()
    positions: list[int] = []
    cms: list[float] = []
    map_chromosome = None
    with args.map.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed map line {line_number}")
            chromosome = fields[0].removeprefix("chr")
            if map_chromosome is None:
                map_chromosome = chromosome
            elif chromosome != map_chromosome:
                raise ValueError("Genetic map contains multiple chromosomes")
            positions.append(int(fields[3]))
            cms.append(float(fields[2]))
    if not positions or any(a >= b for a, b in zip(positions, positions[1:])):
        raise ValueError("Genetic-map positions must be non-empty and strictly increasing")
    if any(a > b for a, b in zip(cms, cms[1:])):
        raise ValueError("Genetic-map cM positions must be monotonic")

    counts = {"within": 0, "below": 0, "above": 0}
    previous_bp = -1
    previous_cm = -1.0
    variants = 0
    bim_chromosome = None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.bim.open() as source, args.output.open("w") as destination:
        for line_number, line in enumerate(source, start=1):
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"Malformed BIM line {line_number}")
            chromosome = fields[0].removeprefix("chr")
            if bim_chromosome is None:
                bim_chromosome = chromosome
            elif chromosome != bim_chromosome:
                raise ValueError("BIM contains multiple chromosomes")
            bp = int(fields[3])
            if bp < previous_bp:
                raise ValueError("BIM base-pair positions are not sorted")
            cm, location = interpolate(bp, positions, cms)
            if cm < previous_cm:
                raise ValueError("Interpolated cM positions are not monotonic")
            fields[2] = f"{cm:.8f}"
            destination.write("\t".join(fields[:6]) + "\n")
            counts[location] += 1
            variants += 1
            previous_bp = bp
            previous_cm = cm
    if bim_chromosome != map_chromosome:
        raise ValueError(
            f"BIM chromosome {bim_chromosome} does not match map chromosome {map_chromosome}"
        )

    report = {
        "bim": str(args.bim),
        "map": str(args.map),
        "chromosome": bim_chromosome,
        "variants": variants,
        "map_markers": len(positions),
        "variants_below_map": counts["below"],
        "variants_above_map": counts["above"],
        "minimum_cm": cms[0],
        "maximum_cm": cms[-1],
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
