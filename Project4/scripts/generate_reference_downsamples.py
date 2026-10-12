#!/usr/bin/env python3
"""Create deterministic nested Chinese reference subsets and a held-out GWAS split."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--keep", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[500, 2500, 10000])
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--seed", default="SG100K-reference-size-20260922")
    parser.add_argument("--simulation-reference-size", type=int, default=10000)
    return parser.parse_args()


def order_key(seed: str, replicate: int, fid: str, iid: str) -> bytes:
    value = f"{seed}\t{replicate}\t{fid}\t{iid}".encode()
    return hashlib.sha256(value).digest()


def write_keep(path: Path, rows: list[tuple[str, str]]) -> None:
    path.write_text("".join(f"{fid}\t{iid}\n" for fid, iid in rows))


def main() -> None:
    args = parse_args()
    rows: list[tuple[str, str]] = []
    with args.keep.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Malformed keep row {line_number}: {args.keep}")
            rows.append((fields[0], fields[1]))
    if not rows or len(set(rows)) != len(rows):
        raise ValueError("Input keep file is empty or contains duplicate sample IDs")
    if max(args.sizes + [args.simulation_reference_size]) >= len(rows):
        raise ValueError("Requested subset must be smaller than the full unrelated sample")

    args.output.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    first_order: list[tuple[str, str]] | None = None
    for replicate in range(1, args.replicates + 1):
        ordered = sorted(rows, key=lambda row: order_key(args.seed, replicate, *row))
        if replicate == 1:
            first_order = ordered
        previous: set[tuple[str, str]] = set()
        for size in sorted(set(args.sizes)):
            selected = ordered[:size]
            selected_set = set(selected)
            if previous and not previous.issubset(selected_set):
                raise AssertionError("Subsets are not nested within replicate")
            previous = selected_set
            name = f"chinese_n{size}_r{replicate}.keep"
            write_keep(args.output / name, selected)
            manifest.append(
                {
                    "sample_set": name.removesuffix(".keep"),
                    "sample_size": size,
                    "replicate": replicate,
                    "purpose": "reference_size_sensitivity",
                }
            )

    assert first_order is not None
    simulation_reference = first_order[: args.simulation_reference_size]
    reference_ids = set(simulation_reference)
    simulation_gwas = [row for row in rows if row not in reference_ids]
    write_keep(args.output / "simulation_reference.keep", simulation_reference)
    write_keep(args.output / "simulation_gwas.keep", simulation_gwas)
    manifest.extend(
        [
            {
                "sample_set": "simulation_reference",
                "sample_size": len(simulation_reference),
                "replicate": 1,
                "purpose": "disjoint_simulation_reference",
            },
            {
                "sample_set": "simulation_gwas",
                "sample_size": len(simulation_gwas),
                "replicate": 1,
                "purpose": "held_out_simulation_gwas",
            },
        ]
    )

    with (args.output / "downsample_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    overlap = reference_ids.intersection(simulation_gwas)
    (args.output / "simulation_split_qc.tsv").write_text(
        "metric\tvalue\n"
        f"full_unrelated_chinese\t{len(rows)}\n"
        f"reference_samples\t{len(simulation_reference)}\n"
        f"held_out_gwas_samples\t{len(simulation_gwas)}\n"
        f"reference_gwas_overlap\t{len(overlap)}\n"
    )


if __name__ == "__main__":
    main()
