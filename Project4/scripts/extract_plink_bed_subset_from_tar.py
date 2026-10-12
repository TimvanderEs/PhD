#!/usr/bin/env python3
"""Stream a SNP subset from a SNP-major PLINK BED stored in a tar.gz archive.

This avoids extracting a very large whole-genome BED when only a prespecified
variant universe is required. The companion BIM and FAM files must describe
the archived BED exactly. No genotype values are decoded or altered.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tarfile
from pathlib import Path


BED_MAGIC = b"\x6c\x1b\x01"
COPY_CHUNK = 16 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--bed-member", required=True)
    parser.add_argument("--bim", required=True, type=Path)
    parser.add_argument("--fam", required=True, type=Path)
    parser.add_argument("--extract", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def count_lines(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def read_exact(handle, size: int) -> bytes:
    chunks = []
    remaining = size
    while remaining:
        chunk = handle.read(remaining)
        if not chunk:
            raise EOFError(f"Unexpected end of BED member with {remaining} bytes pending")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def discard_exact(handle, size: int) -> None:
    remaining = size
    while remaining:
        chunk = handle.read(min(remaining, COPY_CHUNK))
        if not chunk:
            raise EOFError(f"Unexpected end of BED member with {remaining} bytes to skip")
        remaining -= len(chunk)


def main() -> None:
    args = parse_args()
    for path in (args.archive, args.bim, args.fam, args.extract):
        if not path.is_file():
            raise FileNotFoundError(path)

    requested = {
        line.split()[0]
        for line in args.extract.read_text().splitlines()
        if line.strip()
    }
    if not requested:
        raise ValueError("The extraction list is empty")

    selected = []
    duplicate_bim_ids = set()
    seen_bim_ids = set()
    chromosome_counts = {}
    total_variants = 0
    with args.bim.open() as bim_handle:
        for index, line in enumerate(bim_handle):
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"Malformed BIM line {index + 1}")
            chromosome, variant_id = fields[0], fields[1]
            if variant_id in seen_bim_ids:
                duplicate_bim_ids.add(variant_id)
            seen_bim_ids.add(variant_id)
            if variant_id in requested:
                selected.append((index, line))
                chromosome_counts[chromosome] = chromosome_counts.get(chromosome, 0) + 1
            total_variants += 1

    if duplicate_bim_ids.intersection(requested):
        raise ValueError("Requested variant IDs are duplicated in the BIM")
    if not selected:
        raise ValueError("None of the requested variants occur in the BIM")

    sample_count = count_lines(args.fam)
    bytes_per_variant = (sample_count + 3) // 4
    expected_archived_bed_size = 3 + total_variants * bytes_per_variant
    expected_output_bed_size = 3 + len(selected) * bytes_per_variant
    free_bytes = shutil.disk_usage(args.out.parent).free
    if free_bytes < expected_output_bed_size + 2 * 1024**3:
        raise OSError(
            f"Insufficient free space: need at least {expected_output_bed_size + 2 * 1024**3} "
            f"bytes, found {free_bytes}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    bed_partial = args.out.with_suffix(".bed.partial")
    bim_partial = args.out.with_suffix(".bim.partial")
    fam_partial = args.out.with_suffix(".fam.partial")

    with tarfile.open(args.archive, mode="r|gz") as archive:
        bed_member = None
        for member in archive:
            if member.name == args.bed_member:
                bed_member = member
                break
        if bed_member is None:
            raise FileNotFoundError(f"BED member not found: {args.bed_member}")
        if bed_member.size != expected_archived_bed_size:
            raise ValueError(
                f"Archived BED size {bed_member.size} does not match BIM/FAM expectation "
                f"{expected_archived_bed_size}"
            )
        source = archive.extractfile(bed_member)
        if source is None:
            raise OSError("Could not open archived BED member")

        if read_exact(source, 3) != BED_MAGIC:
            raise ValueError("Archived BED is not SNP-major PLINK BED format")

        cursor = 0
        with bed_partial.open("wb") as bed_out:
            bed_out.write(BED_MAGIC)
            for output_index, (source_index, _) in enumerate(selected, start=1):
                discard_exact(source, (source_index - cursor) * bytes_per_variant)
                bed_out.write(read_exact(source, bytes_per_variant))
                cursor = source_index + 1
                if output_index % 100000 == 0:
                    print(
                        f"Copied {output_index}/{len(selected)} variants "
                        f"(source index {source_index})",
                        flush=True,
                    )

    with bim_partial.open("w") as bim_out:
        for _, line in selected:
            bim_out.write(line)
    shutil.copyfile(args.fam, fam_partial)

    if bed_partial.stat().st_size != expected_output_bed_size:
        raise ValueError("Output BED size check failed")
    if count_lines(bim_partial) != len(selected):
        raise ValueError("Output BIM row-count check failed")
    if count_lines(fam_partial) != sample_count:
        raise ValueError("Output FAM row-count check failed")

    os.replace(bed_partial, args.out.with_suffix(".bed"))
    os.replace(bim_partial, args.out.with_suffix(".bim"))
    os.replace(fam_partial, args.out.with_suffix(".fam"))

    report = {
        "archive": str(args.archive),
        "bed_member": args.bed_member,
        "source_bim": str(args.bim),
        "source_fam": str(args.fam),
        "extraction_list": str(args.extract),
        "source_samples": sample_count,
        "source_variants": total_variants,
        "requested_unique_variants": len(requested),
        "matched_variants": len(selected),
        "unmatched_variants": len(requested) - len(selected),
        "bytes_per_variant": bytes_per_variant,
        "source_bed_bytes": expected_archived_bed_size,
        "output_bed_bytes": expected_output_bed_size,
        "chromosome_counts": chromosome_counts,
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
