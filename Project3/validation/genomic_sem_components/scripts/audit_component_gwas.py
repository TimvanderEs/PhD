#!/usr/bin/env python3
"""Audit the six 10k and six 22k HELIOS cognitive-component GWAS files.

The script reads source files without modifying them and writes a compact TSV
manifest containing file checksums, trait identities, row/column counts, and
sample-size ranges. Four of the 22k sources are tar archives; the remaining two
can be supplied as either tar archives or plain TSV files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import sys
import tarfile
from pathlib import Path


TENK_HEADER = [
    "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq",
    "MaxFreq", "Effect", "StdErr", "P-value", "Direction", "HetISq",
    "HetChiSq", "HetDf", "HetPVal", "N", "CHR", "BP",
]

TWENTYTWO_HEADER = [
    "CHR", "POS", "ID", "A1", "A2",
    "FREQ_chinese", "BETA_chinese", "SE_chinese", "P_chinese", "N_chinese",
    "FREQ_indian", "BETA_indian", "SE_indian", "P_indian", "N_indian",
    "FREQ_malay", "BETA_malay", "SE_malay", "P_malay", "N_malay",
    "FREQ_trans", "BETA_trans", "SE_trans", "P_trans", "N_trans",
    "HetISq", "HetChiSq", "HetDf", "HetPVal",
]

TRAITS = {
    "dc7r4_pairing7_guesses": "Pairing guesses",
    "dc7r5_quiz_score": "Quiz score",
    "dc7r6_react_avg": "Reaction time",
    "dc7r8_stroopbox_avg": "Stroop box time",
    "dc7r10_stroopink_avg": "Stroop ink time",
    "dc7r12_wm_score": "Working-memory score",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def trait_from_name(name: str) -> tuple[str, str]:
    lowered = name.lower()
    for key, label in TRAITS.items():
        if key in lowered:
            return key, label
    raise ValueError(f"Could not infer trait from {name}")


def audit_rows(handle: io.TextIOBase, expected_header: list[str], n_column: str) -> dict[str, str]:
    reader = csv.reader(handle, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration as exc:
        raise ValueError("Empty data file") from exc
    header_ok = header == expected_header
    if n_column not in header:
        raise ValueError(f"Missing required sample-size column {n_column}")
    n_index = header.index(n_column)
    expected_columns = len(header)
    rows = 0
    bad_width = 0
    n_min_value: float | None = None
    n_max_value: float | None = None
    chr_values: set[str] = set()
    for row in reader:
        rows += 1
        if len(row) != expected_columns:
            bad_width += 1
            continue
        if row[n_index] not in {"", "NA", "NaN", "nan"}:
            n_value = float(row[n_index])
            n_min_value = n_value if n_min_value is None else min(n_min_value, n_value)
            n_max_value = n_value if n_max_value is None else max(n_max_value, n_value)
        chr_values.add(row[0] if header[0] == "CHR" else row[header.index("CHR")])
    n_min = "NA" if n_min_value is None else f"{n_min_value:g}"
    n_max = "NA" if n_max_value is None else f"{n_max_value:g}"
    return {
        "header_ok": str(header_ok).upper(),
        "columns": str(expected_columns),
        "rows": str(rows),
        "bad_width_rows": str(bad_width),
        "n_min": n_min,
        "n_max": n_max,
        "chromosomes": ",".join(sorted(chr_values, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))),
    }


def audit_10k(path: Path) -> dict[str, str]:
    trait_id, trait_label = trait_from_name(path.name)
    print(f"Auditing 10k: {path.name}", file=sys.stderr, flush=True)
    with gzip.open(path, "rt", newline="") as handle:
        stats = audit_rows(handle, TENK_HEADER, "N")
    return {
        "freeze": "10k",
        "trait_id": trait_id,
        "trait": trait_label,
        "source": str(path),
        "container": "gzip",
        "sha256": sha256(path),
        "release": "legacy 10k trans-ancestry meta-analysis",
        "readme_total_n": "NA",
        **stats,
    }


def read_tar_members(path: Path) -> tuple[str, str, dict[str, str]]:
    # The provider archives use gzip compression despite the .tar suffix.
    with tarfile.open(path, "r:*") as archive:
        members = archive.getmembers()
        tsv_members = [member for member in members if member.isfile() and member.name.endswith(".tsv")]
        readme_members = [member for member in members if member.isfile() and member.name.endswith(".readme")]
        if len(tsv_members) != 1 or len(readme_members) != 1:
            raise ValueError(f"Expected one TSV and one readme in {path.name}")
        readme_handle = archive.extractfile(readme_members[0])
        if readme_handle is None:
            raise ValueError(f"Could not read {readme_members[0].name}")
        readme = readme_handle.read().decode("utf-8", errors="replace")
        match = re.search(r"^Total\s*:\s*(\d+)\s*$", readme, flags=re.MULTILINE)
        readme_total = match.group(1) if match else "NA"
        release_match = re.search(r"Genomic Data Version\s*:\s*(.+?)\s*$", readme, flags=re.MULTILINE)
        release = release_match.group(1).strip() if release_match else "NA"
        tsv_handle = archive.extractfile(tsv_members[0])
        if tsv_handle is None:
            raise ValueError(f"Could not read {tsv_members[0].name}")
        with io.TextIOWrapper(tsv_handle, encoding="utf-8", newline="") as text_handle:
            stats = audit_rows(text_handle, TWENTYTWO_HEADER, "N_trans")
    return readme_total, release, stats


def audit_22k_tar(path: Path) -> dict[str, str]:
    trait_id, trait_label = trait_from_name(path.name)
    print(f"Auditing 22k: {path.name}", file=sys.stderr, flush=True)
    readme_total, release, stats = read_tar_members(path)
    return {
        "freeze": "22k",
        "trait_id": trait_id,
        "trait": trait_label,
        "source": str(path),
        "container": "tar",
        "sha256": sha256(path),
        "release": release,
        "readme_total_n": readme_total,
        **stats,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tenk-dir", required=True, type=Path)
    parser.add_argument("--twentytwo-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    tenk_files = sorted(args.tenk_dir.glob("*.gz"))
    twentytwo_files = sorted(args.twentytwo_dir.glob("*.tar"))
    if len(tenk_files) != 6:
        raise SystemExit(f"Expected six 10k gzip files, found {len(tenk_files)}")
    if len(twentytwo_files) != 4:
        raise SystemExit(f"Expected four newly downloaded 22k tar files, found {len(twentytwo_files)}")

    records = [audit_10k(path) for path in tenk_files]
    records.extend(audit_22k_tar(path) for path in twentytwo_files)
    fieldnames = [
        "freeze", "trait_id", "trait", "source", "container", "sha256",
        "release", "readme_total_n", "header_ok", "columns", "rows",
        "bad_width_rows", "n_min", "n_max", "chromosomes",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
