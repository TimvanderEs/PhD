#!/usr/bin/env python3
"""Convert a HELIOS component GWAS to signed HapMap3-munging input.

The input may be a gzip-compressed 10k METAL file or a provider ``.tar``
archive containing one 22k TSV. Effects are optionally multiplied by -1 so
that all six traits can be oriented toward higher cognitive performance.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import sys
import tarfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


def first_present(fieldnames: list[str], candidates: list[str]) -> str | None:
    return next((candidate for candidate in candidates if candidate in fieldnames), None)


@contextmanager
def open_input(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", newline="") as handle:
            yield handle
        return
    if path.suffix == ".tar":
        with tarfile.open(path, "r:*") as archive:
            members = [member for member in archive.getmembers() if member.isfile() and member.name.endswith(".tsv")]
            if len(members) != 1:
                raise ValueError(f"Expected one TSV in {path}; found {len(members)}")
            binary_handle = archive.extractfile(members[0])
            if binary_handle is None:
                raise ValueError(f"Could not read {members[0].name}")
            with io.TextIOWrapper(binary_handle, encoding="utf-8", newline="") as handle:
                yield handle
        return
    with path.open(newline="") as handle:
        yield handle


def load_rosetta(path: Path) -> dict[str, tuple[str, str, str]]:
    lookup: dict[str, tuple[str, str, str]] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            value = (row["RSID"], row["CHR"], row["BP_37"])
            prefix = f"chr{row['CHR']}:{row['BP_38']}"
            a1, a2 = row["A1"].upper(), row["A2"].upper()
            lookup[f"{prefix}:{a1}:{a2}"] = value
            lookup[f"{prefix}:{a2}:{a1}"] = value
    return lookup


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--rosetta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sign", required=True, choices=(-1, 1), type=int)
    args = parser.parse_args()

    lookup = load_rosetta(args.rosetta)
    print(f"Loaded {len(lookup):,} oriented marker keys", file=sys.stderr)

    counts = {
        "input": 0,
        "mapped": 0,
        "not_in_rosetta": 0,
        "missing_n": 0,
        "missing_stat": 0,
        "duplicate": 0,
        "bad_parse": 0,
    }
    seen: set[tuple[str, str, str]] = set()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open_input(args.input) as source, gzip.open(args.output, "wt", newline="") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        fields = reader.fieldnames or []
        marker_column = first_present(fields, ["SNP", "ID", "MarkerName"])
        a1_column = first_present(fields, ["A1", "Allele1"])
        a2_column = first_present(fields, ["A2", "Allele2"])
        n_column = first_present(fields, ["N", "N_trans"])
        z_column = first_present(fields, ["Z"])
        beta_column = first_present(fields, ["BETA", "Effect", "BETA_trans"])
        se_column = first_present(fields, ["SE", "StdErr", "SE_trans"])
        p_column = first_present(fields, ["P", "P.value", "P-value", "P_trans"])
        required = [marker_column, a1_column, a2_column, n_column]
        if any(column is None for column in required) or (z_column is None and (beta_column is None or se_column is None)):
            raise SystemExit(f"Unsupported columns in {args.input}: {fields}")

        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "N", "Z", "P"])

        for row in reader:
            counts["input"] += 1
            marker = row.get(marker_column or "", "")
            mapped = lookup.get(marker)
            if mapped is None:
                counts["not_in_rosetta"] += 1
                continue
            rsid, chrom, bp37 = mapped
            a1 = row.get(a1_column or "", "").upper()
            a2 = row.get(a2_column or "", "").upper()
            deduplication_key = (rsid, a1, a2)
            if deduplication_key in seen:
                counts["duplicate"] += 1
                continue
            seen.add(deduplication_key)

            n_value = row.get(n_column or "", "")
            if n_value in {"", "NA", "NaN", "nan"}:
                counts["missing_n"] += 1
                continue
            try:
                sample_size = float(n_value)
                if z_column and row.get(z_column, "") not in {"", "NA", "NaN", "nan"}:
                    z_value = float(row[z_column])
                else:
                    beta = float(row.get(beta_column or "", ""))
                    standard_error = float(row.get(se_column or "", ""))
                    if standard_error == 0:
                        counts["missing_stat"] += 1
                        continue
                    z_value = beta / standard_error
                z_value *= args.sign
                p_raw = row.get(p_column or "", "") if p_column else ""
                p_value = float(p_raw) if p_raw not in {"", "NA", "NaN", "nan"} else math.erfc(abs(z_value) / math.sqrt(2))
            except (TypeError, ValueError, ZeroDivisionError):
                counts["bad_parse"] += 1
                continue
            if p_value == 0:
                p_value = 1e-300
            formatted_n = int(sample_size) if sample_size.is_integer() else sample_size
            writer.writerow([rsid, chrom, bp37, a1, a2, formatted_n, z_value, p_value])
            counts["mapped"] += 1

    for label, value in counts.items():
        print(f"{label}: {value}", file=sys.stderr)
    if counts["mapped"] < 2_000_000:
        raise SystemExit("Too few mapped variants; refusing to continue")


if __name__ == "__main__":
    main()
