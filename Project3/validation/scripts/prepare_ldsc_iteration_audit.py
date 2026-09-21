#!/usr/bin/env python3
"""Prepare checksum-linked HELIOS GWAS freezes for matched LDSC sensitivity runs.

The script streams METAL summary statistics, aligns effect directions to a
frozen LD-score allele list, removes strand-ambiguous variants, and writes
standard LDSC ``SNP N Z A1 A2`` inputs under three explicitly named QC rules.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
from pathlib import Path


COMPLEMENT = str.maketrans("ACGT", "TGCA")
AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def column(header: list[str], *names: str) -> int | None:
    normalized = {name.lower().replace("-", "").replace(".", ""): i for i, name in enumerate(header)}
    for name in names:
        key = name.lower().replace("-", "").replace(".", "")
        if key in normalized:
            return normalized[key]
    return None


def load_reference(path: Path) -> dict[str, tuple[str, str]]:
    reference: dict[str, tuple[str, str]] = {}
    with path.open() as stream:
        header = stream.readline().split()
        snp_i, a1_i, a2_i = (header.index(name) for name in ("SNP", "A1", "A2"))
        for line in stream:
            fields = line.split()
            if len(fields) <= max(snp_i, a1_i, a2_i):
                continue
            reference[fields[snp_i]] = (fields[a1_i].upper(), fields[a2_i].upper())
    return reference


def aligned_z(beta: float, se: float, a1: str, a2: str, ref1: str, ref2: str) -> float | None:
    if se <= 0 or not math.isfinite(beta) or not math.isfinite(se):
        return None
    z = beta / se
    if (a1, a2) == (ref1, ref2):
        return z
    if (a1, a2) == (ref2, ref1):
        return -z
    ca1, ca2 = a1.translate(COMPLEMENT), a2.translate(COMPLEMENT)
    if (ca1, ca2) == (ref1, ref2):
        return z
    if (ca1, ca2) == (ref2, ref1):
        return -z
    return None


def prepare(label: str, path: Path, reference: dict[str, tuple[str, str]], output: Path) -> dict[str, object]:
    destinations = {
        "reference_maf01": output / f"{label}.reference_maf01.sumstats.gz",
        "archived_qc_afdiff050": output / f"{label}.archived_qc_afdiff050.sumstats.gz",
        "strict_qc_afdiff005": output / f"{label}.strict_qc_afdiff005.sumstats.gz",
    }
    handles = {name: gzip.open(path, "wt") for name, path in destinations.items()}
    for stream in handles.values():
        stream.write("SNP\tN\tZ\tA1\tA2\n")

    counts = {
        "input": 0,
        "reference_overlap": 0,
        "valid_aligned": 0,
        "ambiguous": 0,
        "allele_mismatch": 0,
        "maf_fail": 0,
        "archived_qc": 0,
        "strict_qc": 0,
    }
    seen: set[str] = set()
    duplicate_count = 0

    try:
        with open_text(path) as stream:
            header = stream.readline().split()
            indices = {
                "snp": column(header, "SNP", "MarkerName", "ID"),
                "a1": column(header, "A1", "Allele1"),
                "a2": column(header, "A2", "Allele2", "Allele0"),
                "beta": column(header, "BETA", "Effect"),
                "se": column(header, "SE", "StdErr"),
                "n": column(header, "N"),
                "maf": column(header, "MAF", "Freq1", "A1FREQ"),
                "direction": column(header, "Direction"),
                "het_df": column(header, "HetDf"),
                "het_p": column(header, "HetPVal"),
                "min_freq": column(header, "MinFreq"),
                "max_freq": column(header, "MaxFreq"),
            }
            required = ("snp", "a1", "a2", "beta", "se", "n", "maf")
            missing = [name for name in required if indices[name] is None]
            if missing:
                raise ValueError(f"{path}: missing columns {missing}; header={header}")
            maximum = max(index for index in indices.values() if index is not None)

            for line in stream:
                fields = line.split()
                if len(fields) <= maximum:
                    continue
                counts["input"] += 1
                snp = fields[indices["snp"]]  # type: ignore[index]
                if snp in seen:
                    duplicate_count += 1
                    continue
                seen.add(snp)
                if snp not in reference:
                    continue
                counts["reference_overlap"] += 1

                a1 = fields[indices["a1"]].upper()  # type: ignore[index]
                a2 = fields[indices["a2"]].upper()  # type: ignore[index]
                if a1 not in "ACGT" or a2 not in "ACGT" or frozenset((a1, a2)) in AMBIGUOUS:
                    counts["ambiguous"] += 1
                    continue
                try:
                    beta = float(fields[indices["beta"]])  # type: ignore[index]
                    se = float(fields[indices["se"]])  # type: ignore[index]
                    n = float(fields[indices["n"]])  # type: ignore[index]
                    maf = float(fields[indices["maf"]])  # type: ignore[index]
                except ValueError:
                    continue
                if maf > 0.5:
                    maf = 1.0 - maf
                if not (math.isfinite(n) and n > 0 and math.isfinite(maf) and maf >= 0.01):
                    counts["maf_fail"] += 1
                    continue
                ref1, ref2 = reference[snp]
                z = aligned_z(beta, se, a1, a2, ref1, ref2)
                if z is None:
                    counts["allele_mismatch"] += 1
                    continue
                counts["valid_aligned"] += 1
                row = f"{snp}\t{n:g}\t{z:.12g}\t{ref1}\t{ref2}\n"
                handles["reference_maf01"].write(row)

                try:
                    direction = fields[indices["direction"]] if indices["direction"] is not None else ""
                    studies = sum(char in "+-" for char in direction)
                    het_df = float(fields[indices["het_df"]]) if indices["het_df"] is not None else studies - 1
                    het_p = float(fields[indices["het_p"]]) if indices["het_p"] is not None else 1.0
                    min_freq = float(fields[indices["min_freq"]]) if indices["min_freq"] is not None else maf
                    max_freq = float(fields[indices["max_freq"]]) if indices["max_freq"] is not None else maf
                    af_diff = max_freq - min_freq
                except ValueError:
                    continue
                common_qc = studies >= 2 and het_df >= 1 and het_p > 0.05
                if common_qc and af_diff < 0.50:
                    counts["archived_qc"] += 1
                    handles["archived_qc_afdiff050"].write(row)
                if common_qc and af_diff < 0.05:
                    counts["strict_qc"] += 1
                    handles["strict_qc_afdiff005"].write(row)
    finally:
        for stream in handles.values():
            stream.close()

    return {
        "label": label,
        "source": str(path),
        "sha256": sha256(path),
        "duplicates_removed": duplicate_count,
        **counts,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-alleles", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input", action="append", required=True, metavar="LABEL=PATH")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference = load_reference(args.reference_alleles)
    records = []
    for item in args.input:
        label, raw_path = item.split("=", 1)
        records.append(prepare(label, Path(raw_path), reference, args.output_dir))
    with (args.output_dir / "preparation_manifest.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(f"Prepared {len(records)} freezes against {len(reference):,} reference variants")


if __name__ == "__main__":
    main()
