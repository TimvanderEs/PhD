#!/usr/bin/env python3
"""Convert a PLINK2 additive linear-GWAS table to LDSC summary-stat format."""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    opener = gzip.open if args.glm.suffix == ".gz" else open
    written = 0
    with opener(args.glm, "rt") as source, gzip.open(args.output, "wt") as target:
        header = source.readline().lstrip("#").split()
        positions = {name: header.index(name) for name in ("ID", "A1", "REF", "OBS_CT", "BETA", "SE")}
        test_index = header.index("TEST") if "TEST" in header else None
        target.write("SNP\tA1\tA2\tN\tZ\n")
        for line in source:
            fields = line.split()
            if test_index is not None and fields[test_index] != "ADD":
                continue
            try:
                beta = float(fields[positions["BETA"]])
                se = float(fields[positions["SE"]])
                z = beta / se
            except (ValueError, ZeroDivisionError):
                continue
            if not math.isfinite(z):
                continue
            target.write(
                f"{fields[positions['ID']]}\t{fields[positions['A1']]}\t"
                f"{fields[positions['REF']]}\t{fields[positions['OBS_CT']]}\t{z:.12g}\n"
            )
            written += 1
    if written == 0:
        raise ValueError(f"No additive variants written from {args.glm}")


if __name__ == "__main__":
    main()
