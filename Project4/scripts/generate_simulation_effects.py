#!/usr/bin/env python3
"""Select common SNPs and assign reproducible Gaussian effects."""

from __future__ import annotations

import argparse
import random
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", type=Path, required=True)
    parser.add_argument("--extract", type=Path, required=True)
    parser.add_argument("--causal", type=int, default=10000)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    allowed = {line.split()[0] for line in args.extract.open() if line.split()}
    variants: list[tuple[str, str]] = []
    with args.bim.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 6 and fields[1] in allowed:
                variants.append((fields[1], fields[4]))
    if len(variants) < args.causal:
        raise ValueError(f"Only {len(variants)} eligible SNPs for {args.causal} causal effects")
    rng = random.Random(args.seed)
    selected = rng.sample(variants, args.causal)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write("ID\tA1\tBETA\n")
        for marker, allele in selected:
            handle.write(f"{marker}\t{allele}\t{rng.gauss(0.0, 1.0):.12g}\n")


if __name__ == "__main__":
    main()
