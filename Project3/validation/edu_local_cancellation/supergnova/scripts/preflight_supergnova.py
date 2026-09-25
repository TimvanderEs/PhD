#!/usr/bin/env python3
"""Run SUPERGNOVA's own preparation routine and record the retained analysis set."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--supergnova-dir", required=True, type=Path)
    parser.add_argument("--bfile", required=True)
    parser.add_argument("--partition", required=True)
    parser.add_argument("--sumstats1", required=True)
    parser.add_argument("--sumstats2", required=True)
    parser.add_argument("--n1", required=True, type=float)
    parser.add_argument("--n2", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.supergnova_dir))
    from prep import prep  # pylint: disable=import-error,import-outside-toplevel

    data, regions, n1, n2 = prep(
        args.bfile,
        args.partition,
        args.sumstats1,
        args.sumstats2,
        args.n1,
        args.n2,
    )
    chromosome_counts = data.groupby("CHR").size().astype(int).to_dict()
    report = {
        "retained_variants": int(len(data)),
        "regions": int(len(regions)),
        "n1": float(n1),
        "n2": float(n2),
        "mean_z1_squared": float(np.mean(np.square(data["Z_x"]))),
        "mean_z2_squared": float(np.mean(np.square(data["Z_y"]))),
        "z_product_mean": float(np.mean(data["Z_x"] * data["Z_y"])),
        "chromosome_variant_counts": {str(int(key)): value for key, value in chromosome_counts.items()},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
