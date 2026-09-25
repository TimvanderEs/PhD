#!/usr/bin/env python3
"""Interpolate GRCh37 centimorgan positions into a chromosome-specific BIM."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", required=True, type=Path)
    parser.add_argument("--map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bim = pd.read_csv(
        args.bim,
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        dtype={"CHR": "int64", "SNP": "string", "BP": "int64", "A1": "string", "A2": "string"},
    )
    genetic_map = pd.read_csv(
        args.map,
        sep=r"\s+",
        header=None,
        names=["CHR", "MARKER", "CM", "BP"],
        usecols=[0, 2, 3],
        dtype={"CHR": "int64", "CM": "float64", "BP": "int64"},
    ).drop_duplicates("BP", keep="last")
    genetic_map = genetic_map.sort_values("BP")

    chromosomes = bim["CHR"].unique()
    if len(chromosomes) != 1:
        raise ValueError(f"Expected one chromosome in BIM, found {chromosomes.tolist()}")
    chromosome = int(chromosomes[0])
    if not (genetic_map["CHR"] == chromosome).all():
        raise ValueError("BIM and genetic map chromosome labels do not match")
    if not genetic_map["BP"].is_monotonic_increasing:
        raise ValueError("Genetic-map base-pair positions are not sorted")
    if not genetic_map["CM"].is_monotonic_increasing:
        raise ValueError("Genetic-map centimorgan positions are not monotonic")

    below = int((bim["BP"] < genetic_map["BP"].iloc[0]).sum())
    above = int((bim["BP"] > genetic_map["BP"].iloc[-1]).sum())
    bim["CM"] = np.interp(
        bim["BP"].to_numpy(),
        genetic_map["BP"].to_numpy(),
        genetic_map["CM"].to_numpy(),
    )
    if not bim["CM"].is_monotonic_increasing:
        raise ValueError("Interpolated BIM centimorgan positions are not monotonic")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    bim.to_csv(args.output, sep="\t", header=False, index=False, float_format="%.8f")
    report = {
        "chromosome": chromosome,
        "variants": len(bim),
        "map_markers": len(genetic_map),
        "variants_below_map": below,
        "variants_above_map": above,
        "minimum_cm": float(bim["CM"].min()),
        "maximum_cm": float(bim["CM"].max()),
        "nonzero_cm": int((bim["CM"] != 0).sum()),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
