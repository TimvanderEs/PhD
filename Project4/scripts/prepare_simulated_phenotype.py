#!/usr/bin/env python3
"""Create a covariate-orthogonal quantitative phenotype with exact realized h2."""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path


def center(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    return [value - mean for value in values]


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def unit(values: list[float]) -> list[float]:
    values = center(values)
    norm = math.sqrt(dot(values, values))
    if norm == 0:
        raise ValueError("Cannot normalize a zero-variance vector")
    return [value / norm for value in values]


def residualize(values: list[float], basis: list[list[float]]) -> list[float]:
    output = center(values)
    for vector in basis:
        coefficient = dot(output, vector)
        output = [value - coefficient * component for value, component in zip(output, vector)]
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--eigenvec", type=Path, required=True)
    parser.add_argument("--h2", type=float, default=0.2)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.h2 < 1:
        raise ValueError("h2 must lie strictly between zero and one")

    scores: dict[str, tuple[str, float]] = {}
    with args.score.open() as handle:
        header = handle.readline().split()
        score_column = "SCORE" if "SCORE" in header else header[-1]
        iid_index = header.index("IID")
        fid_index = header.index("FID")
        score_index = header.index(score_column)
        for line in handle:
            fields = line.split()
            scores[fields[iid_index]] = (fields[fid_index], float(fields[score_index]))

    with args.eigenvec.open() as handle:
        header = handle.readline().lstrip("#").split()
        iid_index = header.index("IID")
        pc_indices = [index for index, value in enumerate(header) if value.startswith("PC")]
        ids: list[tuple[str, str]] = []
        pcs: list[list[float]] = []
        genetic: list[float] = []
        for line in handle:
            fields = line.split()
            iid = fields[iid_index]
            if iid not in scores:
                raise ValueError(f"PCA sample missing from score file: {iid}")
            fid, score = scores[iid]
            ids.append((fid, iid))
            pcs.append([float(fields[index]) for index in pc_indices])
            genetic.append(score)
    if len(ids) != len(scores):
        raise ValueError("Score and PCA sample sets differ")

    # Modified Gram-Schmidt gives an orthonormal PC basis even after text rounding.
    basis: list[list[float]] = []
    for column in range(len(pcs[0])):
        vector = residualize([row[column] for row in pcs], basis)
        basis.append(unit(vector))
    g = unit(residualize(genetic, basis))

    rng = random.Random(args.seed)
    noise = [rng.gauss(0.0, 1.0) for _ in ids]
    e = unit(residualize(noise, basis + [g]))
    scale = math.sqrt(len(ids))
    g = [value * scale for value in g]
    e = [value * scale for value in e]
    phenotype = [math.sqrt(args.h2) * gv + math.sqrt(1 - args.h2) * ev for gv, ev in zip(g, e)]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write("FID\tIID\tPHENO\n")
        for (fid, iid), value in zip(ids, phenotype):
            handle.write(f"{fid}\t{iid}\t{value:.12g}\n")

    n = len(ids)
    var_g = sum(value * value for value in g) / n
    var_e = sum(value * value for value in e) / n
    cov_ge = sum(x * y for x, y in zip(g, e)) / n
    var_y = sum(value * value for value in phenotype) / n
    args.report.write_text(
        json.dumps(
            {
                "samples": n,
                "pcs": len(basis),
                "target_h2": args.h2,
                "genetic_variance_before_weight": var_g,
                "environmental_variance_before_weight": var_e,
                "genetic_environmental_covariance": cov_ge,
                "phenotype_variance": var_y,
                "realized_h2": args.h2 * var_g / var_y,
                "seed": args.seed,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
