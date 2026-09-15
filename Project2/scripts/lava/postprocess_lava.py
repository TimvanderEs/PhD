#!/usr/bin/env python3
"""Apply the manuscript LAVA eligibility rule and within-pair BH correction."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path


PAIR_ORDER = (
    ("EA", "CF"),
    ("MDD", "EA"),
    ("MDD", "CF"),
    ("MDD", "SCZ"),
    ("EA", "SCZ"),
    ("SCZ", "CF"),
)
PAIR_BY_SET = {frozenset(pair): pair for pair in PAIR_ORDER}
PAIR_RANK = {"--".join(pair): index for index, pair in enumerate(PAIR_ORDER)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter LAVA bivariate output to loci where both traits meet "
            "alpha/n-regions, then apply Benjamini-Hochberg correction within "
            "each trait pair. Numerically zero P values are retained."
        )
    )
    parser.add_argument("--univ", required=True, type=Path)
    parser.add_argument("--bivar", required=True, type=Path)
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--n-regions", required=True, type=int)
    parser.add_argument("--trait-aliases", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument(
        "--expected-total",
        type=int,
        help="Optional assertion for the final eligible bivariate row count.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return list(reader.fieldnames), list(reader)


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")


def finite_probability(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {label}: {value!r}") from error
    if not math.isfinite(parsed) or parsed < 0 or parsed > 1:
        raise ValueError(f"Invalid {label}: {value!r}")
    return parsed


def load_aliases(path: Path) -> dict[str, str]:
    fields, rows = read_tsv(path)
    require_columns(path, fields, {"source_trait", "publication_trait"})
    aliases: dict[str, str] = {}
    for row in rows:
        source = row["source_trait"].strip()
        publication = row["publication_trait"].strip()
        if not source or not publication:
            raise ValueError(f"Blank trait alias in {path}")
        if source in aliases and aliases[source] != publication:
            raise ValueError(f"Conflicting alias for {source!r}")
        aliases[source] = publication
    return aliases


def canonical_pair(first: str, second: str) -> tuple[str, str]:
    if first == second:
        raise ValueError(f"A bivariate row repeats trait {first!r}")
    return PAIR_BY_SET.get(frozenset((first, second)), tuple(sorted((first, second))))


def bh_adjust(rows: list[dict[str, object]]) -> list[float]:
    count = len(rows)
    order = sorted(range(count), key=lambda index: float(rows[index]["rho_p"]))
    adjusted = [1.0] * count
    running = 1.0
    for position in range(count - 1, -1, -1):
        index = order[position]
        rank = position + 1
        running = min(running, float(rows[index]["rho_p"]) * count / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def format_number(value: object) -> str:
    if isinstance(value, float):
        return format(value, ".15g")
    return str(value)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_number(row.get(field, "")) for field in fields})


def main() -> int:
    args = parse_args()
    if args.n_regions <= 0:
        raise ValueError("--n-regions must be positive")
    if not 0 < args.alpha <= 1:
        raise ValueError("--alpha must be in (0, 1]")

    aliases = load_aliases(args.trait_aliases)
    univ_fields, univ_rows = read_tsv(args.univ)
    bivar_fields, bivar_rows = read_tsv(args.bivar)
    require_columns(args.univ, univ_fields, {"locus", "phen", "p"})
    require_columns(
        args.bivar,
        bivar_fields,
        {"locus", "chr", "start", "stop", "n.snps", "phen1", "phen2", "rho", "p"},
    )

    threshold = args.alpha / args.n_regions
    univ_p: dict[tuple[str, str], float] = {}
    traits_by_locus: dict[str, list[str]] = defaultdict(list)
    for row in univ_rows:
        key = (row["locus"], row["phen"])
        if key in univ_p:
            raise ValueError(f"Duplicate univariate locus/trait row: {key}")
        if row["phen"] not in aliases:
            raise ValueError(f"No publication alias for univariate trait {row['phen']!r}")
        univ_p[key] = finite_probability(row["p"], "univariate P value")
        traits_by_locus[row["locus"]].append(row["phen"])

    observed_keys: set[tuple[str, frozenset[str]]] = set()
    eligible_rows: list[dict[str, object]] = []
    loose_only = 0
    for row in bivar_rows:
        source_pair = frozenset((row["phen1"], row["phen2"]))
        key = (row["locus"], source_pair)
        if key in observed_keys:
            raise ValueError(f"Duplicate bivariate locus/pair row: {key}")
        observed_keys.add(key)
        for trait in source_pair:
            if trait not in aliases:
                raise ValueError(f"No publication alias for bivariate trait {trait!r}")
            if (row["locus"], trait) not in univ_p:
                raise ValueError(f"Missing univariate result for locus {row['locus']} and {trait}")
        pair_p = [univ_p[(row["locus"], trait)] for trait in source_pair]
        if not all(value <= threshold for value in pair_p):
            loose_only += 1
            continue
        publication_p = {
            aliases[row["phen1"]]: univ_p[(row["locus"], row["phen1"])],
            aliases[row["phen2"]]: univ_p[(row["locus"], row["phen2"])],
        }
        first, second = canonical_pair(aliases[row["phen1"]], aliases[row["phen2"]])
        eligible_rows.append(
            {
                "ancestry": args.ancestry,
                "trait_pair": f"{first}--{second}",
                "trait1": first,
                "trait2": second,
                "lava_locus_id": row["locus"],
                "chr": row["chr"],
                "start_bp": row["start"],
                "end_bp": row["stop"],
                "rho": float(row["rho"]),
                "rho_p": finite_probability(row["p"], "bivariate P value"),
                "univ_p_trait1": publication_p[first],
                "univ_p_trait2": publication_p[second],
                "n_snps": row["n.snps"],
            }
        )

    expected_keys: set[tuple[str, frozenset[str]]] = set()
    for locus, traits in traits_by_locus.items():
        passing = sorted(trait for trait in traits if univ_p[(locus, trait)] <= threshold)
        expected_keys.update((locus, frozenset(pair)) for pair in itertools.combinations(passing, 2))
    missing_keys = expected_keys.difference(observed_keys)
    if missing_keys:
        examples = "; ".join(
            f"{locus}:{','.join(sorted(pair))}" for locus, pair in sorted(missing_keys)[:5]
        )
        raise ValueError(
            f"Bivariate file is missing {len(missing_keys)} pairs eligible at the stated threshold; "
            f"rerun LAVA before post-processing. Examples: {examples}"
        )

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in eligible_rows:
        grouped[str(row["trait_pair"])].append(row)
    for rows in grouped.values():
        for row, adjusted in zip(rows, bh_adjust(rows)):
            row["rho_q_FDR_within_pair"] = adjusted
            row["significant_FDR"] = adjusted < args.alpha
            row["direction"] = "positive" if float(row["rho"]) > 0 else "negative" if float(row["rho"]) < 0 else "zero"

    eligible_rows.sort(
        key=lambda row: (
            PAIR_RANK.get(str(row["trait_pair"]), len(PAIR_RANK)),
            int(str(row["lava_locus_id"])),
        )
    )
    if args.expected_total is not None and len(eligible_rows) != args.expected_total:
        raise ValueError(
            f"Expected {args.expected_total} eligible rows but obtained {len(eligible_rows)}"
        )

    summary_rows: list[dict[str, object]] = []
    for pair in sorted(grouped, key=lambda value: (PAIR_RANK.get(value, len(PAIR_RANK)), value)):
        rows = grouped[pair]
        significant = [row for row in rows if bool(row["significant_FDR"])]
        summary_rows.append(
            {
                "ancestry": args.ancestry,
                "trait_pair": pair,
                "n_bivariate_tests": len(rows),
                "n_FDR_significant": len(significant),
                "n_FDR_positive": sum(row["direction"] == "positive" for row in significant),
                "n_FDR_negative": sum(row["direction"] == "negative" for row in significant),
                "FDR_family": "ancestry_by_trait_pair",
            }
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_fields = [
        "ancestry", "trait_pair", "trait1", "trait2", "lava_locus_id", "chr",
        "start_bp", "end_bp", "rho", "rho_p", "rho_q_FDR_within_pair",
        "significant_FDR", "direction", "univ_p_trait1", "univ_p_trait2", "n_snps",
    ]
    write_tsv(args.out_dir / "lava_bivariate_eligible.tsv", detail_fields, eligible_rows)
    write_tsv(
        args.out_dir / "lava_pair_summary.tsv",
        [
            "ancestry", "trait_pair", "n_bivariate_tests", "n_FDR_significant",
            "n_FDR_positive", "n_FDR_negative", "FDR_family",
        ],
        summary_rows,
    )
    qc_rows = [
        {"metric": "ancestry", "value": args.ancestry},
        {"metric": "n_regions", "value": args.n_regions},
        {"metric": "univariate_alpha", "value": args.alpha},
        {"metric": "univariate_threshold", "value": threshold},
        {"metric": "raw_bivariate_rows", "value": len(bivar_rows)},
        {"metric": "eligible_bivariate_rows", "value": len(eligible_rows)},
        {"metric": "rows_removed_by_exact_univariate_threshold", "value": loose_only},
        {"metric": "missing_eligible_bivariate_rows", "value": len(missing_keys)},
        {"metric": "FDR_family", "value": "ancestry_by_trait_pair"},
    ]
    write_tsv(args.out_dir / "lava_postprocess_qc.tsv", ["metric", "value"], qc_rows)
    print(
        f"PASS: {args.ancestry} retained {len(eligible_rows)} eligible rows; "
        f"BH correction applied within {len(grouped)} trait pairs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
