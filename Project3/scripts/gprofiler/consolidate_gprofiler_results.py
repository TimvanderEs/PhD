#!/usr/bin/env python3
"""Validate and consolidate the archived publication g:Profiler CSV exports.

The script never calls the live g:Profiler service. This is intentional: the
original web exports did not record a database release, and live annotation
databases change. The archived query/background lists and CSVs are therefore
the reproducible record of the analysis as run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path


SCRIPT_VERSION = "1.0.0"
REQUIRED_RAW_COLUMNS = {
    "source",
    "term_name",
    "term_id",
    "adjusted_p_value",
    "negative_log10_of_adjusted_p_value",
    "term_size",
    "query_size",
    "intersection_size",
    "effective_domain_size",
    "intersections",
}
SOURCE_ORDER = [
    "GO:BP",
    "GO:MF",
    "GO:CC",
    "REAC",
    "KEGG",
    "WP",
    "HP",
    "HPA",
    "CORUM",
    "TF",
    "MIRNA",
]
MODEL_ORDER = {
    "FOURTRAIT": 0,
    "MDD_3TRAIT": 1,
    "SCZ_3TRAIT": 2,
    "EAS_4TRAIT": 0,
    "EAS_MDD_EA_CF": 1,
    "EAS_SCZ_EA_CF": 2,
}
DIRECTION_ORDER = {
    "concordant": 0,
    "discordant": 1,
    "dual": 2,
    "all_mapped_genes": 3,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=SCRIPT_VERSION)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path.cwd(),
        help="Root used to resolve paths stored in the run manifest.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("inputs/gprofiler/run_manifest.tsv"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/gprofiler/processed"),
    )
    return parser.parse_args()


def sha256_file(path: Path, decompress: bool = False) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if decompress else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gene_list_details(path: Path) -> tuple[int, int, str, set[str]]:
    genes = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines()]
    genes = [gene for gene in genes if gene]
    return len(genes), len(set(genes)), sha256_file(path), set(genes)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Manifest contains no runs: {path}")
    required = {
        "run_id",
        "ancestry",
        "model",
        "direction_class",
        "query_file",
        "query_genes",
        "query_sha256",
        "background_scope",
        "background_file",
        "background_genes",
        "background_sha256",
        "raw_file",
        "raw_rows",
        "significant_rows",
        "original_csv_sha256",
        "compressed_sha256",
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(sorted(missing))}")
    run_ids = [row["run_id"] for row in rows]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("Manifest run_id values are not unique")
    return rows


def broad_category(source: str) -> str:
    return {
        "GO:BP": "GO_Biological_Process",
        "GO:MF": "GO_Molecular_Function",
        "GO:CC": "GO_Cellular_Component",
        "HP": "Human_Phenotype",
        "WP": "Pathway",
        "REAC": "Pathway",
        "KEGG": "Pathway",
    }.get(source, source or "Unknown")


def source_label(source: str) -> str:
    return {
        "GO:BP": "GO biological process",
        "GO:MF": "GO molecular function",
        "GO:CC": "GO cellular component",
        "REAC": "Reactome",
        "KEGG": "KEGG pathway",
        "WP": "WikiPathways",
        "HP": "Human phenotype",
        "HPA": "Human Protein Atlas",
        "CORUM": "CORUM",
        "TF": "Transcription factor target",
        "MIRNA": "miRNA target",
    }.get(source, source)


def as_int(value: str) -> int:
    return int(float(value))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    manifest_path = args.manifest
    if not manifest_path.is_absolute():
        manifest_path = root / manifest_path
    outdir = args.outdir
    if not outdir.is_absolute():
        outdir = root / outdir

    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(2**31 - 1)

    manifest = read_manifest(manifest_path)
    validation_rows: list[dict] = []
    run_summaries: list[dict] = []
    source_summaries: list[dict] = []
    significant_terms: list[dict] = []

    for meta in manifest:
        raw_path = root / meta["raw_file"]
        query_path = root / meta["query_file"]
        if not raw_path.is_file() or not query_path.is_file():
            raise FileNotFoundError(
                f"Missing raw/query file for {meta['run_id']}: {raw_path} / {query_path}"
            )

        query_n, query_unique, query_sha, query_set = gene_list_details(query_path)
        background_n = 0
        background_unique = 0
        background_sha = ""
        background_contains_query = "NA"
        if meta["background_file"]:
            background_path = root / meta["background_file"]
            if not background_path.is_file():
                raise FileNotFoundError(background_path)
            background_n, background_unique, background_sha, background_set = gene_list_details(
                background_path
            )
            background_contains_query = str(query_set.issubset(background_set))

        original_sha = sha256_file(raw_path, decompress=raw_path.suffix == ".gz")
        compressed_sha = sha256_file(raw_path)
        checks = {
            "query_count_match": query_n == int(meta["query_genes"]),
            "query_unique": query_n == query_unique,
            "query_sha256_match": query_sha == meta["query_sha256"],
            "background_count_match": (
                not meta["background_file"]
                or background_n == int(meta["background_genes"])
            ),
            "background_unique": (
                not meta["background_file"] or background_n == background_unique
            ),
            "background_sha256_match": (
                not meta["background_file"]
                or background_sha == meta["background_sha256"]
            ),
            "original_csv_sha256_match": original_sha == meta["original_csv_sha256"],
            "compressed_sha256_match": compressed_sha == meta["compressed_sha256"],
        }
        failed = [name for name, value in checks.items() if not value]
        if failed:
            raise ValueError(f"{meta['run_id']} failed: {', '.join(failed)}")

        opener = gzip.open if raw_path.suffix == ".gz" else open
        n_rows = 0
        n_sig = 0
        min_p = math.inf
        source_stats: dict[str, dict[str, float | int]] = defaultdict(
            lambda: {"n_returned_terms": 0, "n_significant_terms": 0, "min_adjusted_p": math.inf}
        )
        query_sizes: set[int] = set()
        domain_sizes: set[int] = set()
        with opener(raw_path, "rt", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            missing = REQUIRED_RAW_COLUMNS.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"{raw_path} is missing columns: {', '.join(sorted(missing))}"
                )
            for raw in reader:
                n_rows += 1
                p_value = float(raw["adjusted_p_value"])
                min_p = min(min_p, p_value)
                query_sizes.add(as_int(raw["query_size"]))
                domain_sizes.add(as_int(raw["effective_domain_size"]))
                source = raw["source"]
                stats = source_stats[source]
                stats["n_returned_terms"] = int(stats["n_returned_terms"]) + 1
                stats["min_adjusted_p"] = min(float(stats["min_adjusted_p"]), p_value)
                if p_value < float(meta["user_threshold"]):
                    n_sig += 1
                    stats["n_significant_terms"] = int(stats["n_significant_terms"]) + 1
                    significant_terms.append(
                        {
                            "ancestry": meta["ancestry"],
                            "model": meta["model"],
                            "direction_class": meta["direction_class"],
                            "source": source,
                            "source_label": source_label(source),
                            "term_id": raw["term_id"],
                            "term_name": raw["term_name"],
                            "adjusted_p_value": raw["adjusted_p_value"],
                            "negative_log10_adjP": -math.log10(p_value),
                            "term_size": raw["term_size"],
                            "query_size": raw["query_size"],
                            "intersection_size": raw["intersection_size"],
                            "effective_domain_size": raw["effective_domain_size"],
                            "intersection_genes_compact": raw["intersections"].replace(
                                ",", ";"
                            ),
                            "broad_category": broad_category(source),
                            "source_file": raw_path.name.removesuffix(".gz"),
                            "_p": p_value,
                        }
                    )

        row_match = n_rows == int(meta["raw_rows"])
        significant_match = n_sig == int(meta["significant_rows"])
        if not row_match or not significant_match:
            raise ValueError(
                f"{meta['run_id']} count mismatch: rows={n_rows}, significant={n_sig}"
            )

        validation_rows.append(
            {
                "run_id": meta["run_id"],
                **{name: str(value) for name, value in checks.items()},
                "raw_row_count_match": str(row_match),
                "significant_row_count_match": str(significant_match),
                "background_contains_query": background_contains_query,
                "status": "PASS",
            }
        )
        run_summaries.append(
            {
                "run_id": meta["run_id"],
                "ancestry": meta["ancestry"],
                "model": meta["model"],
                "direction_class": meta["direction_class"],
                "input_query_genes": query_n,
                "input_background_genes": background_n if meta["background_file"] else "NA",
                "background_scope": meta["background_scope"],
                "n_returned_terms": n_rows,
                "n_significant_terms": n_sig,
                "min_adjusted_p": min_p,
                "query_sizes_across_sources": ";".join(map(str, sorted(query_sizes))),
                "effective_domain_sizes_across_sources": ";".join(
                    map(str, sorted(domain_sizes))
                ),
                "correction_method": meta["correction_method"],
                "user_threshold": meta["user_threshold"],
                "database_release": meta["database_release"],
                "raw_file": meta["raw_file"],
                "original_csv_sha256": original_sha,
            }
        )
        for source in SOURCE_ORDER:
            stats = source_stats.get(source)
            if stats is None:
                continue
            source_summaries.append(
                {
                    "run_id": meta["run_id"],
                    "ancestry": meta["ancestry"],
                    "model": meta["model"],
                    "direction_class": meta["direction_class"],
                            "source": source,
                            "source_label": source_label(source),
                    "n_returned_terms": stats["n_returned_terms"],
                    "n_significant_terms": stats["n_significant_terms"],
                    "min_adjusted_p": stats["min_adjusted_p"],
                }
            )

    significant_terms.sort(
        key=lambda row: (
            0 if row["ancestry"] == "EAS" else 1,
            MODEL_ORDER.get(row["model"], 99),
            DIRECTION_ORDER.get(row["direction_class"], 99),
            row["_p"],
            row["source"],
            row["term_id"],
        )
    )
    significant_fields = [
        "ancestry",
        "model",
        "direction_class",
        "source",
        "source_label",
        "term_id",
        "term_name",
        "adjusted_p_value",
        "negative_log10_adjP",
        "term_size",
        "query_size",
        "intersection_size",
        "effective_domain_size",
        "intersection_genes_compact",
        "broad_category",
        "source_file",
    ]
    write_tsv(outdir / "significant_terms.tsv", significant_terms, significant_fields)
    write_tsv(
        outdir / "run_summary.tsv",
        run_summaries,
        list(run_summaries[0]),
    )
    write_tsv(
        outdir / "source_summary.tsv",
        source_summaries,
        list(source_summaries[0]),
    )
    write_tsv(
        outdir / "manifest_validation.tsv",
        validation_rows,
        list(validation_rows[0]),
    )
    print(f"Validated {len(manifest)} archived g:Profiler runs")
    print(f"Significant terms (g:SCS-adjusted P < 0.05): {len(significant_terms)}")
    print(f"Wrote: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
