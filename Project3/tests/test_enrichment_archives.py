#!/usr/bin/env python3
"""Integrity checks for the public FUMA/MAGMA and g:Profiler archives."""

from __future__ import annotations

import configparser
import csv
import hashlib
import math
import subprocess
import sys
import tempfile
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def data_rows(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#")) - 1


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_fuma_checksums(root: Path) -> None:
    checksum_file = root / "provenance/fuma_magma_raw_files.sha256"
    rows = [line.rstrip("\n").split("  ", 1) for line in checksum_file.open() if line.strip()]
    assert len(rows) == 52
    for expected, stored_path in rows:
        path = root.parent / stored_path
        assert path.is_file(), path
        assert sha256_file(path) == expected, path


def check_gprofiler(root: Path) -> None:
    committed = root / "results/gprofiler/processed"
    with tempfile.TemporaryDirectory(prefix="project3_gprofiler_") as tmp:
        generated = Path(tmp)
        subprocess.run(
            [
                sys.executable,
                str(root / "scripts/gprofiler/consolidate_gprofiler_results.py"),
                "--project-root",
                str(root),
                "--outdir",
                str(generated),
            ],
            check=True,
        )
        significant = read_tsv(generated / "significant_terms.tsv")
        assert len(significant) == 61, len(significant)
        assert len(read_tsv(generated / "run_summary.tsv")) == 12
        assert {row["status"] for row in read_tsv(generated / "manifest_validation.tsv")} == {
            "PASS"
        }
        for filename in (
            "significant_terms.tsv",
            "run_summary.tsv",
            "source_summary.tsv",
            "manifest_validation.tsv",
        ):
            assert (generated / filename).read_bytes() == (committed / filename).read_bytes(), (
                f"Committed g:Profiler output is stale: {filename}"
            )


def check_fuma_manifest(root: Path, manifest_name: str, expected_runs: int) -> None:
    manifest = read_tsv(root / f"inputs/fuma_magma/{manifest_name}_run_manifest.tsv")
    assert len(manifest) == expected_runs
    summary_rows = read_tsv(
        root / f"results/fuma_magma/processed/{manifest_name}/01_run_level_summary.tsv"
    )
    expected_magma_rows = {
        (row["ancestry"], row["analysis_label"]): int(row["n_rows_total"])
        for row in summary_rows
        if row["result_type"] == "competitive_gene_set"
    }
    seen: set[tuple[str, str]] = set()
    for row in manifest:
        key = (row["ancestry"], row["analysis_label"])
        assert key not in seen, key
        seen.add(key)
        run_dir = root / row["run_dir"]
        required = (
            "params.config",
            "magma.gsa.out",
            "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out",
            "magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out",
        )
        assert all((run_dir / filename).is_file() for filename in required), run_dir
        assert data_rows(run_dir / "magma.gsa.out") == expected_magma_rows[key]

        config = configparser.ConfigParser(interpolation=None)
        config.read(run_dir / "params.config")
        assert config["params"]["pop"] == row["ancestry"]
        assert config["version"]["MAGMA"] == "v1.08"
        expected_fuma = "v1.8.2" if manifest_name == "pleio" else "v1.5.2"
        assert config["version"]["FUMA"] == expected_fuma
        assert config["params"]["refpanel"] == "1KG/Phase3"
        assert config["params"]["exMHC"] == "1"
        assert config["params"]["ensembl"] == "v102"


def compare_significant_tables(committed: Path, rebuilt: Path) -> None:
    committed_rows = read_tsv(committed)
    rebuilt_rows = read_tsv(rebuilt)
    key_fields = ("ancestry", "analysis_label", "VARIABLE")
    key = lambda row: tuple(row[field] for field in key_fields)
    left = {key(row): row for row in committed_rows}
    right = {key(row): row for row in rebuilt_rows}
    assert set(left) == set(right), committed.name
    for record_key in left:
        for field in ("p_numeric", "beta_numeric", "bonferroni_threshold", "p_bh"):
            assert math.isclose(
                float(left[record_key][field]),
                float(right[record_key][field]),
                rel_tol=1e-12,
                abs_tol=1e-15,
            ), (committed.name, record_key, field)


def check_fuma_rebuild(root: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="project3_fuma_") as tmp:
        tmp_path = Path(tmp)
        single = tmp_path / "single"
        pleio = tmp_path / "pleio"
        comparison = tmp_path / "comparison"
        extractor = root / "scripts/fuma_magma/extract_fuma_magma_gtex.py"
        comparer = root / "scripts/fuma_magma/compare_magma_single_vs_pleio.py"
        for scope, manifest, outdir in (
            ("single-trait", "single_trait_run_manifest.tsv", single),
            ("multivariate", "pleio_run_manifest.tsv", pleio),
        ):
            subprocess.run(
                [
                    sys.executable,
                    str(extractor),
                    "--run-manifest",
                    str(root / "inputs/fuma_magma" / manifest),
                    "--scope",
                    scope,
                    "--outdir",
                    str(outdir),
                ],
                cwd=root,
                check=True,
            )
        subprocess.run(
            [
                sys.executable,
                str(comparer),
                "--single-dir",
                str(single),
                "--pleio-dir",
                str(pleio),
                "--outdir",
                str(comparison),
            ],
            check=True,
        )

        for group, rebuilt_dir in (("single_trait", single), ("pleio", pleio)):
            committed_dir = root / "results/fuma_magma/processed" / group
            for filename in (
                "03_MAGMA_competitive_gene_sets_bonferroni_positive.tsv",
                "07_GTEx_v8_specific_tissues_bonferroni_positive.tsv",
                "10_GTEx_v8_general_tissues_bonferroni_positive.tsv",
            ):
                compare_significant_tables(
                    committed_dir / filename,
                    rebuilt_dir / filename,
                )

        committed_comparison = root / "results/fuma_magma/processed/comparison"
        for filename, key_fields in (
            ("02_model_level_summary.tsv", ("ancestry", "pleio_model", "result_type")),
            (
                "11_missing_or_incomplete_comparisons.tsv",
                ("ancestry", "pleio_model", "result_type"),
            ),
        ):
            committed_rows = read_tsv(committed_comparison / filename)
            rebuilt_rows = read_tsv(comparison / filename)
            stable = lambda row: tuple(row[field] for field in key_fields)
            assert {stable(row): row for row in committed_rows} == {
                stable(row): row for row in rebuilt_rows
            }, filename
        for filename in (
            "04_PLEIO_only_among_available_constituents_bonferroni.tsv",
            "05_PLEIO_shared_with_constituent_bonferroni.tsv",
            "06_constituent_only_bonferroni.tsv",
        ):
            fields = ("ancestry", "pleio_model", "result_type", "term")
            stable_set = lambda path: {
                tuple(row[field] for field in fields) for row in read_tsv(path)
            }
            assert stable_set(committed_comparison / filename) == stable_set(
                comparison / filename
            ), filename


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parents[1]).resolve()
    check_gprofiler(root)
    check_fuma_checksums(root)
    check_fuma_manifest(root, "pleio", 6)
    check_fuma_manifest(root, "single_trait", 7)
    check_fuma_rebuild(root)

    availability = read_tsv(
        root / "results/fuma_magma/processed/single_trait/12_expected_run_audit.tsv"
    )
    eas_cf = [
        row for row in availability
        if row["ancestry"] == "EAS" and row["analysis_label"] == "CF"
    ]
    assert len(eas_cf) == 1 and eas_cf[0]["observed_magma_gsa"] == "False"
    print("PASS: FUMA/MAGMA and g:Profiler archive integrity")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
