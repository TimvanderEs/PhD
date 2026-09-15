#!/usr/bin/env python3
"""Exercise exact LAVA eligibility and within-pair FDR consolidation."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary_directory:
        out_dir = Path(temporary_directory)
        command = [
            sys.executable,
            str(ROOT / "scripts/lava/postprocess_lava.py"),
            "--univ", str(ROOT / "examples/synthetic/lava/univ.tsv"),
            "--bivar", str(ROOT / "examples/synthetic/lava/bivar.tsv"),
            "--ancestry", "SYNTHETIC",
            "--n-regions", "100",
            "--trait-aliases", str(ROOT / "scripts/lava/config/EUR_trait_aliases.tsv"),
            "--out-dir", str(out_dir),
            "--expected-total", "3",
        ]
        subprocess.run(command, check=True)

        details = read_tsv(out_dir / "lava_bivariate_eligible.tsv")
        assert len(details) == 3
        by_key = {(row["trait_pair"], row["lava_locus_id"]): row for row in details}
        assert float(by_key[("EA--CF", "1")]["univ_p_trait1"]) == 0.0
        assert float(by_key[("EA--CF", "1")]["rho_q_FDR_within_pair"]) == 0.02
        assert float(by_key[("EA--CF", "3")]["rho_q_FDR_within_pair"]) == 0.04
        assert float(by_key[("MDD--SCZ", "2")]["rho_q_FDR_within_pair"]) == 0.02

        qc = {row["metric"]: row["value"] for row in read_tsv(out_dir / "lava_postprocess_qc.tsv")}
        assert qc["univariate_threshold"] == "0.0005"
        assert qc["raw_bivariate_rows"] == "5"
        assert qc["eligible_bivariate_rows"] == "3"
        assert qc["rows_removed_by_exact_univariate_threshold"] == "2"

    print("PASS: LAVA exact-threshold and within-pair FDR post-processing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
