#!/usr/bin/env python3
"""
Audit and extract FUMA/MAGMA competitive gene-set and GTEx v8 gene-property outputs.

Expected tree:
  <BASE>/{EUR,EAS}/<RUN>/

Primary files:
  magma.gsa.out
  magma.gsa.sets.genes.out
  magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out
  magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out

Key safeguards
--------------
* Raw ancestry, directory labels, file paths and MAGMA columns are preserved.
* Canonical model labels are stored in separate audit columns.
* Ambiguous generic directories such as EAS/PLEIO are NOT silently assigned.
* An optional run-map TSV can override labels explicitly.
* Bonferroni and BH-FDR corrections are calculated within each output file.
* Competitive and GTEx significance tables require a positive MAGMA beta.
* Output-directory writability is checked before extraction.
* Only direct FUMA SNP2GENE run directories under BASE/EUR and BASE/EAS are analysed.
* Single-trait and PLEIO scopes both use exact ancestry/path allowlists.
* An explicit run manifest can bypass discovery and lock extraction to reviewed directories.
* Gene2Func, archive, backup and other post-FUMA directories are excluded.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


FILE_TYPES = {
    "magma.gsa.out": "competitive_gene_set",
    "magma.gsa.sets.genes.out": "gene_set_membership",
    "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out": "gtex_v8_specific_tissue",
    "magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out": "gtex_v8_general_tissue",
}

EXPECTED_TYPE = {
    "competitive_gene_set": "SET",
    "gtex_v8_specific_tissue": "COVAR",
    "gtex_v8_general_tissue": "COVAR",
}

CANONICAL_LABELS = {
    "EA",
    "CF",
    "MDD",
    "SCZ",
    "FOURTRAIT",
    "MDD_3TRAIT",
    "SCZ_3TRAIT",
    "UNRESOLVED",
}

SCRIPT_VERSION = "2.6.0"


# Exact FUMA SNP2GENE directories to include for the single-trait audit.
# This deliberately avoids inferring scope from params.config or filenames.
SINGLE_TRAIT_RUN_ALLOWLIST = {
    ("EUR", "EDU"): "EA",
    ("EUR", "G"): "CF",
    ("EUR", "MDD"): "MDD",
    ("EUR", "SCZ"): "SCZ",
    ("EAS", "EDU"): "EA",
    ("EAS", "MDD"): "MDD",
    ("EAS", "SCZ"): "SCZ",
}


# Exact FUMA SNP2GENE directories for the model-wide PLEIO audit.
# EAS/PLEIO remains unresolved until a reviewed --run-map assigns it.
PLEIO_RUN_ALLOWLIST = {
    ("EUR", "PLEIO"): "FOURTRAIT",
    ("EUR", "FUMA_PLEIO_MDD_EA_CF"): "MDD_3TRAIT",
    ("EUR", "FUMA_PLEIO_SCZ_EA_CF"): "SCZ_3TRAIT",
    ("EAS", "FUMA_EAS_PLEIO_4trait"): "FOURTRAIT",
    ("EAS", "PLEIO"): "UNRESOLVED",
}


EXPECTED_RUNS = [
    "EA",
    "CF",
    "MDD",
    "SCZ",
    "FOURTRAIT",
    "MDD_3TRAIT",
    "SCZ_3TRAIT",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {SCRIPT_VERSION}",
    )
    parser.add_argument(
        "--base",
        required=False,
        type=Path,
        default=None,
        help=(
            "Base FUMA directory containing EUR and EAS subdirectories. "
            "Not required when --run-manifest is supplied."
        ),
    )
    parser.add_argument(
        "--run-manifest",
        type=Path,
        default=None,
        help=(
            "Explicit TSV of FUMA SNP2GENE run directories. Required columns: "
            "ancestry, canonical_label_final, run_dir. This bypasses directory "
            "discovery and path allowlists."
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Writable output directory for audit and extracted TSV files.",
    )
    parser.add_argument(
        "--run-map",
        type=Path,
        default=None,
        help=(
            "Optional TSV overriding canonical labels. Required columns: ancestry, "
            "run_relative_path, canonical_label_manual."
        ),
    )
    parser.add_argument(
        "--alpha",
        default=0.05,
        type=float,
        help="Family-wise/FDR alpha.",
    )
    parser.add_argument(
        "--top-n",
        default=20,
        type=int,
        help="Number of lowest-P rows retained per run in top-results tables.",
    )
    parser.add_argument(
        "--scope",
        choices=["all-snp2gene", "single-trait", "multivariate"],
        default="all-snp2gene",
        help=(
            "Restrict extraction to all direct SNP2GENE runs, single-trait runs "
            "(EA/CF/MDD/SCZ), or multivariate PLEIO runs."
        ),
    )
    return parser.parse_args()


def normalise(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")


def contains_token(s: str, tokens: Iterable[str]) -> bool:
    parts = set(s.split("_"))
    return any(token in parts for token in tokens)


def infer_label_from_text(raw_run: str, context_text: str = "") -> tuple[str, str, str]:
    """Infer only when the label is reasonably explicit."""
    s = normalise(raw_run)
    context = normalise(f"{raw_run} {context_text}")

    # Explicit four-trait naming.
    if re.search(r"(^|_)(4|four)_?trait($|_)", context) or "fourtrait" in context:
        return "FOURTRAIT", "inferred_explicit", "explicit 4-trait token"

    # Explicit component naming. These observed directory names do not contain "3trait".
    has_ea = contains_token(context, ["ea", "edu", "education"])
    has_cf = contains_token(context, ["cf", "g", "helios", "cognition", "cognitive"])
    has_mdd = contains_token(context, ["mdd", "depression"])
    has_scz = contains_token(context, ["scz", "schizophrenia"])

    if has_mdd and has_ea and has_cf and not has_scz:
        return "MDD_3TRAIT", "inferred_components", "contains MDD+EA+CF components"
    if has_scz and has_ea and has_cf and not has_mdd:
        return "SCZ_3TRAIT", "inferred_components", "contains SCZ+EA+CF components"

    exact = {
        "ea": "EA",
        "edu": "EA",
        "education": "EA",
        "educational_attainment": "EA",
        "cf": "CF",
        "g": "CF",
        "cognitive_function": "CF",
        "cognition": "CF",
        "helios": "CF",
        "helios_g": "CF",
        "mdd": "MDD",
        "depression": "MDD",
        "scz": "SCZ",
        "schizophrenia": "SCZ",
    }
    if s in exact:
        return exact[s], "inferred_exact", f"exact directory alias: {s}"

    if s == "pleio" or "pleio" in s:
        return "UNRESOLVED", "ambiguous", "generic PLEIO label"

    return "UNRESOLVED", "ambiguous", "no safe label rule"


def read_small_text(path: Path, max_chars: int = 100_000) -> str:
    if not path.exists():
        return ""
    try:
        return path.read_text(errors="replace")[:max_chars]
    except Exception:
        return ""


def compact_text(text: str, max_chars: int = 4000) -> str:
    return re.sub(r"\s+", " ", text).strip()[:max_chars]


def discover_runs(base: Path) -> tuple[list[dict], list[dict]]:
    """
    Discover only direct child directories of BASE/EUR and BASE/EAS that contain
    FUMA SNP2GENE MAGMA outputs. Nested Gene2Func/post-FUMA folders and archive
    trees are deliberately excluded from analysis.
    """
    excluded_name_tokens = {
        "gene2func", "gene_2_func", "archive", "backup", "duplicate",
        "old", "tmp", "temp", "post_fuma", "postfuma"
    }

    all_run_dirs = {
        p.parent
        for filename in FILE_TYPES
        for p in base.glob(f"**/{filename}")
        if p.is_file()
    }

    standard_run_dirs: set[Path] = set()
    skipped: list[dict] = []

    for ancestry in ("EUR", "EAS"):
        ancestry_root = base / ancestry
        if not ancestry_root.exists():
            continue
        for child in sorted(ancestry_root.iterdir()):
            if not child.is_dir():
                continue
            norm_name = normalise(child.name)
            if any(token in norm_name for token in excluded_name_tokens):
                skipped.append({
                    "source_directory": str(child),
                    "relative_path": str(child.relative_to(base)),
                    "reason": "excluded non-SNP2GENE/post-FUMA directory name",
                })
                continue
            if any((child / filename).is_file() for filename in FILE_TYPES):
                standard_run_dirs.add(child)

    for path in sorted(all_run_dirs.difference(standard_run_dirs)):
        rel = str(path.relative_to(base))
        if not any(row["source_directory"] == str(path) for row in skipped):
            skipped.append({
                "source_directory": str(path),
                "relative_path": rel,
                "reason": "not a direct SNP2GENE run directory under BASE/EUR or BASE/EAS",
            })

    runs: list[dict] = []
    for run_dir in sorted(standard_run_dirs):
        rel = run_dir.relative_to(base)
        ancestry = rel.parts[0]
        run_relative = str(run_dir.relative_to(base / ancestry))
        raw_run = run_dir.name
        top_files = sorted(p.name for p in run_dir.iterdir() if p.is_file())
        params_text = read_small_text(run_dir / "params.config")
        context = " ".join(top_files) + " " + params_text
        label, status, reason = infer_label_from_text(raw_run, context)
        runs.append(
            {
                "ancestry": ancestry,
                "run_label_raw": raw_run,
                "run_relative_path": run_relative,
                "run_key": f"{ancestry}/{run_relative}",
                "run_dir": run_dir,
                "canonical_label_auto": label,
                "canonical_label_status_auto": status,
                "canonical_label_reason_auto": reason,
                "canonical_label_final": label,
                "canonical_label_status_final": status,
                "canonical_label_reason_final": reason,
                "top_level_files": "|".join(top_files),
                "params_config_compact": compact_text(params_text),
            }
        )

    # Contextual EUR rule: generic PLEIO plus explicit MDD/SCZ three-trait runs
    # identifies the remaining generic run as FOURTRAIT. Do not apply this to EAS
    # where an explicit four-trait directory already exists.
    by_ancestry: dict[str, list[dict]] = {}
    for run in runs:
        by_ancestry.setdefault(run["ancestry"], []).append(run)

    for ancestry, group in by_ancestry.items():
        known = {
            r["canonical_label_auto"]
            for r in group
            if r["canonical_label_auto"] != "UNRESOLVED"
        }
        generic = [
            r
            for r in group
            if normalise(r["run_label_raw"]) == "pleio"
            and r["canonical_label_auto"] == "UNRESOLVED"
        ]
        if (
            len(generic) == 1
            and "FOURTRAIT" not in known
            and {"MDD_3TRAIT", "SCZ_3TRAIT"}.issubset(known)
        ):
            r = generic[0]
            r["canonical_label_auto"] = "FOURTRAIT"
            r["canonical_label_status_auto"] = "inferred_sibling_context"
            r["canonical_label_reason_auto"] = (
                "generic PLEIO is the remaining model beside explicit MDD_3TRAIT "
                "and SCZ_3TRAIT runs"
            )
            r["canonical_label_final"] = r["canonical_label_auto"]
            r["canonical_label_status_final"] = r["canonical_label_status_auto"]
            r["canonical_label_reason_final"] = r["canonical_label_reason_auto"]

    return runs, skipped



def runs_from_manifest(manifest_path: Path) -> tuple[list[dict], list[dict]]:
    """Construct an explicit audited run list from a user-reviewed TSV manifest."""
    if not manifest_path.exists():
        raise FileNotFoundError(f"Run manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    required = {"ancestry", "canonical_label_final", "run_dir"}
    missing = required.difference(manifest.columns)
    if missing:
        raise ValueError(
            f"Run manifest is missing columns: {', '.join(sorted(missing))}"
        )

    runs: list[dict] = []
    seen_keys: set[tuple[str, str]] = set()

    for i, row in manifest.iterrows():
        ancestry = row["ancestry"].strip().upper()
        label = row["canonical_label_final"].strip().upper()
        run_dir = Path(row["run_dir"]).expanduser().resolve()

        if ancestry not in {"EUR", "EAS"}:
            raise ValueError(
                f"Manifest row {i + 2}: ancestry must be EUR or EAS, got {ancestry!r}"
            )
        if label not in CANONICAL_LABELS:
            raise ValueError(
                f"Manifest row {i + 2}: invalid canonical label {label!r}; "
                f"allowed: {sorted(CANONICAL_LABELS)}"
            )
        if label == "UNRESOLVED":
            raise ValueError(
                f"Manifest row {i + 2}: explicit manifests may not use UNRESOLVED"
            )
        if not run_dir.is_dir():
            raise FileNotFoundError(
                f"Manifest row {i + 2}: run directory not found: {run_dir}"
            )

        required_primary = [
            "magma.gsa.out",
            "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out",
            "magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out",
        ]
        absent = [name for name in required_primary if not (run_dir / name).is_file()]
        if absent:
            raise FileNotFoundError(
                f"Manifest row {i + 2}: {run_dir} is missing required files: "
                + ", ".join(absent)
            )

        key = (ancestry, label)
        if key in seen_keys:
            raise ValueError(
                f"Duplicate ancestry/model entry in manifest: {ancestry}/{label}"
            )
        seen_keys.add(key)

        raw_run = run_dir.name
        run_relative = row.get("run_relative_path", "").strip() or raw_run
        top_files = sorted(p.name for p in run_dir.iterdir() if p.is_file())
        params_text = read_small_text(run_dir / "params.config")

        runs.append(
            {
                "ancestry": ancestry,
                "run_label_raw": raw_run,
                "run_relative_path": run_relative,
                "run_key": f"{ancestry}/{label}",
                "run_dir": run_dir,
                "canonical_label_auto": label,
                "canonical_label_status_auto": "explicit_manifest",
                "canonical_label_reason_auto": "user-reviewed explicit run manifest",
                "canonical_label_final": label,
                "canonical_label_status_final": "explicit_manifest",
                "canonical_label_reason_final": "user-reviewed explicit run manifest",
                "top_level_files": "|".join(top_files),
                "params_config_compact": compact_text(params_text),
            }
        )

    return runs, []



def apply_manual_map(runs: list[dict], run_map_path: Path | None) -> None:
    if run_map_path is None:
        return
    if not run_map_path.exists():
        raise FileNotFoundError(f"Run-map file not found: {run_map_path}")

    mapping = pd.read_csv(run_map_path, sep="\t", dtype=str).fillna("")
    required = {"ancestry", "run_relative_path", "canonical_label_manual"}
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(
            f"Run-map is missing columns: {', '.join(sorted(missing))}"
        )

    lookup: dict[tuple[str, str], str] = {}
    for _, row in mapping.iterrows():
        label = row["canonical_label_manual"].strip().upper()
        if not label:
            continue
        if label not in CANONICAL_LABELS:
            raise ValueError(
                f"Invalid canonical label {label!r}; allowed: {sorted(CANONICAL_LABELS)}"
            )
        key = (row["ancestry"].strip(), row["run_relative_path"].strip())
        if key in lookup and lookup[key] != label:
            raise ValueError(f"Conflicting run-map labels for {key}")
        lookup[key] = label

    for run in runs:
        key = (run["ancestry"], run["run_relative_path"])
        if key in lookup:
            run["canonical_label_final"] = lookup[key]
            run["canonical_label_status_final"] = "manual_override"
            run["canonical_label_reason_final"] = str(run_map_path)


def check_outdir(outdir: Path) -> tuple[bool, str]:
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        probe = outdir / f".write_test_{os.getpid()}"
        probe.write_text("ok\n")
        probe.unlink()
        return True, ""
    except Exception as exc:
        return False, str(exc)


def read_magma_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            engine="python",
            comment="#",
            dtype=str,
            keep_default_na=False,
            na_values=["NA", "NaN", "nan", "."],
        )
    except Exception as exc:
        raise RuntimeError(f"Could not parse {path}: {exc}") from exc

    if df.empty:
        return df
    df = df.dropna(axis=1, how="all")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def find_column(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {str(c).strip().lower(): str(c) for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def bh_adjust(p: pd.Series) -> pd.Series:
    values = pd.to_numeric(p, errors="coerce").to_numpy(dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid_idx = np.where(np.isfinite(values))[0]
    if len(valid_idx) == 0:
        return pd.Series(result, index=p.index)

    valid_p = values[valid_idx]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result[valid_idx[order]] = adjusted
    return pd.Series(result, index=p.index)


def add_metadata(df: pd.DataFrame, run: dict, file_type: str, path: Path) -> pd.DataFrame:
    out = df.copy()
    metadata = {
        "ancestry": run["ancestry"],
        "run_label_raw": run["run_label_raw"],
        "run_relative_path": run["run_relative_path"],
        "run_key": run["run_key"],
        "canonical_label_auto": run["canonical_label_auto"],
        "canonical_label_status_auto": run["canonical_label_status_auto"],
        "canonical_label_final": run["canonical_label_final"],
        "canonical_label_status_final": run["canonical_label_status_final"],
        "result_type": file_type,
        "source_file": str(path),
    }
    for key, value in reversed(list(metadata.items())):
        out.insert(0, key, value)
    return out


def add_statistics(
    df: pd.DataFrame, alpha: float, expected_type: str | None
) -> tuple[pd.DataFrame, dict]:
    out = df.copy()
    p_col = find_column(out.columns, ["P", "P_VALUE", "PVALUE", "PVAL"])
    beta_col = find_column(out.columns, ["BETA"])
    type_col = find_column(out.columns, ["TYPE"])

    if expected_type and type_col:
        out["expected_magma_type"] = expected_type
        out["included_in_correction"] = (
            out[type_col].astype(str).str.upper() == expected_type
        )
    else:
        out["expected_magma_type"] = expected_type or ""
        out["included_in_correction"] = True

    summary = {
        "p_column": p_col or "",
        "beta_column": beta_col or "",
        "type_column": type_col or "",
        "expected_magma_type": expected_type or "",
        "n_rows_total": int(len(out)),
        "n_rows_in_correction": int(out["included_in_correction"].sum()),
        "n_valid_p": 0,
        "bonferroni_threshold": np.nan,
        "n_bonferroni_significant_positive": 0,
        "n_fdr_significant_positive": 0,
        "min_p": np.nan,
        "min_p_positive_beta": np.nan,
        "status": "parsed",
    }

    if p_col is None:
        summary["status"] = "parsed_no_p_column"
        return out, summary

    out["p_numeric"] = pd.to_numeric(out[p_col], errors="coerce")
    correction_mask = out["included_in_correction"] & out["p_numeric"].notna()
    n_tests = int(correction_mask.sum())
    threshold = alpha / n_tests if n_tests else np.nan

    out["n_tests_in_file"] = n_tests
    out["bonferroni_threshold"] = threshold
    out["p_bh"] = np.nan
    if n_tests:
        out.loc[correction_mask, "p_bh"] = bh_adjust(
            out.loc[correction_mask, "p_numeric"]
        ).to_numpy()
    out["significant_bonferroni"] = correction_mask & out["p_numeric"].le(
        threshold
    )
    out["significant_fdr"] = correction_mask & out["p_bh"].le(alpha)
    out["p_rank_within_file"] = np.nan
    if n_tests:
        out.loc[correction_mask, "p_rank_within_file"] = out.loc[
            correction_mask, "p_numeric"
        ].rank(method="min", ascending=True)

    if beta_col is not None:
        out["beta_numeric"] = pd.to_numeric(out[beta_col], errors="coerce")
        out["positive_beta"] = out["beta_numeric"].gt(0)
    else:
        out["beta_numeric"] = np.nan
        out["positive_beta"] = True

    out["significant_bonferroni_positive"] = (
        out["significant_bonferroni"] & out["positive_beta"]
    )
    out["significant_fdr_positive"] = out["significant_fdr"] & out["positive_beta"]

    positive_mask = correction_mask & out["positive_beta"]
    summary.update(
        {
            "n_valid_p": n_tests,
            "bonferroni_threshold": threshold,
            "n_bonferroni_significant_positive": int(
                out["significant_bonferroni_positive"].sum()
            ),
            "n_fdr_significant_positive": int(
                out["significant_fdr_positive"].sum()
            ),
            "min_p": (
                float(out.loc[correction_mask, "p_numeric"].min())
                if n_tests
                else np.nan
            ),
            "min_p_positive_beta": (
                float(out.loc[positive_mask, "p_numeric"].min())
                if positive_mask.any()
                else np.nan
            ),
        }
    )
    return out, summary


def safe_concat(frames: list[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", quoting=csv.QUOTE_MINIMAL)


def top_results(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    if df.empty or "p_numeric" not in df.columns:
        return pd.DataFrame()
    eligible = df[df["included_in_correction"].fillna(False)].copy()
    eligible = eligible.sort_values(
        ["ancestry", "run_key", "p_numeric"], na_position="last"
    )
    return eligible.groupby(["ancestry", "run_key"], group_keys=False).head(top_n)


def recurrence_table(frames: list[pd.DataFrame]) -> pd.DataFrame:
    combined = safe_concat(frames)
    if combined.empty or "VARIABLE" not in combined.columns:
        return pd.DataFrame()

    sig = combined[combined["significant_bonferroni_positive"].fillna(False)].copy()
    if sig.empty:
        return pd.DataFrame(
            columns=[
                "result_type",
                "VARIABLE",
                "n_significant_runs",
                "ancestries",
                "run_keys",
                "canonical_labels",
                "minimum_p",
            ]
        )

    rows = []
    for (result_type, variable), grp in sig.groupby(["result_type", "VARIABLE"]):
        rows.append(
            {
                "result_type": result_type,
                "VARIABLE": variable,
                "n_significant_runs": grp["run_key"].nunique(),
                "ancestries": ";".join(sorted(set(grp["ancestry"].astype(str)))),
                "run_keys": ";".join(sorted(set(grp["run_key"].astype(str)))),
                "canonical_labels": ";".join(
                    sorted(set(grp["canonical_label_final"].astype(str)))
                ),
                "minimum_p": pd.to_numeric(grp["p_numeric"], errors="coerce").min(),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["n_significant_runs", "minimum_p"], ascending=[False, True]
    )


def main() -> int:
    args = parse_args()
    base = args.base.expanduser().resolve() if args.base is not None else None
    run_manifest = (
        args.run_manifest.expanduser().resolve()
        if args.run_manifest is not None
        else None
    )
    outdir = args.outdir.expanduser().resolve()

    if base is None and run_manifest is None:
        print(
            "ERROR: supply either --base or --run-manifest",
            file=sys.stderr,
        )
        return 2
    if base is not None and run_manifest is not None:
        print(
            "ERROR: use either --base or --run-manifest, not both",
            file=sys.stderr,
        )
        return 2
    if base is not None and not base.exists():
        print(f"ERROR: base directory does not exist: {base}", file=sys.stderr)
        return 2

    writable, message = check_outdir(outdir)
    if not writable:
        fallback = Path.home() / outdir.name
        print(f"ERROR: cannot write to output directory: {outdir}", file=sys.stderr)
        print(f"Reason: {message}", file=sys.stderr)
        print(
            f"Use a writable location, for example: --outdir {fallback}",
            file=sys.stderr,
        )
        return 3

    try:
        explicit_manifest_mode = run_manifest is not None
        if explicit_manifest_mode:
            runs, skipped_runs = runs_from_manifest(run_manifest)
            if args.run_map is not None:
                raise ValueError(
                    "--run-map is not used with an explicit --run-manifest"
                )
        else:
            runs, skipped_runs = discover_runs(base)
            apply_manual_map(runs, args.run_map)

        multivariate_labels = {"FOURTRAIT", "MDD_3TRAIT", "SCZ_3TRAIT", "UNRESOLVED"}

        if explicit_manifest_mode:
            labels = {r["canonical_label_final"] for r in runs}
            single_labels = {"EA", "CF", "MDD", "SCZ"}
            pleio_labels = {"FOURTRAIT", "MDD_3TRAIT", "SCZ_3TRAIT"}

            if args.scope == "single-trait" and not labels.issubset(single_labels):
                raise ValueError(
                    "single-trait scope manifest contains non-single labels: "
                    + ", ".join(sorted(labels.difference(single_labels)))
                )
            if args.scope == "multivariate" and not labels.issubset(pleio_labels):
                raise ValueError(
                    "multivariate scope manifest contains non-PLEIO labels: "
                    + ", ".join(sorted(labels.difference(pleio_labels)))
                )

        elif args.scope == "single-trait":
            # Strict path allowlist: never infer single-trait scope from params.config.
            selected = []
            for r in runs:
                key = (r["ancestry"], r["run_relative_path"])
                if key not in SINGLE_TRAIT_RUN_ALLOWLIST:
                    continue
                expected_label = SINGLE_TRAIT_RUN_ALLOWLIST[key]
                r["canonical_label_final"] = expected_label
                r["canonical_label_status_final"] = "strict_path_allowlist"
                r["canonical_label_reason_final"] = (
                    f"exact single-trait directory allowlist: {r['ancestry']}/{r['run_relative_path']}"
                )
                selected.append(r)
            runs = selected

            forbidden = [
                r["run_key"]
                for r in runs
                if any(token in normalise(r["run_relative_path"])
                       for token in ("pleio", "4trait", "fourtrait"))
            ]
            if forbidden:
                raise RuntimeError(
                    "multivariate directory entered strict single-trait scope: "
                    + ", ".join(forbidden)
                )
        elif args.scope == "multivariate":
            # Strict path allowlist: analyse only the five known SNP2GENE PLEIO directories.
            selected = []
            for r in runs:
                key = (r["ancestry"], r["run_relative_path"])
                if key not in PLEIO_RUN_ALLOWLIST:
                    continue

                expected_label = PLEIO_RUN_ALLOWLIST[key]

                # A reviewed manual mapping is allowed only for the generic EAS/PLEIO run.
                if expected_label == "UNRESOLVED":
                    manual_label = r.get("canonical_label_final", "UNRESOLVED")
                    if manual_label in {"MDD_3TRAIT", "SCZ_3TRAIT"}:
                        r["canonical_label_status_final"] = "manual_reviewed"
                        r["canonical_label_reason_final"] = (
                            "reviewed manual mapping for generic EAS/PLEIO directory"
                        )
                    else:
                        r["canonical_label_final"] = "UNRESOLVED"
                        r["canonical_label_status_final"] = "strict_path_unresolved"
                        r["canonical_label_reason_final"] = (
                            "exact PLEIO directory allowlist; model identity requires --run-map"
                        )
                else:
                    r["canonical_label_final"] = expected_label
                    r["canonical_label_status_final"] = "strict_path_allowlist"
                    r["canonical_label_reason_final"] = (
                        f"exact PLEIO directory allowlist: "
                        f"{r['ancestry']}/{r['run_relative_path']}"
                    )
                selected.append(r)

            runs = selected

            forbidden = [
                r["run_key"]
                for r in runs
                if (r["ancestry"], r["run_relative_path"]) in SINGLE_TRAIT_RUN_ALLOWLIST
            ]
            if forbidden:
                raise RuntimeError(
                    "single-trait directory entered strict multivariate scope: "
                    + ", ".join(forbidden)
                )
    except Exception as exc:
        print(f"ERROR while constructing run map: {exc}", file=sys.stderr)
        return 4

    if not runs:
        location = run_manifest if run_manifest is not None else base
        print(f"ERROR: no target MAGMA runs found from {location}", file=sys.stderr)
        return 5

    # Mapping template: retain auto labels and leave manual column blank for deliberate review.
    mapping_template = pd.DataFrame(
        [
            {
                "ancestry": r["ancestry"],
                "run_relative_path": r["run_relative_path"],
                "run_label_raw": r["run_label_raw"],
                "canonical_label_auto": r["canonical_label_auto"],
                "canonical_label_status_auto": r["canonical_label_status_auto"],
                "canonical_label_reason_auto": r["canonical_label_reason_auto"],
                "canonical_label_manual": "",
                "notes": (
                    "REVIEW REQUIRED"
                    if r["canonical_label_final"] == "UNRESOLVED"
                    else ""
                ),
            }
            for r in runs
        ]
    )
    write_tsv(mapping_template, outdir / "00_run_mapping_template.tsv")

    run_audit = pd.DataFrame(
        [
            {
                k: v
                for k, v in r.items()
                if k not in {"run_dir"}
            }
            for r in runs
        ]
    )
    write_tsv(run_audit, outdir / "00_run_audit.tsv")
    selected_manifest = pd.DataFrame(
        [
            {
                "ancestry": r["ancestry"],
                "run_relative_path": r["run_relative_path"],
                "run_label_raw": r["run_label_raw"],
                "canonical_label_final": r["canonical_label_final"],
                "canonical_label_status_final": r["canonical_label_status_final"],
                "canonical_label_reason_final": r["canonical_label_reason_final"],
                "run_dir": str(r["run_dir"]),
            }
            for r in runs
        ]
    )
    write_tsv(selected_manifest, outdir / "00a_selected_run_manifest.tsv")
    write_tsv(
        pd.DataFrame(skipped_runs),
        outdir / "00b_skipped_nonstandard_run_dirs.tsv",
    )

    inventory_rows: list[dict] = []
    summaries: list[dict] = []
    extracted: dict[str, list[pd.DataFrame]] = {
        "competitive_gene_set": [],
        "gene_set_membership": [],
        "gtex_v8_specific_tissue": [],
        "gtex_v8_general_tissue": [],
    }

    for run in runs:
        run_dir: Path = run["run_dir"]
        for filename, file_type in FILE_TYPES.items():
            path = run_dir / filename
            inv = {
                "ancestry": run["ancestry"],
                "run_label_raw": run["run_label_raw"],
                "run_relative_path": run["run_relative_path"],
                "run_key": run["run_key"],
                "canonical_label_auto": run["canonical_label_auto"],
                "canonical_label_status_auto": run["canonical_label_status_auto"],
                "canonical_label_final": run["canonical_label_final"],
                "canonical_label_status_final": run["canonical_label_status_final"],
                "file_type": file_type,
                "filename": filename,
                "source_file": str(path),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
                "mtime_iso": (
                    pd.Timestamp(path.stat().st_mtime, unit="s", tz="UTC").isoformat()
                    if path.exists()
                    else ""
                ),
                "has_params_config": (run_dir / "params.config").exists(),
                "has_magma_genes_out": (run_dir / "magma.genes.out").exists(),
                "parse_status": "missing",
                "n_rows": 0,
                "columns": "",
            }

            if path.exists():
                try:
                    raw = read_magma_table(path)
                    inv["parse_status"] = "parsed"
                    inv["n_rows"] = len(raw)
                    inv["columns"] = "|".join(map(str, raw.columns))
                    meta = add_metadata(raw, run, file_type, path)

                    if file_type == "gene_set_membership":
                        extracted[file_type].append(meta)
                        stats = {
                            "p_column": "",
                            "beta_column": "",
                            "type_column": "",
                            "expected_magma_type": "",
                            "n_rows_total": len(meta),
                            "n_rows_in_correction": 0,
                            "n_valid_p": 0,
                            "bonferroni_threshold": np.nan,
                            "n_bonferroni_significant_positive": 0,
                            "n_fdr_significant_positive": 0,
                            "min_p": np.nan,
                            "min_p_positive_beta": np.nan,
                            "status": "membership_only",
                        }
                    else:
                        enriched, stats = add_statistics(
                            meta, args.alpha, EXPECTED_TYPE.get(file_type)
                        )
                        extracted[file_type].append(enriched)

                    summaries.append(
                        {
                            "ancestry": run["ancestry"],
                            "run_label_raw": run["run_label_raw"],
                            "run_relative_path": run["run_relative_path"],
                            "run_key": run["run_key"],
                            "canonical_label_auto": run["canonical_label_auto"],
                            "canonical_label_status_auto": run[
                                "canonical_label_status_auto"
                            ],
                            "canonical_label_final": run["canonical_label_final"],
                            "canonical_label_status_final": run[
                                "canonical_label_status_final"
                            ],
                            "result_type": file_type,
                            "source_file": str(path),
                            **stats,
                        }
                    )
                except Exception as exc:
                    inv["parse_status"] = f"ERROR: {exc}"
            inventory_rows.append(inv)

    inventory = pd.DataFrame(inventory_rows)
    if not inventory.empty:
        counts = (
            inventory[
                (inventory["file_type"] == "competitive_gene_set")
                & inventory["exists"]
            ]
            .groupby(
                ["ancestry", "canonical_label_final"], dropna=False
            )
            .size()
            .rename("n_run_dirs_per_ancestry_canonical")
            .reset_index()
        )
        inventory = inventory.merge(
            counts,
            on=["ancestry", "canonical_label_final"],
            how="left",
        )

    write_tsv(inventory, outdir / "00_file_inventory.tsv")
    if not inventory.empty:
        write_tsv(
            inventory[~inventory["exists"].fillna(False)].copy(),
            outdir / "00c_missing_source_files.tsv",
        )
    else:
        write_tsv(pd.DataFrame(), outdir / "00c_missing_source_files.tsv")

    summary_df = pd.DataFrame(summaries)
    write_tsv(summary_df, outdir / "01_run_level_summary.tsv")

    comp = safe_concat(extracted["competitive_gene_set"])
    membership = safe_concat(extracted["gene_set_membership"])
    gtex_specific = safe_concat(extracted["gtex_v8_specific_tissue"])
    gtex_general = safe_concat(extracted["gtex_v8_general_tissue"])

    write_tsv(comp, outdir / "02_MAGMA_competitive_gene_sets_all.tsv")
    if not comp.empty:
        write_tsv(
            comp[comp["significant_bonferroni_positive"].fillna(False)].copy(),
            outdir / "03_MAGMA_competitive_gene_sets_bonferroni_positive.tsv",
        )
        write_tsv(
            comp[comp["significant_fdr_positive"].fillna(False)].copy(),
            outdir / "04_MAGMA_competitive_gene_sets_FDR_positive.tsv",
        )
        write_tsv(
            top_results(comp, args.top_n),
            outdir / "04b_MAGMA_competitive_gene_sets_top_per_run.tsv",
        )

    write_tsv(membership, outdir / "05_MAGMA_gene_set_membership_all.tsv")

    # Restrict membership rows to Bonferroni-significant competitive sets when a
    # compatible set-name column is present. The all-membership table is always retained.
    significant_membership = pd.DataFrame()
    if not membership.empty and not comp.empty and "VARIABLE" in comp.columns:
        membership_set_col = find_column(
            membership.columns, ["VARIABLE", "SET", "GENE_SET", "FULL_NAME"]
        )
        if membership_set_col is not None:
            sig_sets = comp[
                comp["significant_bonferroni_positive"].fillna(False)
            ][["ancestry", "run_key", "VARIABLE"]].drop_duplicates()
            significant_membership = membership.merge(
                sig_sets,
                left_on=["ancestry", "run_key", membership_set_col],
                right_on=["ancestry", "run_key", "VARIABLE"],
                how="inner",
                suffixes=("", "_significant_set"),
            )
    write_tsv(
        significant_membership,
        outdir / "05b_MAGMA_bonferroni_significant_set_gene_membership.tsv",
    )

    write_tsv(gtex_specific, outdir / "06_GTEx_v8_specific_tissues_all.tsv")
    if not gtex_specific.empty:
        write_tsv(
            gtex_specific[
                gtex_specific["significant_bonferroni_positive"].fillna(False)
            ].copy(),
            outdir / "07_GTEx_v8_specific_tissues_bonferroni_positive.tsv",
        )
        write_tsv(
            gtex_specific[
                gtex_specific["significant_fdr_positive"].fillna(False)
            ].copy(),
            outdir / "08_GTEx_v8_specific_tissues_FDR_positive.tsv",
        )
        write_tsv(
            top_results(gtex_specific, args.top_n),
            outdir / "08b_GTEx_v8_specific_tissues_top_per_run.tsv",
        )

    write_tsv(gtex_general, outdir / "09_GTEx_v8_general_tissues_all.tsv")
    if not gtex_general.empty:
        write_tsv(
            gtex_general[
                gtex_general["significant_bonferroni_positive"].fillna(False)
            ].copy(),
            outdir / "10_GTEx_v8_general_tissues_bonferroni_positive.tsv",
        )
        write_tsv(
            gtex_general[
                gtex_general["significant_fdr_positive"].fillna(False)
            ].copy(),
            outdir / "11_GTEx_v8_general_tissues_FDR_positive.tsv",
        )
        write_tsv(
            top_results(gtex_general, args.top_n),
            outdir / "11b_GTEx_v8_general_tissues_top_per_run.tsv",
        )

    if args.scope == "single-trait":
        expected_labels = ["EA", "CF", "MDD", "SCZ"]
    elif args.scope == "multivariate":
        expected_labels = ["FOURTRAIT", "MDD_3TRAIT", "SCZ_3TRAIT"]
    else:
        expected_labels = EXPECTED_RUNS

    expected = pd.MultiIndex.from_product(
        [["EUR", "EAS"], expected_labels],
        names=["ancestry", "canonical_label_final"],
    ).to_frame(index=False)
    observed = (
        inventory[
            (inventory["file_type"] == "competitive_gene_set")
            & inventory["exists"]
            & (inventory["canonical_label_final"] != "UNRESOLVED")
        ][["ancestry", "canonical_label_final"]]
        .drop_duplicates()
        .assign(observed_magma_gsa=True)
    )
    expected_audit = expected.merge(
        observed,
        on=["ancestry", "canonical_label_final"],
        how="left",
    )
    expected_audit["observed_magma_gsa"] = expected_audit[
        "observed_magma_gsa"
    ].eq(True)
    write_tsv(expected_audit, outdir / "12_expected_run_audit.tsv")

    unresolved_df = run_audit[
        run_audit["canonical_label_final"].eq("UNRESOLVED")
    ].copy()
    write_tsv(unresolved_df, outdir / "12b_unresolved_run_labels.tsv")

    recurrence = recurrence_table([comp, gtex_specific, gtex_general])
    write_tsv(recurrence, outdir / "13_bonferroni_cross_run_recurrence.tsv")

    unresolved = [r for r in runs if r["canonical_label_final"] == "UNRESOLVED"]
    print(f"Wrote outputs to: {outdir}")
    print(f"Extraction scope: {args.scope}")
    print(f"Discovered direct SNP2GENE run directories in scope: {len(runs)}")
    if args.scope == "single-trait":
        print("Strict single-trait directories included:")
        for r in sorted(runs, key=lambda x: (x["ancestry"], x["run_relative_path"])):
            print(
                f"  - {r['ancestry']}/{r['run_relative_path']} "
                f"-> {r['canonical_label_final']}"
            )
    elif args.scope == "multivariate":
        print("Strict PLEIO SNP2GENE directories included:")
        for r in sorted(runs, key=lambda x: (x["ancestry"], x["run_relative_path"])):
            print(
                f"  - {r['ancestry']}/{r['run_relative_path']} "
                f"-> {r['canonical_label_final']}"
            )
    if skipped_runs:
        print(
            f"Skipped nonstandard/archive run directories: {len(skipped_runs)} "
            "(see 00b_skipped_nonstandard_run_dirs.tsv)"
        )
    if not inventory.empty:
        print(
            inventory.groupby(["file_type", "exists"])
            .size()
            .rename("n")
            .reset_index()
            .to_string(index=False)
        )
    if unresolved:
        print("\nNOTE: unresolved raw PLEIO labels remain; extraction is complete, but manuscript model naming still requires review:", file=sys.stderr)
        for run in unresolved:
            print(f"  - {run['run_key']}", file=sys.stderr)
        print(
            "Use 00_run_mapping_template.tsv only if canonical manuscript labels are needed; "
            "raw SNP2GENE results have already been extracted.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
