#!/usr/bin/env python3
"""
Compare model-wide FUMA SNP2GENE MAGMA/GTEx results from PLEIO models with
their constituent single-trait analyses.

Inputs are the output directories produced by extract_fuma_magma_gtex_v2_5.py
(or a compatible earlier extractor).

The comparison is term-level and ancestry-stratified. It compares:
  FOURTRAIT    versus EA, CF, MDD, SCZ
  MDD_3TRAIT   versus EA, CF, MDD
  SCZ_3TRAIT   versus EA, CF, SCZ

Important:
- EAS lacks a single-trait CF SNP2GENE run in the current project tree, so EAS
  classifications are explicitly marked incomplete.
- A generic EAS/PLEIO run labelled UNRESOLVED is retained for audit but is not
  assigned constituent traits until the extractor is rerun with a reviewed map.
- "PLEIO-only" means significant in PLEIO and not significant in any AVAILABLE
  constituent single-trait result. It does not prove biological novelty.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


SCRIPT_VERSION = "1.0.0"

RESULT_FILES = {
    "competitive_gene_set": "02_MAGMA_competitive_gene_sets_all.tsv",
    "gtex_v8_specific_tissue": "06_GTEx_v8_specific_tissues_all.tsv",
    "gtex_v8_general_tissue": "09_GTEx_v8_general_tissues_all.tsv",
}

MODEL_CONSTITUENTS = {
    "FOURTRAIT": ["EA", "CF", "MDD", "SCZ"],
    "MDD_3TRAIT": ["EA", "CF", "MDD"],
    "SCZ_3TRAIT": ["EA", "CF", "SCZ"],
}

PLEIO_MODELS = list(MODEL_CONSTITUENTS)
SINGLE_TRAITS = ["EA", "CF", "MDD", "SCZ"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    p.add_argument("--single-dir", required=True, type=Path)
    p.add_argument("--pleio-dir", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
        quoting=csv.QUOTE_MINIMAL,
    )


def read_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def bool_series(x: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False)
    return (
        x.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )


def normalise_result_table(df: pd.DataFrame, source_class: str) -> pd.DataFrame:
    required = {
        "ancestry",
        "canonical_label_final",
        "VARIABLE",
        "p_numeric",
        "beta_numeric",
        "significant_bonferroni_positive",
        "significant_fdr_positive",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"{source_class} table is missing columns: {', '.join(sorted(missing))}"
        )

    out = df.copy()
    out["source_class"] = source_class
    out["ancestry"] = out["ancestry"].astype(str)
    out["canonical_label_final"] = out["canonical_label_final"].astype(str)
    out["term"] = out["VARIABLE"].astype(str)

    for col in [
        "p_numeric",
        "beta_numeric",
        "BETA_STD",
        "SE",
        "NGENES",
        "bonferroni_threshold",
        "p_bh",
    ]:
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    out["significant_bonferroni_positive"] = bool_series(
        out["significant_bonferroni_positive"]
    )
    out["significant_fdr_positive"] = bool_series(
        out["significant_fdr_positive"]
    )

    if "included_in_correction" in out:
        out = out[bool_series(out["included_in_correction"])].copy()

    key = ["ancestry", "canonical_label_final", "term"]
    duplicate = out.duplicated(key, keep=False)
    if duplicate.any():
        examples = out.loc[duplicate, key].drop_duplicates().head(10)
        raise ValueError(
            f"Duplicate term keys in {source_class} table:\n"
            + examples.to_string(index=False)
        )
    return out


def compact_lookup(df: pd.DataFrame) -> pd.DataFrame:
    wanted = [
        "term",
        "FULL_NAME",
        "TYPE",
        "NGENES",
        "beta_numeric",
        "BETA_STD",
        "SE",
        "p_numeric",
        "p_bh",
        "bonferroni_threshold",
        "significant_bonferroni_positive",
        "significant_fdr_positive",
    ]
    cols = [c for c in wanted if c in df.columns]
    return df[cols].copy().set_index("term", drop=False)


def prefix_record(
    row: pd.Series | None,
    prefix: str,
    tested: bool,
) -> Dict[str, object]:
    fields = {
        f"{prefix}_tested": tested,
        f"{prefix}_n_genes": np.nan,
        f"{prefix}_beta": np.nan,
        f"{prefix}_beta_std": np.nan,
        f"{prefix}_se": np.nan,
        f"{prefix}_p": np.nan,
        f"{prefix}_p_bh": np.nan,
        f"{prefix}_bonf_threshold": np.nan,
        f"{prefix}_bonf_sig_positive": False,
        f"{prefix}_fdr_sig_positive": False,
    }
    if row is None:
        return fields

    mapping = {
        f"{prefix}_n_genes": "NGENES",
        f"{prefix}_beta": "beta_numeric",
        f"{prefix}_beta_std": "BETA_STD",
        f"{prefix}_se": "SE",
        f"{prefix}_p": "p_numeric",
        f"{prefix}_p_bh": "p_bh",
        f"{prefix}_bonf_threshold": "bonferroni_threshold",
        f"{prefix}_bonf_sig_positive": "significant_bonferroni_positive",
        f"{prefix}_fdr_sig_positive": "significant_fdr_positive",
    }
    for output_col, input_col in mapping.items():
        if input_col in row.index:
            value = row[input_col]
            if output_col.endswith("_sig_positive"):
                value = bool(value) if pd.notna(value) else False
            fields[output_col] = value
    return fields


def classify(pleio_sig: bool, any_single_sig: bool) -> str:
    if pleio_sig and any_single_sig:
        return "PLEIO_shared_with_available_constituent"
    if pleio_sig and not any_single_sig:
        return "PLEIO_only_among_available_constituents"
    if not pleio_sig and any_single_sig:
        return "constituent_only"
    return "neither"


def compare_one_model(
    single: pd.DataFrame,
    pleio: pd.DataFrame,
    ancestry: str,
    model: str,
    result_type: str,
) -> tuple[pd.DataFrame, dict, List[dict]]:
    constituents = MODEL_CONSTITUENTS[model]

    psub = pleio[
        (pleio["ancestry"] == ancestry)
        & (pleio["canonical_label_final"] == model)
    ].copy()
    if psub.empty:
        return pd.DataFrame(), {
            "ancestry": ancestry,
            "pleio_model": model,
            "result_type": result_type,
            "pleio_run_available": False,
            "comparison_complete": False,
            "missing_constituent_traits": ";".join(constituents),
        }, []

    p_lookup = compact_lookup(psub)

    single_lookups: Dict[str, pd.DataFrame] = {}
    available_traits: List[str] = []
    missing_traits: List[str] = []
    for trait in constituents:
        ssub = single[
            (single["ancestry"] == ancestry)
            & (single["canonical_label_final"] == trait)
        ].copy()
        if ssub.empty:
            missing_traits.append(trait)
        else:
            available_traits.append(trait)
            single_lookups[trait] = compact_lookup(ssub)

    terms = set(p_lookup.index)
    for lookup in single_lookups.values():
        terms.update(lookup.index)

    rows: List[dict] = []
    for term in sorted(terms):
        prow = p_lookup.loc[term] if term in p_lookup.index else None
        record: Dict[str, object] = {
            "ancestry": ancestry,
            "pleio_model": model,
            "result_type": result_type,
            "term": term,
            "full_name": (
                prow.get("FULL_NAME", term)
                if prow is not None
                else term
            ),
            "expected_constituent_traits": ";".join(constituents),
            "available_constituent_traits": ";".join(available_traits),
            "missing_constituent_traits": ";".join(missing_traits),
            "comparison_complete": len(missing_traits) == 0,
        }
        record.update(prefix_record(prow, "pleio", tested=prow is not None))

        single_bonf_traits: List[str] = []
        single_fdr_traits: List[str] = []
        for trait in constituents:
            lookup = single_lookups.get(trait)
            srow = lookup.loc[term] if lookup is not None and term in lookup.index else None
            record.update(prefix_record(srow, trait.lower(), tested=srow is not None))
            if srow is not None:
                if bool(srow.get("significant_bonferroni_positive", False)):
                    single_bonf_traits.append(trait)
                if bool(srow.get("significant_fdr_positive", False)):
                    single_fdr_traits.append(trait)

                pbeta = record.get("pleio_beta")
                sbeta = record.get(f"{trait.lower()}_beta")
                record[f"{trait.lower()}_same_beta_sign"] = (
                    bool(np.sign(pbeta) == np.sign(sbeta))
                    if pd.notna(pbeta) and pd.notna(sbeta) and pbeta != 0 and sbeta != 0
                    else np.nan
                )
            else:
                record[f"{trait.lower()}_same_beta_sign"] = np.nan

        pleio_bonf = bool(record["pleio_bonf_sig_positive"])
        pleio_fdr = bool(record["pleio_fdr_sig_positive"])
        any_single_bonf = len(single_bonf_traits) > 0
        any_single_fdr = len(single_fdr_traits) > 0

        record.update(
            {
                "n_available_constituent_traits": len(available_traits),
                "n_missing_constituent_traits": len(missing_traits),
                "n_constituent_bonf_sig": len(single_bonf_traits),
                "constituent_traits_bonf_sig": ";".join(single_bonf_traits),
                "any_constituent_bonf_sig": any_single_bonf,
                "classification_bonferroni": classify(pleio_bonf, any_single_bonf),
                "n_constituent_fdr_sig": len(single_fdr_traits),
                "constituent_traits_fdr_sig": ";".join(single_fdr_traits),
                "any_constituent_fdr_sig": any_single_fdr,
                "classification_fdr": classify(pleio_fdr, any_single_fdr),
            }
        )
        rows.append(record)

    comp = pd.DataFrame(rows)

    summary = {
        "ancestry": ancestry,
        "pleio_model": model,
        "result_type": result_type,
        "pleio_run_available": True,
        "comparison_complete": len(missing_traits) == 0,
        "available_constituent_traits": ";".join(available_traits),
        "missing_constituent_traits": ";".join(missing_traits),
        "n_terms_union": len(comp),
        "n_pleio_bonf_sig": int(comp["pleio_bonf_sig_positive"].sum()),
        "n_any_constituent_bonf_sig": int(comp["any_constituent_bonf_sig"].sum()),
        "n_pleio_only_available_bonf": int(
            (
                comp["classification_bonferroni"]
                == "PLEIO_only_among_available_constituents"
            ).sum()
        ),
        "n_pleio_shared_bonf": int(
            (
                comp["classification_bonferroni"]
                == "PLEIO_shared_with_available_constituent"
            ).sum()
        ),
        "n_constituent_only_bonf": int(
            (comp["classification_bonferroni"] == "constituent_only").sum()
        ),
        "n_pleio_fdr_sig": int(comp["pleio_fdr_sig_positive"].sum()),
        "n_any_constituent_fdr_sig": int(comp["any_constituent_fdr_sig"].sum()),
        "n_pleio_only_available_fdr": int(
            (
                comp["classification_fdr"]
                == "PLEIO_only_among_available_constituents"
            ).sum()
        ),
        "n_pleio_shared_fdr": int(
            (
                comp["classification_fdr"]
                == "PLEIO_shared_with_available_constituent"
            ).sum()
        ),
        "n_constituent_only_fdr": int(
            (comp["classification_fdr"] == "constituent_only").sum()
        ),
    }

    overlaps: List[dict] = []
    pleio_bonf_set = set(
        comp.loc[comp["pleio_bonf_sig_positive"], "term"]
    )
    pleio_fdr_set = set(
        comp.loc[comp["pleio_fdr_sig_positive"], "term"]
    )

    for trait in constituents:
        if trait not in available_traits:
            overlaps.append(
                {
                    "ancestry": ancestry,
                    "pleio_model": model,
                    "result_type": result_type,
                    "constituent_trait": trait,
                    "single_trait_available": False,
                }
            )
            continue

        t = trait.lower()
        single_bonf_set = set(
            comp.loc[comp[f"{t}_bonf_sig_positive"], "term"]
        )
        single_fdr_set = set(
            comp.loc[comp[f"{t}_fdr_sig_positive"], "term"]
        )

        def overlap_metrics(a: set, b: set, suffix: str) -> dict:
            intersection = a.intersection(b)
            union = a.union(b)
            return {
                f"n_pleio_sig_{suffix}": len(a),
                f"n_single_sig_{suffix}": len(b),
                f"n_exact_overlap_{suffix}": len(intersection),
                f"jaccard_{suffix}": len(intersection) / len(union) if union else np.nan,
                f"fraction_pleio_sig_shared_{suffix}": (
                    len(intersection) / len(a) if a else np.nan
                ),
                f"fraction_single_sig_recovered_{suffix}": (
                    len(intersection) / len(b) if b else np.nan
                ),
                f"overlap_terms_{suffix}": ";".join(sorted(intersection)),
            }

        row = {
            "ancestry": ancestry,
            "pleio_model": model,
            "result_type": result_type,
            "constituent_trait": trait,
            "single_trait_available": True,
        }
        row.update(overlap_metrics(pleio_bonf_set, single_bonf_set, "bonf"))
        row.update(overlap_metrics(pleio_fdr_set, single_fdr_set, "fdr"))
        overlaps.append(row)

    return comp, summary, overlaps


def pleio_recurrence(pleio: pd.DataFrame, result_type: str) -> pd.DataFrame:
    rows: List[dict] = []
    for (ancestry, term), group in pleio.groupby(["ancestry", "term"], sort=True):
        tested_models = sorted(group["canonical_label_final"].astype(str).unique())
        bonf_models = sorted(
            group.loc[
                group["significant_bonferroni_positive"],
                "canonical_label_final",
            ].astype(str).unique()
        )
        fdr_models = sorted(
            group.loc[
                group["significant_fdr_positive"],
                "canonical_label_final",
            ].astype(str).unique()
        )
        rows.append(
            {
                "ancestry": ancestry,
                "result_type": result_type,
                "term": term,
                "n_pleio_models_tested": len(tested_models),
                "pleio_models_tested": ";".join(tested_models),
                "n_pleio_models_bonf_sig": len(bonf_models),
                "pleio_models_bonf_sig": ";".join(bonf_models),
                "n_pleio_models_fdr_sig": len(fdr_models),
                "pleio_models_fdr_sig": ";".join(fdr_models),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    single_dir = args.single_dir.resolve()
    pleio_dir = args.pleio_dir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    audit_rows: List[dict] = []
    comparison_frames: List[pd.DataFrame] = []
    summaries: List[dict] = []
    overlap_rows: List[dict] = []
    recurrence_frames: List[pd.DataFrame] = []

    try:
        single_manifest = read_required(single_dir / "00a_selected_run_manifest.tsv")
        pleio_manifest = read_required(pleio_dir / "00a_selected_run_manifest.tsv")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    # Input audit.
    for source, manifest in [("single_trait", single_manifest), ("pleio", pleio_manifest)]:
        for _, row in manifest.iterrows():
            audit_rows.append(
                {
                    "source": source,
                    "ancestry": row.get("ancestry"),
                    "run_relative_path": row.get("run_relative_path"),
                    "canonical_label_final": row.get("canonical_label_final"),
                    "run_dir": row.get("run_dir"),
                }
            )
    write_tsv(pd.DataFrame(audit_rows), outdir / "00_input_run_audit.tsv")

    single_labels = set(single_manifest["canonical_label_final"].astype(str))
    bad_single = sorted(single_labels.difference(SINGLE_TRAITS))
    if bad_single:
        print(
            "ERROR: single-trait directory contains non-single labels: "
            + ", ".join(bad_single),
            file=sys.stderr,
        )
        return 3

    pleio_labels = set(pleio_manifest["canonical_label_final"].astype(str))
    bad_pleio = sorted(pleio_labels.difference(set(PLEIO_MODELS + ["UNRESOLVED"])))
    if bad_pleio:
        print(
            "ERROR: PLEIO directory contains unexpected labels: "
            + ", ".join(bad_pleio),
            file=sys.stderr,
        )
        return 4

    unresolved = pleio_manifest[
        pleio_manifest["canonical_label_final"].astype(str).eq("UNRESOLVED")
    ].copy()
    write_tsv(unresolved, outdir / "00b_unresolved_pleio_runs.tsv")

    for result_type, filename in RESULT_FILES.items():
        try:
            single_raw = normalise_result_table(
                read_required(single_dir / filename),
                source_class="single_trait",
            )
            pleio_raw = normalise_result_table(
                read_required(pleio_dir / filename),
                source_class="pleio",
            )
        except Exception as exc:
            print(f"ERROR loading {result_type}: {exc}", file=sys.stderr)
            return 5

        single = single_raw[single_raw["canonical_label_final"].isin(SINGLE_TRAITS)].copy()
        pleio_known = pleio_raw[pleio_raw["canonical_label_final"].isin(PLEIO_MODELS)].copy()

        for ancestry in ["EUR", "EAS"]:
            for model in PLEIO_MODELS:
                comp, summary, overlaps = compare_one_model(
                    single=single,
                    pleio=pleio_known,
                    ancestry=ancestry,
                    model=model,
                    result_type=result_type,
                )
                summaries.append(summary)
                overlap_rows.extend(overlaps)
                if not comp.empty:
                    comparison_frames.append(comp)

        recurrence_frames.append(pleio_recurrence(pleio_raw, result_type))

    all_comp = (
        pd.concat(comparison_frames, ignore_index=True, sort=False)
        if comparison_frames
        else pd.DataFrame()
    )
    write_tsv(all_comp, outdir / "01_all_term_comparisons.tsv")
    write_tsv(pd.DataFrame(summaries), outdir / "02_model_level_summary.tsv")
    write_tsv(pd.DataFrame(overlap_rows), outdir / "03_exact_overlap_by_trait.tsv")

    if not all_comp.empty:
        write_tsv(
            all_comp[
                all_comp["classification_bonferroni"]
                == "PLEIO_only_among_available_constituents"
            ].copy(),
            outdir / "04_PLEIO_only_among_available_constituents_bonferroni.tsv",
        )
        write_tsv(
            all_comp[
                all_comp["classification_bonferroni"]
                == "PLEIO_shared_with_available_constituent"
            ].copy(),
            outdir / "05_PLEIO_shared_with_constituent_bonferroni.tsv",
        )
        write_tsv(
            all_comp[
                all_comp["classification_bonferroni"] == "constituent_only"
            ].copy(),
            outdir / "06_constituent_only_bonferroni.tsv",
        )
        write_tsv(
            all_comp[
                all_comp["classification_fdr"]
                == "PLEIO_only_among_available_constituents"
            ].copy(),
            outdir / "07_PLEIO_only_among_available_constituents_FDR.tsv",
        )
        write_tsv(
            all_comp[
                all_comp["classification_fdr"]
                == "PLEIO_shared_with_available_constituent"
            ].copy(),
            outdir / "08_PLEIO_shared_with_constituent_FDR.tsv",
        )

    recurrence = (
        pd.concat(recurrence_frames, ignore_index=True, sort=False)
        if recurrence_frames
        else pd.DataFrame()
    )
    write_tsv(recurrence, outdir / "09_PLEIO_cross_model_recurrence.tsv")
    if not recurrence.empty:
        write_tsv(
            recurrence[recurrence["n_pleio_models_bonf_sig"] > 0].copy(),
            outdir / "10_PLEIO_bonferroni_significant_recurrence.tsv",
        )

    missing = pd.DataFrame(summaries)
    if not missing.empty:
        missing = missing[
            (~missing["pleio_run_available"].fillna(False))
            | (~missing["comparison_complete"].fillna(False))
        ].copy()
    write_tsv(missing, outdir / "11_missing_or_incomplete_comparisons.tsv")

    readme = f"""# MAGMA single-trait versus PLEIO comparison

Generated by `compare_magma_single_vs_pleio.py` version {SCRIPT_VERSION}.

## Classification

- `PLEIO_shared_with_available_constituent`: the term is significant in the
  PLEIO model and in at least one available constituent single-trait analysis.
- `PLEIO_only_among_available_constituents`: significant in PLEIO but not in
  any available constituent single-trait analysis.
- `constituent_only`: not significant in PLEIO but significant in at least one
  available constituent.
- `neither`: no corrected positive enrichment in either source.

The word **available** is essential. EAS currently has no single-trait CF
SNP2GENE result, so EAS comparisons involving CF are incomplete.

## Interpretation safeguards

1. PLEIO-only status is a descriptive significance-pattern comparison, not proof
   that a pathway is unique, causal, or newly created by multivariate analysis.
2. MAGMA P values and beta estimates should not be ranked directly across models
   as though they were measured on an identical statistical scale. Emphasise
   corrected-significance patterns, exact term overlap, and effect direction.
3. Gene-set databases are highly redundant; collapse related terms into themes
   after examining overlapping member genes.
4. These are model-wide PLEIO SNP2GENE results. They are distinct from the later
   directional Concordant/Discordant/DUAL Gene2Func and g:Profiler analyses.
5. Generic `EAS/PLEIO` results remain outside model-specific comparisons until
   assigned with a reviewed run map.
"""
    (outdir / "README_comparison.md").write_text(readme)

    print(f"Wrote comparison outputs to: {outdir}")
    print(f"Single-trait labels: {', '.join(sorted(single_labels))}")
    print(f"PLEIO labels: {', '.join(sorted(pleio_labels))}")
    if not unresolved.empty:
        print(
            "NOTE: unresolved PLEIO run(s) were retained for audit but excluded "
            "from model-specific comparisons:"
        )
        for _, row in unresolved.iterrows():
            print(f"  - {row['ancestry']}/{row['run_relative_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
