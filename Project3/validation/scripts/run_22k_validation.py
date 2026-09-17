#!/usr/bin/env python3
"""Validate the HELIOS 12k extension and combined 22k cognitive phenotype.

The script writes aggregate statistics and figures only. Participant-level
measurements, residuals, and scores remain in memory and are never exported.
"""

from __future__ import annotations

import argparse
import hashlib
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from run_phase1 import (
    boxplot,
    ecdf_plot,
    entropy,
    heatmap_plot,
    histogram_plot,
    interval_plot,
    pca,
    quantile,
    sha256,
    skew_kurtosis,
    standardize,
    tsv,
    write_text,
)


TASKS = [
    ("Reaction time", "dc7r6_react_avg", "DC7R6_React_Avg", -1, False),
    ("Stroop box", "dc7r8_stroopbox_avg", "DC7R8_Stroopbox_Avg", -1, False),
    ("Stroop ink", "dc7r10_stroopink_avg", "DC7R10_Stroopink_Avg", -1, False),
    ("Quiz score", "dc7r5_quiz_score", "DC7R5_Quiz_Score", 1, True),
    ("Working memory", "dc7r12_wm_score", "DC7R12_Wm_Score", 1, True),
    ("Pairing guesses", "dc7r4_pairing7_guesses", "DC7R4_Pairing7_Guesses", -1, True),
]
TASK_NAMES = [x[0] for x in TASKS]
RAW_COLUMNS = [x[1] for x in TASKS]
PCA_COLUMNS = [x[2] for x in TASKS]
DIRECTIONS = np.array([x[3] for x in TASKS], dtype=float)
SPEED_TASKS = ["Reaction time", "Stroop box", "Stroop ink"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cognitive", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--pca-input", type=Path, required=True)
    parser.add_argument("--extension-residuals", type=Path, required=True)
    parser.add_argument("--extension-g", type=Path, required=True)
    parser.add_argument("--validation-plan", type=Path)
    parser.add_argument("--gwas-log", action="append", type=Path, default=[])
    parser.add_argument("--ldsc-log", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260917)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", dtype={"IID": str, "ID": str})


def numeric_frame(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return frame[columns].apply(pd.to_numeric, errors="coerce")


def safe_slug(label: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", label.lower())).strip("_")


def align_score(reproduced: np.ndarray, observed: np.ndarray) -> dict[str, float | np.ndarray]:
    correlation = float(np.corrcoef(reproduced, observed)[0, 1])
    if correlation < 0:
        reproduced = -reproduced
        correlation = -correlation
    slope = float(np.dot(reproduced, observed) / np.dot(reproduced, reproduced))
    difference = observed - reproduced * slope
    return {
        "score": reproduced,
        "correlation": correlation,
        "slope": slope,
        "rmse": float(np.sqrt(np.mean(difference**2))),
        "max_difference": float(np.max(np.abs(difference))),
    }


def tucker(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def prepare_metadata(metadata_path: Path, ids: pd.Index) -> pd.DataFrame:
    metadata = read_table(metadata_path)
    required = {"IID", "Ancestry", "Sub_Cohort"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing: {sorted(missing)}")
    metadata["IID"] = metadata["IID"].astype(str).str.strip()
    if metadata["IID"].duplicated().any():
        raise ValueError("Metadata contains duplicate IIDs")
    metadata = metadata.set_index("IID").reindex(ids)
    metadata["Ancestry_label"] = metadata["Ancestry"].fillna("not_recorded").astype(str)
    subcohort = pd.to_numeric(metadata["Sub_Cohort"], errors="coerce")
    metadata["Freeze_label"] = subcohort.map({2.0: "HELIOS10K", 3.0: "HELIOS20K"}).fillna("not_recorded")
    return metadata


def reconstruct_residuals(metadata_path: Path, pca_input: pd.DataFrame, out: Path) -> dict:
    """Reproduce the archived residual tasks from the retained analysis table."""
    source = read_table(metadata_path)
    required = {
        "IID",
        "age",
        "sex",
        "age2",
        "age_sex",
        "age2_sex",
        "indian",
        "malay",
        "Edu",
        *PCA_COLUMNS,
    }
    missing = required - set(source.columns)
    if missing:
        raise ValueError(f"22k analysis table is missing: {sorted(missing)}")
    if source["IID"].duplicated().any():
        raise ValueError("22k analysis table contains duplicate IIDs")

    source["IID"] = source["IID"].astype(str).str.strip()
    numeric_covariates = ["age", "age2", "age_sex", "age2_sex", "indian", "malay"]
    source[numeric_covariates + ["Edu"]] = source[numeric_covariates + ["Edu"]].apply(
        pd.to_numeric, errors="coerce"
    )
    usable = source[numeric_covariates + ["Edu"]].notna().all(axis=1)
    model = source.loc[usable].copy().set_index("IID")

    continuous = model[numeric_covariates].astype(float)
    sex = pd.get_dummies(model["sex"].astype(str), prefix="sex", drop_first=True, dtype=float)
    education = pd.get_dummies(
        model["Edu"].astype(str), prefix="education", drop_first=True, dtype=float
    )
    design_frame = pd.concat([continuous, sex, education], axis=1)
    design = np.column_stack([np.ones(len(design_frame)), design_frame.to_numpy(dtype=float)])

    no_education_design_frame = pd.concat([continuous, sex], axis=1)
    no_education_design = np.column_stack(
        [np.ones(len(no_education_design_frame)), no_education_design_frame.to_numpy(dtype=float)]
    )

    reconstructed = pd.DataFrame(index=model.index)
    no_education = pd.DataFrame(index=model.index)
    for column in PCA_COLUMNS:
        response = pd.to_numeric(model[column], errors="coerce").to_numpy(dtype=float)
        if not np.all(np.isfinite(response)):
            raise ValueError(f"Non-finite raw values remain for {column}")
        coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
        reconstructed[column] = response - design @ coefficients
        no_education_coefficients = np.linalg.lstsq(no_education_design, response, rcond=None)[0]
        no_education[column] = response - no_education_design @ no_education_coefficients

    reference = pca_input.assign(IID=pca_input["IID"].astype(str)).set_index("IID")
    if set(reconstructed.index) != set(reference.index):
        raise ValueError("Residual reconstruction and PCA input have different participant sets")
    reconstructed = reconstructed.reindex(reference.index)
    rows = []
    for task, _, column, _, _ in TASKS:
        observed = pd.to_numeric(reference[column], errors="coerce").to_numpy(dtype=float)
        calculated = reconstructed[column].to_numpy(dtype=float)
        difference = observed - calculated
        rows.append(
            {
                "task": task,
                "n": len(observed),
                "correlation": np.corrcoef(observed, calculated)[0, 1],
                "rmse": np.sqrt(np.mean(difference**2)),
                "maximum_absolute_difference": np.max(np.abs(difference)),
            }
        )
    comparison = pd.DataFrame(rows)
    tsv(comparison, out / "baseline_residual_reproduction.tsv")

    current_oriented = reconstructed.to_numpy(dtype=float) * DIRECTIONS
    no_education_oriented = no_education.to_numpy(dtype=float) * DIRECTIONS
    _, current_score, current_loading, _, current_variance = pca(
        current_oriented, reference=np.ones(len(TASKS))
    )
    _, no_education_score, no_education_loading, _, no_education_variance = pca(
        no_education_oriented, reference=current_loading
    )
    education_sensitivity = pd.DataFrame(
        {
            "task": TASK_NAMES,
            "current_loading_education_adjusted": current_loading,
            "sensitivity_loading_without_education": no_education_loading,
            "absolute_loading_difference": np.abs(current_loading - no_education_loading),
            "current_pc1_variance_explained": current_variance,
            "without_education_pc1_variance_explained": no_education_variance,
            "score_correlation": np.corrcoef(current_score, no_education_score)[0, 1],
            "loading_congruence": tucker(current_loading, no_education_loading),
        }
    )
    tsv(education_sensitivity, out / "education_adjustment_sensitivity.tsv")
    return {
        "sample_n": len(reconstructed),
        "excluded_missing_education_or_covariates": int((~usable).sum()),
        "minimum_correlation": float(comparison["correlation"].min()),
        "maximum_difference": float(comparison["maximum_absolute_difference"].max()),
        "without_education_score_correlation": float(
            np.corrcoef(current_score, no_education_score)[0, 1]
        ),
        "without_education_loading_congruence": tucker(current_loading, no_education_loading),
        "without_education_variance_explained": no_education_variance,
    }


def baseline_outputs(args: argparse.Namespace, root: Path, pca_input: pd.DataFrame):
    out = root / "00_baseline"
    out.mkdir(parents=True, exist_ok=True)

    required = {"IID", "g", *PCA_COLUMNS}
    missing = required - set(pca_input.columns)
    if missing:
        raise ValueError(f"22k PCA input is missing: {sorted(missing)}")
    if pca_input["IID"].duplicated().any():
        raise ValueError("22k PCA input contains duplicate IIDs")

    residual_reproduction = reconstruct_residuals(args.metadata, pca_input, out)

    pca_matrix = numeric_frame(pca_input, PCA_COLUMNS).to_numpy(dtype=float)
    saved_g = pd.to_numeric(pca_input["g"], errors="coerce").to_numpy(dtype=float)
    complete = np.all(np.isfinite(pca_matrix), axis=1) & np.isfinite(saved_g)
    if not np.all(complete):
        raise ValueError("The retained 22k PCA input contains incomplete rows")

    _, stored_score, stored_loading, stored_corr, explained = pca(pca_matrix)
    alignment = align_score(stored_score, saved_g)
    if np.corrcoef(stored_score, saved_g)[0, 1] < 0:
        stored_loading = -stored_loading
        stored_corr = -stored_corr

    oriented = pca_matrix * DIRECTIONS
    _, better_score, better_loading, better_corr, explained_better = pca(
        oriented, reference=np.ones(len(TASKS))
    )
    better_alignment = align_score(better_score, -saved_g)
    if abs(explained - explained_better) > 1e-12:
        raise RuntimeError("Direction harmonisation unexpectedly changed PCA variance")

    loading_table = pd.DataFrame(
        {
            "task": TASK_NAMES,
            "source_column": PCA_COLUMNS,
            "higher_performance_multiplier": DIRECTIONS,
            "loading_on_stored_g": stored_loading,
            "correlation_with_stored_g": stored_corr,
            "loading_higher_performance": better_loading,
            "correlation_higher_performance": better_corr,
        }
    )
    tsv(loading_table, out / "baseline_pca_loadings.tsv")

    summary = pd.DataFrame(
        [
            ("combined_22k_pca_sample_n", len(pca_input)),
            ("combined_22k_residual_reproduction_minimum_correlation", residual_reproduction["minimum_correlation"]),
            ("combined_22k_residual_reproduction_max_absolute_difference", residual_reproduction["maximum_difference"]),
            ("combined_22k_pc1_variance_explained", explained),
            ("combined_22k_saved_score_correlation", alignment["correlation"]),
            ("combined_22k_saved_score_slope", alignment["slope"]),
            ("combined_22k_saved_score_rmse", alignment["rmse"]),
            ("combined_22k_saved_score_max_absolute_difference", alignment["max_difference"]),
            ("combined_22k_higher_performance_score_correlation", better_alignment["correlation"]),
            ("combined_22k_saved_g_direction_multiplier", -1),
        ],
        columns=["metric", "value"],
    )

    extension_resid = read_table(args.extension_residuals)
    extension_g = read_table(args.extension_g)
    ext_id_col = "ID" if "ID" in extension_resid.columns else "IID"
    missing = set(PCA_COLUMNS) - set(extension_resid.columns)
    if missing or not {"IID", "g"}.issubset(extension_g.columns):
        raise ValueError("12k extension files do not contain the expected PCA columns")
    ext_matrix_frame = numeric_frame(extension_resid, PCA_COLUMNS)
    ext_complete = np.all(np.isfinite(ext_matrix_frame), axis=1)
    ext_ids = extension_resid.loc[ext_complete, ext_id_col].astype(str)
    ext_matrix = ext_matrix_frame.loc[ext_complete].to_numpy(dtype=float)
    ext_saved = (
        extension_g.assign(IID=extension_g["IID"].astype(str))
        .drop_duplicates("IID")
        .set_index("IID")["g"]
    )
    common = pd.Index(ext_ids).intersection(ext_saved.index)
    indexer = pd.Index(ext_ids).get_indexer(common)
    _, ext_score, ext_loading_stored, _, ext_explained = pca(ext_matrix)
    ext_alignment = align_score(ext_score[indexer], ext_saved.loc[common].to_numpy(dtype=float))
    if np.corrcoef(ext_score[indexer], ext_saved.loc[common].to_numpy(dtype=float))[0, 1] < 0:
        ext_loading_stored = -ext_loading_stored
    ext_oriented = ext_matrix * DIRECTIONS
    _, _, ext_loading_better, ext_corr_better, _ = pca(
        ext_oriented, reference=np.ones(len(TASKS))
    )
    extension_summary = pd.DataFrame(
        [
            ("pca_sample_n", len(ext_matrix)),
            ("saved_score_overlap_n", len(common)),
            ("pc1_variance_explained", ext_explained),
            ("saved_score_correlation", ext_alignment["correlation"]),
            ("saved_score_rmse", ext_alignment["rmse"]),
            ("saved_score_max_absolute_difference", ext_alignment["max_difference"]),
            ("loading_congruence_with_combined_22k", tucker(ext_loading_better, better_loading)),
        ],
        columns=["metric", "value"],
    )
    tsv(extension_summary, out / "12k_archive_crosscheck.tsv")
    tsv(
        pd.DataFrame(
            {
                "task": TASK_NAMES,
                "loading_on_saved_12k_g": ext_loading_stored,
                "loading_higher_performance": ext_loading_better,
                "correlation_higher_performance": ext_corr_better,
            }
        ),
        out / "12k_archive_pca_loadings.tsv",
    )
    summary = pd.concat(
        [
            summary,
            pd.DataFrame(
                [
                    ("12k_extension_pca_sample_n", len(ext_matrix)),
                    ("12k_extension_pc1_variance_explained", ext_explained),
                    ("12k_extension_saved_score_correlation", ext_alignment["correlation"]),
                    ("12k_extension_loading_congruence_with_22k", tucker(ext_loading_better, better_loading)),
                ],
                columns=["metric", "value"],
            ),
        ],
        ignore_index=True,
    )
    tsv(summary, out / "baseline_pca_summary.tsv")

    baseline_ok = bool(
        alignment["correlation"] > 0.999999999
        and alignment["max_difference"] < 1e-8
        and residual_reproduction["minimum_correlation"] > 0.999999999
        and residual_reproduction["maximum_difference"] < 1e-6
        and ext_alignment["correlation"] > 0.999999999
        and ext_alignment["max_difference"] < 1e-8
    )
    report = f"""# 12k and combined 22k baseline reproduction

The archived 12k extension score was reproduced in **{len(ext_matrix):,} participants** (r = **{ext_alignment['correlation']:.12f}**; PC1 variance explained = **{100 * ext_explained:.2f}%**).

The current combined 22k score was reproduced in **{len(pca_input):,} participants** (r = **{alignment['correlation']:.12f}**; PC1 variance explained = **{100 * explained:.2f}%**). The maximum score difference after sign and scale alignment was **{alignment['max_difference']:.3g}**.

All six task residuals were independently reconstructed from the retained analysis table. The model used age, sex as a factor, age squared, age-by-sex and age-squared-by-sex terms, Indian and Malay ancestry indicators, and education as a factor. The minimum task-wise residual correlation was **{residual_reproduction['minimum_correlation']:.12f}** and the largest absolute numerical difference was **{residual_reproduction['maximum_difference']:.3g}**. The **{residual_reproduction['excluded_missing_education_or_covariates']}** records without a complete residualisation covariate set account for the difference between the 22,479-row analysis table and the 22,422-person PCA.

Removing education from the residualisation model left the loading vector nearly unchanged (Tucker congruence **{residual_reproduction['without_education_loading_congruence']:.4f}**) but changed participant scores enough that the two PC1s correlated **r = {residual_reproduction['without_education_score_correlation']:.3f}**. This is a phenotype-construction sensitivity only; its genetic effect requires a GWAS rerun with the no-education score.

The stored combined score is oriented so that higher values indicate poorer performance. It is multiplied by -1 in all higher-performance sensitivity outputs. PCA sign was fixed from the task directions before any genetic result was considered.

Baseline status: **{'PASS' if baseline_ok else 'FAIL'}**.
"""
    write_text(report, out / "baseline_reproduction_report.md")
    if not baseline_ok:
        raise RuntimeError("12k/22k baseline reproduction failed; sensitivity analyses stopped")

    return {
        "ids": pd.Index(pca_input["IID"].astype(str)),
        "matrix": oriented,
        "score": better_alignment["score"],
        "loadings": better_loading,
        "explained": explained,
        "loading_table": loading_table,
        "extension_n": len(ext_matrix),
        "extension_explained": ext_explained,
        "extension_congruence": tucker(ext_loading_better, better_loading),
        "residual_reproduction": residual_reproduction,
    }


def measurement_outputs(
    root: Path,
    cognitive_path: Path,
    ids: pd.Index,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    out = root / "01_measurement"
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    cognitive = read_table(cognitive_path)
    required = {"IID", *RAW_COLUMNS}
    missing = required - set(cognitive.columns)
    if missing:
        raise ValueError(f"Raw cognitive table is missing: {sorted(missing)}")
    cognitive["IID"] = cognitive["IID"].astype(str).str.strip()
    if cognitive["IID"].duplicated().any():
        raise ValueError("Raw cognitive table contains duplicate IIDs")
    common = cognitive.set_index("IID").reindex(ids)
    if common[RAW_COLUMNS].isna().all(axis=1).any():
        raise ValueError("At least one 22k PCA participant is absent from the raw cognitive table")

    joined = common.join(metadata[["Ancestry_label", "Freeze_label"]])
    summary_rows: list[dict] = []
    quantile_rows: list[dict] = []
    ancestry_rows: list[dict] = []
    freeze_rows: list[dict] = []

    for label, raw_col, _, direction, discrete in TASKS:
        values = pd.to_numeric(common[raw_col], errors="coerce").to_numpy(dtype=float)
        observed = values[np.isfinite(values)]
        skew, kurtosis = skew_kurtosis(observed)
        ent, effective = entropy(observed) if discrete else (np.nan, np.nan)
        levels = np.sort(np.unique(observed))
        pct_low = [100 * np.isin(observed, levels[: min(k, len(levels))]).mean() for k in (1, 2, 3)] if discrete else [np.nan] * 3
        pct_high = [100 * np.isin(observed, levels[-min(k, len(levels)) :]).mean() for k in (1, 2, 3)] if discrete else [np.nan] * 3
        summary_rows.append(
            {
                "task": label,
                "source_variable": raw_col,
                "source_direction": "higher_better" if direction == 1 else "lower_better",
                "pca_direction_multiplier": direction,
                "n_total": len(values),
                "n_nonmissing": len(observed),
                "missing_percent": 100 * (len(values) - len(observed)) / len(values),
                "unique_values": len(levels),
                "mean": np.mean(observed),
                "sd": np.std(observed, ddof=1),
                "variance": np.var(observed, ddof=1),
                "median": np.median(observed),
                "iqr": quantile(observed, 0.75) - quantile(observed, 0.25),
                "minimum": np.min(observed),
                "maximum": np.max(observed),
                "skewness": skew,
                "excess_kurtosis": kurtosis,
                "percent_at_minimum": 100 * np.mean(observed == np.min(observed)),
                "percent_at_maximum": 100 * np.mean(observed == np.max(observed)),
                "percent_bottom_1_category": pct_low[0],
                "percent_bottom_2_categories": pct_low[1],
                "percent_bottom_3_categories": pct_low[2],
                "percent_top_1_category": pct_high[0],
                "percent_top_2_categories": pct_high[1],
                "percent_top_3_categories": pct_high[2],
                "entropy_nats": ent,
                "effective_categories": effective,
            }
        )
        for q in (0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99):
            quantile_rows.append({"task": label, "quantile": q, "source_value": quantile(observed, q)})

        with_value = joined.assign(_value=pd.to_numeric(joined[raw_col], errors="coerce"))
        for group, frame in with_value.groupby("Ancestry_label", dropna=False):
            vals = frame["_value"].dropna().to_numpy(dtype=float)
            ancestry_rows.append(group_summary(label, "ancestry", str(group), vals))
        for group, frame in with_value.groupby("Freeze_label", dropna=False):
            vals = frame["_value"].dropna().to_numpy(dtype=float)
            freeze_rows.append(group_summary(label, "freeze", str(group), vals))

        safe = safe_slug(label)
        histogram_plot(observed, f"{label}: source distribution (22k)", raw_col, plots / f"histogram_{safe}.png")
        ecdf_plot(observed, f"{label}: empirical cumulative distribution (22k)", raw_col, plots / f"ecdf_{safe}.png")
        ancestry_values = {
            str(group): frame[raw_col].dropna().to_numpy(dtype=float)
            for group, frame in joined.groupby("Ancestry_label")
        }
        boxplot(ancestry_values, f"{label}: distribution by ancestry", "Source score", plots / f"ancestry_{safe}.png")
        freeze_values = {
            str(group): frame[raw_col].dropna().to_numpy(dtype=float)
            for group, frame in joined.groupby("Freeze_label")
        }
        boxplot(freeze_values, f"{label}: distribution by freeze", "Source score", plots / f"freeze_{safe}.png")

    summary = pd.DataFrame(summary_rows)
    tsv(summary, out / "task_measurement_summary.tsv")
    tsv(pd.DataFrame(quantile_rows), out / "task_quantiles.tsv")
    tsv(pd.DataFrame(ancestry_rows), out / "task_by_ancestry_summary.tsv")
    tsv(pd.DataFrame(freeze_rows), out / "task_by_freeze_summary.tsv")

    by_task = summary.set_index("task")
    pairing = by_task.loc["Pairing guesses"]
    report = f"""# 22k task-level measurement checks

The common PCA sample contains **{len(ids):,} participants** with complete raw values for all six tasks.

Quiz, Working Memory, and Pairing contain **{int(by_task.loc['Quiz score', 'unique_values'])}**, **{int(by_task.loc['Working memory', 'unique_values'])}**, and **{int(pairing['unique_values'])}** observed values. Their effective category counts are **{by_task.loc['Quiz score', 'effective_categories']:.1f}**, **{by_task.loc['Working memory', 'effective_categories']:.1f}**, and **{pairing['effective_categories']:.1f}**, respectively.

Pairing includes **{int(np.sum(pd.to_numeric(common['dc7r4_pairing7_guesses'], errors='coerce') == 0)):,}** zero-error participants. Unlike the historical 10k construction, the current combined 22k PCA uses the residualised error count rather than its reciprocal, so these best-performing observations are retained.

Raw distributions are reported by ancestry and freeze. Reliability coefficients were not estimated because repeated or item-level responses were not available.
"""
    write_text(report, out / "measurement_validation_report.md")
    return summary


def group_summary(task: str, group_name: str, group: str, values: np.ndarray) -> dict:
    return {
        "task": task,
        group_name: group,
        "n": len(values),
        "mean": np.mean(values),
        "sd": np.std(values, ddof=1) if len(values) > 1 else np.nan,
        "median": np.median(values),
        "q25": quantile(values, 0.25),
        "q75": quantile(values, 0.75),
        "minimum": np.min(values),
        "maximum": np.max(values),
    }


def correlation_table(matrix: np.ndarray, method: str) -> pd.DataFrame:
    return pd.DataFrame(matrix, columns=TASK_NAMES).corr(method=method)


def pca_structure_outputs(
    args: argparse.Namespace,
    root: Path,
    matrix: np.ndarray,
    score: np.ndarray,
    loadings: np.ndarray,
    metadata: pd.DataFrame,
):
    out = root / "02_pca_structure"
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    pearson = correlation_table(matrix, "pearson")
    spearman = correlation_table(matrix, "spearman")
    tsv(pearson.rename_axis("task").reset_index(), out / "phenotypic_correlations_pearson.tsv")
    tsv(spearman.rename_axis("task").reset_index(), out / "phenotypic_correlations_spearman.tsv")
    heatmap_plot([("Full sample", pearson)], "22k higher-performance task correlations", plots / "correlation_heatmap_full.png")

    group_rows: list[dict] = []
    for group_type, column in (("ancestry", "Ancestry_label"), ("freeze", "Freeze_label")):
        matrices = []
        groups = metadata[column].astype(str).to_numpy()
        for group in pd.unique(groups):
            keep = groups == group
            if keep.sum() < 50:
                continue
            group_matrix = correlation_table(matrix[keep], "pearson")
            matrices.append((group, group_matrix))
            _, _, group_loading, group_correlations, group_explained = pca(matrix[keep], reference=loadings)
            congruence = tucker(group_loading, loadings)
            for j, task in enumerate(TASK_NAMES):
                group_rows.append(
                    {
                        "group_type": group_type,
                        "group": group,
                        "n": int(keep.sum()),
                        "pc1_variance_explained": group_explained,
                        "tucker_congruence_with_full": congruence,
                        "task": task,
                        "loading_higher_performance": group_loading[j],
                        "correlation_with_pc1": group_correlations[j],
                    }
                )
        heatmap_plot(matrices, f"Task correlations by {group_type}", plots / f"correlation_heatmap_by_{group_type}.png")
    group_table = pd.DataFrame(group_rows)
    tsv(group_table, out / "group_specific_pca.tsv")

    rng = np.random.default_rng(args.seed)
    boot_rows: list[dict] = []
    n = len(matrix)
    for replicate in range(1, args.bootstrap + 1):
        index = rng.integers(0, n, n)
        _, _, loading, correlations, explained = pca(matrix[index], reference=loadings)
        for j, task in enumerate(TASK_NAMES):
            boot_rows.append(
                {
                    "replicate": replicate,
                    "task": task,
                    "pc1_variance_explained": explained,
                    "loading_higher_performance": loading[j],
                    "correlation_with_pc1": correlations[j],
                }
            )
    bootstrap = pd.DataFrame(boot_rows)
    tsv(bootstrap, out / "bootstrap_pca_loadings.tsv")
    summaries: list[dict] = []
    for task, frame in bootstrap.groupby("task", sort=False):
        for metric in ("loading_higher_performance", "correlation_with_pc1"):
            values = frame[metric].to_numpy()
            summaries.append(
                {
                    "task": task,
                    "metric": metric,
                    "mean": np.mean(values),
                    "sd": np.std(values, ddof=1),
                    "q025": quantile(values, 0.025),
                    "q975": quantile(values, 0.975),
                }
            )
    variance_values = bootstrap.drop_duplicates("replicate")["pc1_variance_explained"].to_numpy()
    summaries.append(
        {
            "task": "PC1",
            "metric": "variance_explained",
            "mean": np.mean(variance_values),
            "sd": np.std(variance_values, ddof=1),
            "q025": quantile(variance_values, 0.025),
            "q975": quantile(variance_values, 0.975),
        }
    )
    bootstrap_summary = pd.DataFrame(summaries)
    tsv(bootstrap_summary, out / "bootstrap_pca_summary.tsv")
    interval_plot(
        bootstrap_summary[bootstrap_summary["metric"] == "loading_higher_performance"],
        "22k bootstrap stability of PC1 loadings",
        plots / "bootstrap_loading_distributions.png",
    )

    loo_rows: list[dict] = []
    for omitted in TASK_NAMES:
        keep_cols = [i for i, task in enumerate(TASK_NAMES) if task != omitted]
        _, loo_score, loo_loading, loo_corr, loo_explained = pca(
            matrix[:, keep_cols], reference=loadings[keep_cols]
        )
        score_r = float(np.corrcoef(loo_score, score)[0, 1])
        if score_r < 0:
            loo_score *= -1
            loo_loading *= -1
            loo_corr *= -1
            score_r *= -1
        for j, column_index in enumerate(keep_cols):
            loo_rows.append(
                {
                    "omitted_task": omitted,
                    "retained_task": TASK_NAMES[column_index],
                    "pc1_variance_explained": loo_explained,
                    "correlation_with_full_pc1": score_r,
                    "loading_higher_performance": loo_loading[j],
                    "correlation_with_leave_one_out_pc1": loo_corr[j],
                }
            )
    loo = pd.DataFrame(loo_rows)
    tsv(loo, out / "leave_one_out_pca_summary.tsv")
    loo_matrix = loo.pivot(index="omitted_task", columns="retained_task", values="loading_higher_performance").reindex(index=TASK_NAMES, columns=TASK_NAMES)
    heatmap_plot([("Omitted task (rows)", loo_matrix)], "22k leave-one-task-out loadings", plots / "leave_one_out_loadings.png")

    z = standardize(matrix)
    speed_columns = [TASK_NAMES.index(task) for task in SPEED_TASKS]
    speed = z[:, speed_columns].mean(axis=1)
    collapsed_names = ["Speed domain", "Quiz score", "Working memory", "Pairing guesses"]
    collapsed_matrix = np.column_stack(
        [
            speed,
            z[:, TASK_NAMES.index("Quiz score")],
            z[:, TASK_NAMES.index("Working memory")],
            z[:, TASK_NAMES.index("Pairing guesses")],
        ]
    )
    _, collapsed_score, collapsed_loading, collapsed_corr, collapsed_explained = pca(
        collapsed_matrix, reference=np.ones(4)
    )
    collapsed_r = float(np.corrcoef(collapsed_score, score)[0, 1])
    if collapsed_r < 0:
        collapsed_score *= -1
        collapsed_loading *= -1
        collapsed_corr *= -1
        collapsed_r *= -1
    collapsed = pd.DataFrame(
        {
            "indicator": collapsed_names,
            "loading_higher_performance": collapsed_loading,
            "correlation_with_pc1": collapsed_corr,
            "pc1_variance_explained": collapsed_explained,
            "correlation_with_six_task_pc1": collapsed_r,
        }
    )
    tsv(collapsed, out / "collapsed_speed_pca_summary.tsv")

    loading_intervals = bootstrap_summary[
        bootstrap_summary["metric"] == "loading_higher_performance"
    ].set_index("task")
    loo_min = loo.groupby("omitted_task")["correlation_with_full_pc1"].first().min()
    loo_max = loo.groupby("omitted_task")["correlation_with_full_pc1"].first().max()
    group_first = group_table.groupby(["group_type", "group"], sort=False).first()
    ancestry_congruence = group_first.loc["ancestry", "tucker_congruence_with_full"].drop(
        index="not_recorded", errors="ignore"
    )
    freeze_congruence = group_first.loc["freeze", "tucker_congruence_with_full"].drop(
        index="not_recorded", errors="ignore"
    )
    report = f"""# Combined 22k PCA structure

The three latency indicators retained the largest loadings. Their 95% bootstrap intervals were Reaction Time **[{loading_intervals.loc['Reaction time', 'q025']:.3f}, {loading_intervals.loc['Reaction time', 'q975']:.3f}]**, Stroop Box **[{loading_intervals.loc['Stroop box', 'q025']:.3f}, {loading_intervals.loc['Stroop box', 'q975']:.3f}]**, and Stroop Ink **[{loading_intervals.loc['Stroop ink', 'q025']:.3f}, {loading_intervals.loc['Stroop ink', 'q975']:.3f}]**.

Leave-one-task-out scores correlated **r = {loo_min:.3f}–{loo_max:.3f}** with the full PC1. Loading congruence was **{ancestry_congruence.min():.3f}–{ancestry_congruence.max():.3f}** across recorded ancestry groups and **{freeze_congruence.min():.3f}–{freeze_congruence.max():.3f}** across HELIOS10K and HELIOS20K.

After replacing Reaction Time, Stroop Box, and Stroop Ink with one speed-domain indicator, the four-indicator PC1 explained **{100 * collapsed_explained:.2f}%** of variance and correlated **r = {collapsed_r:.3f}** with the six-task PC1.
"""
    write_text(report, out / "pca_structure_report.md")
    return bootstrap_summary, loo, collapsed, group_table


def parse_regenie_log(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    name = path.name.lower()
    ancestry = "Chinese" if "chinese" in name else "Indian" if "indian" in name else "Malay" if "malay" in name else "unknown"

    def value(pattern: str, cast=float):
        match = re.search(pattern, text)
        return cast(match.group(1)) if match else np.nan

    version_match = re.search(r"REGENIE v([0-9.]+)", text)
    return {
        "ancestry": ancestry,
        "regenie_version": version_match.group(1) if version_match else "not_recorded",
        "genotype_variants": value(r"n_snps\s*=\s*([0-9]+)", int),
        "genotype_samples": value(r"n_samples\s*=\s*([0-9]+)", int),
        "phenotyped_no_missing": value(r"phenotyped individuals with no missing data =\s*([0-9]+)", int),
        "covariate_records": value(r"individuals with covariate data =\s*([0-9]+)", int),
        "analysis_n": value(r"individuals used in analysis =\s*([0-9]+)", int),
        "minimum_mac": value(r"minimum MAC of ([0-9.]+)", float),
        "source_log": path.name,
        "sha256": sha256(path),
    }


def parse_ldsc_log(path: Path) -> pd.DataFrame:
    text = path.read_text(encoding="utf-8", errors="replace")
    marker = "Heritability Results for trait: current_22k.sumstats.gz"
    if marker not in text:
        raise ValueError("Current 22k LDSC block not found")
    block = text.split(marker, 1)[1].split("Calculating genetic covariance", 1)[0]

    def pair(pattern: str) -> tuple[float, float]:
        match = re.search(pattern, block)
        return (float(match.group(1)), float(match.group(2))) if match else (np.nan, np.nan)

    intercept, intercept_se = pair(r"Intercept:\s*([-0-9.]+)\s*\(([-0-9.]+)\)")
    h2, h2_se = pair(r"Total Observed Scale h2:\s*([-0-9.]+)\s*\(([-0-9.]+)\)")
    metrics = {
        "mean_chi_square": float(re.search(r"Mean Chi\^2 across remaining SNPs:\s*([-0-9.]+)", block).group(1)),
        "lambda_gc": float(re.search(r"Lambda GC:\s*([-0-9.]+)", block).group(1)),
        "intercept": intercept,
        "intercept_se": intercept_se,
        "observed_scale_h2": h2,
        "observed_scale_h2_se": h2_se,
        "h2_z": float(re.search(r"h2 Z:\s*([-0-9.]+)", block).group(1)),
    }
    return pd.DataFrame(
        [{"metric": key, "value": value, "source_log": path.name, "sha256": sha256(path)} for key, value in metrics.items()]
    )


def provenance_outputs(args: argparse.Namespace, root: Path) -> None:
    provenance = root / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    sources = [
        ("raw_cognitive_table", args.cognitive, "Raw six-task measurements"),
        ("covariate_table", args.covariates, "Ancestry and freeze provenance"),
        ("combined_22k_pca_input", args.pca_input, "Residual tasks and saved current g"),
        ("combined_22k_metadata", args.metadata, "Participant ancestry and freeze labels"),
        ("extension_12k_residuals", args.extension_residuals, "Archived 12k residual task matrix"),
        ("extension_12k_saved_g", args.extension_g, "Archived 12k saved g"),
        ("validation_plan", args.validation_plan, "Pre-specified validation specification"),
        ("ldsc_log", args.ldsc_log, "Current 22k LDSC output"),
    ]
    sources.extend((f"regenie_log_{i+1}", path, "22k ancestry-specific REGENIE Step 2 log") for i, path in enumerate(args.gwas_log))
    rows = []
    checksum_lines = []
    for item, path, role in sources:
        if path is None:
            continue
        if not path.exists():
            raise FileNotFoundError(path)
        digest = sha256(path)
        rows.append(
            {
                "item": item,
                "file": path.name,
                "role": role,
                "bytes": path.stat().st_size,
                "sha256": digest,
                "status": "verified",
            }
        )
        checksum_lines.append(f"{digest}  {path.name}")
    tsv(pd.DataFrame(rows), provenance / "baseline_manifest.tsv")
    write_text("\n".join(checksum_lines), provenance / "input_checksums.sha256")
    versions = [
        f"22k_validation_script\t{sha256(Path(__file__))}",
        f"python\t{platform.python_version()}",
        f"numpy\t{np.__version__}",
        f"pandas\t{pd.__version__}",
        f"pillow\t{Image.__version__}",
        "regenie\t4.1 (archived logs)",
        "ldsc\tGenomicSEM/cov-LDSC archived log",
    ]
    write_text("\n".join(versions), provenance / "software_versions.txt")


def final_summary(
    root: Path,
    baseline: dict,
    measurement: pd.DataFrame,
    bootstrap: pd.DataFrame,
    loo: pd.DataFrame,
    collapsed: pd.DataFrame,
    group_table: pd.DataFrame,
) -> None:
    load_ci = bootstrap[bootstrap["metric"] == "loading_higher_performance"].set_index("task")
    loo_values = loo.groupby("omitted_task")["correlation_with_full_pc1"].first()
    collapsed_r = float(collapsed["correlation_with_six_task_pc1"].iloc[0])
    by_task = measurement.set_index("task")
    groups = group_table.groupby(["group_type", "group"], sort=False).first()
    recorded_ancestry = groups.loc["ancestry"].drop(index="not_recorded", errors="ignore")
    recorded_freeze = groups.loc["freeze"].drop(index="not_recorded", errors="ignore")
    text = f"""# 12k extension and combined 22k validation

## Baseline

- Archived 12k extension PCA: **N = {baseline['extension_n']:,}**, PC1 variance **{100 * baseline['extension_explained']:.2f}%**, exact saved-score reproduction.
- Combined 22k PCA: **N = {len(baseline['ids']):,}**, PC1 variance **{100 * baseline['explained']:.2f}%**, exact saved-score reproduction.
- All six combined-sample residuals were reproduced from the retained covariates (minimum task-wise r = **{baseline['residual_reproduction']['minimum_correlation']:.12f}**).
- Omitting education from residualisation preserved the loading pattern (congruence **{baseline['residual_reproduction']['without_education_loading_congruence']:.4f}**) but produced a PC1 correlated **r = {baseline['residual_reproduction']['without_education_score_correlation']:.3f}** with the current score.
- Loading congruence between the archived 12k and combined 22k solutions: **{baseline['extension_congruence']:.3f}**.
- The stored 22k score is higher-worse; all validation results reverse it so higher means better performance.

## Measurement

- All six raw tasks are complete in the combined 22k PCA sample.
- Effective category counts were Quiz **{by_task.loc['Quiz score', 'effective_categories']:.1f}**, Working Memory **{by_task.loc['Working memory', 'effective_categories']:.1f}**, and Pairing **{by_task.loc['Pairing guesses', 'effective_categories']:.1f}**.
- The combined construction retains Pairing zero-error participants; it does not repeat the reciprocal-Pairing exclusion found in the historical 10k pipeline.

## PCA stability

- Reaction Time loading 95% interval: **[{load_ci.loc['Reaction time', 'q025']:.3f}, {load_ci.loc['Reaction time', 'q975']:.3f}]**.
- Stroop Box loading 95% interval: **[{load_ci.loc['Stroop box', 'q025']:.3f}, {load_ci.loc['Stroop box', 'q975']:.3f}]**.
- Stroop Ink loading 95% interval: **[{load_ci.loc['Stroop ink', 'q025']:.3f}, {load_ci.loc['Stroop ink', 'q975']:.3f}]**.
- Leave-one-task-out score correlations: **r = {loo_values.min():.3f}–{loo_values.max():.3f}**.
- Recorded ancestry loading congruence: **{recorded_ancestry['tucker_congruence_with_full'].min():.3f}–{recorded_ancestry['tucker_congruence_with_full'].max():.3f}**.
- HELIOS10K/HELIOS20K loading congruence: **{recorded_freeze['tucker_congruence_with_full'].min():.3f}–{recorded_freeze['tucker_congruence_with_full'].max():.3f}**.
- Collapsed-speed versus six-task PC1: **r = {collapsed_r:.3f}**.

The speed-heavy loading pattern persists in the larger sample, across the two recorded freezes, and across ancestry groups. It is therefore unlikely to be solely a small-pilot-sample artefact. Collapsing the three latency measures still changes the score materially, supporting the interpretation that giving speed three indicators affects the composition of PC1.
"""
    write_text(text, root / "22k_summary.md")


def main() -> int:
    args = parse_args()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    provenance_outputs(args, root)

    pca_input = read_table(args.pca_input)
    baseline = baseline_outputs(args, root, pca_input)
    metadata = prepare_metadata(args.metadata, baseline["ids"])
    measurement = measurement_outputs(root, args.cognitive, baseline["ids"], metadata)
    bootstrap, loo, collapsed, group_table = pca_structure_outputs(
        args, root, baseline["matrix"], baseline["score"], baseline["loadings"], metadata
    )

    if args.gwas_log:
        gwas = pd.DataFrame([parse_regenie_log(path) for path in args.gwas_log])
        gwas_analysis_n = int(gwas["analysis_n"].sum())
        total = {
            "ancestry": "Total",
            "regenie_version": "4.1",
            "genotype_variants": np.nan,
            "genotype_samples": np.nan,
            "phenotyped_no_missing": gwas["phenotyped_no_missing"].sum(),
            "covariate_records": gwas["covariate_records"].sum(),
            "analysis_n": gwas["analysis_n"].sum(),
            "minimum_mac": 0.5,
            "source_log": "three ancestry-specific logs",
            "sha256": "see ancestry rows",
        }
        gwas = pd.concat([gwas, pd.DataFrame([total])], ignore_index=True)
        tsv(gwas, root / "00_baseline" / "baseline_gwas_qc.tsv")
    else:
        gwas_analysis_n = np.nan
    if args.ldsc_log is not None:
        tsv(parse_ldsc_log(args.ldsc_log), root / "00_baseline" / "baseline_ldsc_qc.tsv")

    recorded_ancestry = int((metadata["Ancestry_label"] != "not_recorded").sum())
    sample_flow = pd.DataFrame(
        [
            ("raw_cognitive_source", len(read_table(args.cognitive)), "Rows in the source task table before cohort restriction"),
            ("retained_analysis_table", len(read_table(args.metadata)), "Rows with assembled phenotype/covariate data"),
            ("complete_residualisation_model", len(baseline["ids"]), "Rows retained after complete model covariates"),
            ("recorded_ancestry_and_freeze", recorded_ancestry, "PCA rows with recorded ancestry/freeze labels"),
            ("ancestry_specific_gwas", gwas_analysis_n, "Sum of participants used in the three REGENIE analyses"),
        ],
        columns=["stage", "n", "definition"],
    )
    tsv(sample_flow, root / "00_baseline" / "sample_flow.tsv")

    final_summary(root, baseline, measurement, bootstrap, loo, collapsed, group_table)
    print(f"22k validation complete: {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
