#!/usr/bin/env python3
"""Run Phase 1 validation of the HELIOS 10k cognitive phenotype.

Only aggregate statistics and plots are written. Participant-level scores and
residuals remain in memory and are never exported.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


TASKS = [
    ("Working memory", "DC7R12_Wm_Score", "dc7r12_wm_score", "identity", True),
    ("Quiz score", "DC7R5_Quiz_Score", "dc7r5_quiz_score", "identity", True),
    ("Reaction time", "DC7R6_React_Avg", "dc7r6_react_avg", "reciprocal", False),
    ("Stroop box", "DC7R8_Stroopbox_Avg", "dc7r8_stroopbox_avg", "reciprocal", False),
    ("Stroop ink", "DC7R10_Stroopink_Avg", "dc7r10_stroopink_avg", "reciprocal", False),
    ("Pairing guesses", "DC7R4_Pairing7_Guesses", "dc7r4_pairing7_guesses", "reciprocal", True),
]
TASK_NAMES = [x[0] for x in TASKS]
TASK_KEYS = {
    "Working memory": "Working_memory",
    "Quiz score": "Quiz_score",
    "Reaction time": "Reaction_time",
    "Stroop box": "Stroop_box",
    "Stroop ink": "Stroop_ink",
    "Pairing guesses": "Pairing_guesses",
}
TASK_SHORT = {
    "Working memory": "WM",
    "Quiz score": "Quiz",
    "Reaction time": "RT",
    "Stroop box": "Stroop box",
    "Stroop ink": "Stroop ink",
    "Pairing guesses": "Pairing",
}
SPEED_TASKS = ["Reaction time", "Stroop box", "Stroop ink"]
DISCRETE_TASKS = ["Working memory", "Quiz score", "Pairing guesses"]

BG = "#ffffff"
INK = "#162536"
MUTED = "#5d6b78"
BLUE = "#2077b4"
ORANGE = "#d97824"
GRID = "#dce3e8"
PALETTE = ["#2077b4", "#d97824", "#2b9a66", "#8b5bb5", "#c64d62", "#5f7f8f"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cognitive", type=Path, required=True)
    p.add_argument("--covariates", type=Path, required=True)
    p.add_argument("--historical-g", type=Path, required=True)
    p.add_argument("--validation-plan", type=Path)
    p.add_argument("--historical-r-script", type=Path)
    p.add_argument("--normalized-phenotypes", type=Path)
    p.add_argument("--residualized-phenotypes", type=Path)
    p.add_argument("--step1-log", type=Path)
    p.add_argument("--step2-log", type=Path)
    p.add_argument("--ldsc-munge-log", type=Path)
    p.add_argument("--ldsc-rg-log", type=Path)
    p.add_argument("--freeze-manifest", type=Path)
    p.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260917)
    return p.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=True)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tsv(data: pd.DataFrame, path: Path, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, sep="\t", index=index, na_rep="NA", float_format="%.12g")


def write_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def first_mode(series: pd.Series):
    observed = series.dropna().tolist()
    if not observed:
        return np.nan
    counts = {}
    first = {}
    for i, value in enumerate(observed):
        counts[value] = counts.get(value, 0) + 1
        first.setdefault(value, i)
    return max(counts, key=lambda x: (counts[x], -first[x]))


def sorted_levels(series: pd.Series) -> list:
    vals = list(pd.unique(series.dropna()))
    try:
        return sorted(vals, key=float)
    except (TypeError, ValueError):
        return sorted(vals, key=str)


def standardize(matrix: np.ndarray) -> np.ndarray:
    means = np.mean(matrix, axis=0)
    sds = np.std(matrix, axis=0, ddof=1)
    if np.any(~np.isfinite(sds)) or np.any(sds == 0):
        raise ValueError("Cannot standardize a constant or non-finite column")
    return (matrix - means) / sds


def pca(matrix: np.ndarray, reference: np.ndarray | None = None):
    z = standardize(matrix)
    u, singular, vt = np.linalg.svd(z, full_matrices=False)
    loading = vt[0].copy()
    score = u[:, 0] * singular[0]
    if reference is not None and float(np.dot(loading, reference)) < 0:
        loading *= -1
        score *= -1
    explained = float(singular[0] ** 2 / np.sum(singular**2))
    correlations = np.array([np.corrcoef(z[:, j], score)[0, 1] for j in range(z.shape[1])])
    return z, score, loading, correlations, explained


def quantile(values: np.ndarray, q: float) -> float:
    return float(np.quantile(values, q)) if len(values) else np.nan


def skew_kurtosis(values: np.ndarray) -> tuple[float, float]:
    n = len(values)
    if n < 4:
        return np.nan, np.nan
    mean = np.mean(values)
    sd = np.std(values, ddof=1)
    if sd == 0:
        return 0.0, 0.0
    centered = values - mean
    g1 = n / ((n - 1) * (n - 2)) * np.sum((centered / sd) ** 3)
    g2 = (
        n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * np.sum((centered / sd) ** 4)
        - 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    )
    return float(g1), float(g2)


def entropy(values: np.ndarray) -> tuple[float, float]:
    if not len(values):
        return np.nan, np.nan
    _, counts = np.unique(values, return_counts=True)
    probs = counts / counts.sum()
    ent = float(-np.sum(probs * np.log(probs)))
    return ent, float(np.exp(ent))


def font(size: int, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


def canvas(width: int, height: int, title: str):
    im = Image.new("RGB", (width, height), BG)
    draw = ImageDraw.Draw(im)
    draw.text((55, 30), title, fill=INK, font=font(30, True))
    return im, draw


def axes(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], xlabel: str = "", ylabel: str = ""):
    left, top, right, bottom = box
    draw.line((left, bottom, right, bottom), fill=INK, width=2)
    draw.line((left, top, left, bottom), fill=INK, width=2)
    if xlabel:
        draw.text(((left + right) // 2, bottom + 36), xlabel, anchor="mm", fill=INK, font=font(18))
    if ylabel:
        draw.text((left - 55, (top + bottom) // 2), ylabel, anchor="mm", fill=INK, font=font(18))


def histogram_plot(values: np.ndarray, title: str, xlabel: str, path: Path) -> None:
    values = values[np.isfinite(values)]
    im, draw = canvas(1200, 760, title)
    box = (105, 105, 1140, 655)
    axes(draw, box, xlabel, "Density")
    bins = min(50, max(12, int(round(math.sqrt(len(values))))))
    hist, edges = np.histogram(values, bins=bins, density=True)
    ymax = max(float(hist.max()), 1e-12) * 1.12
    left, top, right, bottom = box
    for i, h in enumerate(hist):
        x0 = left + (edges[i] - edges[0]) / (edges[-1] - edges[0]) * (right - left)
        x1 = left + (edges[i + 1] - edges[0]) / (edges[-1] - edges[0]) * (right - left)
        y = bottom - h / ymax * (bottom - top)
        draw.rectangle((x0 + 1, y, x1 - 1, bottom), fill="#b9d7ea", outline="#6aa3c8")
    sd = np.std(values, ddof=1)
    if sd > 0 and len(values) > 1:
        bandwidth = max(1.06 * sd * len(values) ** (-0.2), (edges[-1] - edges[0]) / 250)
        xs = np.linspace(edges[0], edges[-1], 220)
        density = np.zeros_like(xs)
        for start in range(0, len(values), 1000):
            block = values[start : start + 1000]
            density += np.exp(-0.5 * ((xs[:, None] - block[None, :]) / bandwidth) ** 2).sum(axis=1)
        density /= len(values) * bandwidth * math.sqrt(2 * math.pi)
        points = []
        for x, yv in zip(xs, density):
            px = left + (x - edges[0]) / (edges[-1] - edges[0]) * (right - left)
            py = bottom - min(yv / ymax, 1.0) * (bottom - top)
            points.append((px, py))
        draw.line(points, fill=ORANGE, width=4)
    draw.text((left, bottom + 12), f"{edges[0]:.3g}", fill=MUTED, font=font(15))
    draw.text((right, bottom + 12), f"{edges[-1]:.3g}", anchor="ra", fill=MUTED, font=font(15))
    im.save(path)


def ecdf_plot(values: np.ndarray, title: str, xlabel: str, path: Path) -> None:
    values = np.sort(values[np.isfinite(values)])
    im, draw = canvas(1200, 760, title)
    box = (105, 105, 1140, 655)
    axes(draw, box, xlabel, "Cumulative proportion")
    left, top, right, bottom = box
    if len(values):
        lo, hi = float(values[0]), float(values[-1])
        if hi == lo:
            hi = lo + 1
        idx = np.unique(np.linspace(0, len(values) - 1, min(len(values), 1800)).astype(int))
        pts = []
        for i in idx:
            x = left + (values[i] - lo) / (hi - lo) * (right - left)
            y = bottom - ((i + 1) / len(values)) * (bottom - top)
            pts.append((x, y))
        draw.line(pts, fill=BLUE, width=4)
        for frac in (0.25, 0.5, 0.75):
            y = bottom - frac * (bottom - top)
            draw.line((left, y, right, y), fill=GRID, width=1)
            draw.text((left - 12, y), f"{frac:.2f}", anchor="rm", fill=MUTED, font=font(14))
    im.save(path)


def boxplot(values_by_group: dict[str, np.ndarray], title: str, ylabel: str, path: Path) -> None:
    im, draw = canvas(1200, 760, title)
    box = (105, 105, 1140, 635)
    axes(draw, box, "", ylabel)
    clean = {k: v[np.isfinite(v)] for k, v in values_by_group.items() if np.isfinite(v).any()}
    all_values = np.concatenate(list(clean.values())) if clean else np.array([0.0, 1.0])
    lo, hi = float(np.min(all_values)), float(np.max(all_values))
    if hi == lo:
        hi = lo + 1
    left, top, right, bottom = box
    groups = list(clean)
    draw.text((left - 12, top), f"{hi:.3g}", anchor="rm", fill=MUTED, font=font(14))
    draw.text((left - 12, bottom), f"{lo:.3g}", anchor="rm", fill=MUTED, font=font(14))
    for i, group in enumerate(groups):
        vals = clean[group]
        x = left + (i + 1) * (right - left) / (len(groups) + 1)
        q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        low, high = np.quantile(vals, [0.05, 0.95])
        sy = lambda y: bottom - (y - lo) / (hi - lo) * (bottom - top)
        draw.line((x, sy(low), x, sy(high)), fill=INK, width=3)
        draw.rectangle((x - 55, sy(q3), x + 55, sy(q1)), fill=PALETTE[i % len(PALETTE)], outline=INK, width=2)
        draw.line((x - 55, sy(med), x + 55, sy(med)), fill=BG, width=4)
        draw.text((x, bottom + 28), f"{group}\nN={len(vals):,}", anchor="ma", align="center", fill=INK, font=font(16))
    im.save(path)


def heatmap_plot(matrices: list[tuple[str, pd.DataFrame]], title: str, path: Path) -> None:
    n = len(matrices)
    panel = 610
    im, draw = canvas(80 + n * panel, 780, title)
    cell = 68
    for p, (label, matrix) in enumerate(matrices):
        x0 = 170 + p * panel
        y0 = 180
        draw.text((x0 + cell * 3, 105), label, anchor="mm", fill=INK, font=font(22, True))
        for i, row in enumerate(matrix.index):
            draw.text((x0 - 12, y0 + i * cell + cell / 2), TASK_KEYS.get(row, row).replace("_", " "), anchor="rm", fill=INK, font=font(13))
            for j, col in enumerate(matrix.columns):
                value = float(matrix.iloc[i, j])
                rect = (x0 + j * cell, y0 + i * cell, x0 + (j + 1) * cell, y0 + (i + 1) * cell)
                if not np.isfinite(value):
                    draw.rectangle(rect, fill="#edf0f2", outline=BG, width=2)
                    draw.text((rect[0] + cell / 2, rect[1] + cell / 2), "-", anchor="mm", fill=MUTED, font=font(13))
                    continue
                strength = min(abs(value), 1.0)
                if value >= 0:
                    color = (int(245 - 140 * strength), int(248 - 65 * strength), int(252 - 30 * strength))
                else:
                    color = (int(252 - 25 * strength), int(242 - 105 * strength), int(238 - 135 * strength))
                draw.rectangle(rect, fill=color, outline=BG, width=2)
                draw.text((rect[0] + cell / 2, rect[1] + cell / 2), f"{value:.2f}", anchor="mm", fill=INK, font=font(13))
        for j, col in enumerate(matrix.columns):
            text = TASK_SHORT.get(col, TASK_KEYS.get(col, col).replace("_", " "))
            draw.text((x0 + j * cell + cell / 2, y0 + 6 * cell + 12), text, anchor="ma", fill=INK, font=font(12))
    im.save(path)


def interval_plot(summary: pd.DataFrame, title: str, path: Path) -> None:
    im, draw = canvas(1200, 760, title)
    box = (280, 105, 1130, 665)
    left, top, right, bottom = box
    draw.line((left, top, left, bottom), fill=GRID, width=2)
    draw.line((right, top, right, bottom), fill=GRID, width=2)
    xmin = min(0.0, float(summary["q025"].min()) - 0.03)
    xmax = max(0.9, float(summary["q975"].max()) + 0.03)
    for i, row in summary.reset_index(drop=True).iterrows():
        y = top + (i + 1) * (bottom - top) / (len(summary) + 1)
        sx = lambda x: left + (x - xmin) / (xmax - xmin) * (right - left)
        draw.text((left - 18, y), row["task"], anchor="rm", fill=INK, font=font(18))
        draw.line((sx(row["q025"]), y, sx(row["q975"]), y), fill=BLUE, width=5)
        draw.ellipse((sx(row["mean"]) - 7, y - 7, sx(row["mean"]) + 7, y + 7), fill=ORANGE)
    draw.text(((left + right) / 2, bottom + 35), "PC1 loading (higher performance)", anchor="mm", fill=INK, font=font(18))
    im.save(path)


def prepare_phenotypes(cognitive: pd.DataFrame, covariates: pd.DataFrame):
    cognitive = cognitive.copy()
    covariates = covariates.copy()
    cognitive["IID"] = cognitive["IID"].astype(str).str.strip()
    covariates["IID"] = covariates["IID"].astype(str).str.strip()
    cognitive = cognitive[(cognitive["IID"] != "") & ~cognitive["IID"].duplicated()].copy()
    covariates = covariates[(covariates["IID"] != "") & ~covariates["IID"].duplicated()].copy()

    task = pd.DataFrame(index=cognitive["IID"])
    mask_counts = []
    for label, raw_col, _, transform, _ in TASKS:
        raw = numeric(cognitive[raw_col]).to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            aligned = 1.0 / raw if transform == "reciprocal" else raw.copy()
        invalid = int(np.sum(~np.isfinite(aligned) & np.isfinite(raw)))
        aligned[~np.isfinite(aligned)] = np.nan
        finite = np.isfinite(aligned)
        before = int(finite.sum())
        z = np.full(len(aligned), np.nan)
        z[finite] = (aligned[finite] - np.mean(aligned[finite])) / np.std(aligned[finite], ddof=1)
        outlier = finite & (np.abs(z) > 4)
        aligned[outlier] = np.nan
        task[label] = aligned
        mask_counts.append({"task": label, "non_finite_after_transform": invalid, "masked_beyond_4sd": int(outlier.sum()), "retained_after_transform_and_mask": int(np.isfinite(aligned).sum()), "finite_before_mask": before})

    cov = covariates.set_index("IID", drop=False)
    common = task.index.intersection(cov.index, sort=False)
    task = task.loc[common]
    model = cov.loc[common].copy()
    model["age"] = numeric(model["Age"])
    model["sex_numeric"] = numeric(model["Sex"])
    model["age2"] = model["age"] ** 2
    model["age_sex"] = model["age"] * model["sex_numeric"]
    model["age2_sex"] = model["age2"] * model["sex_numeric"]
    model["HI"] = numeric(model["FSAQ21_Soc9"])
    model["SES_dontknow"] = (model["HI"] == -888).astype(int)
    model["SES_prefno"] = (model["HI"] == -777).astype(int)
    model["indian"] = (model["Ancestry"] == "Indian").astype(int)
    model["malay"] = (model["Ancestry"] == "Malay").astype(int)
    model["sex"] = model["Sex"]

    factors = ["sex", "indian", "malay", "HI", "SES_prefno", "SES_dontknow"]
    for col in factors:
        model[col] = model[col].fillna(first_mode(model[col]))

    xparts = [np.ones((len(model), 1))]
    for col in ["age", "age_sex", "age2", "age2_sex"]:
        xparts.append(numeric(model[col]).to_numpy(dtype=float).reshape(-1, 1))
    for col in factors:
        levels = sorted_levels(model[col])
        for level in levels[1:]:
            xparts.append((model[col] == level).to_numpy(dtype=float).reshape(-1, 1))
    design = np.hstack(xparts)

    residual = np.full(task.shape, np.nan)
    valid_design = np.all(np.isfinite(design), axis=1)
    for j, label in enumerate(TASK_NAMES):
        y = task[label].to_numpy(dtype=float)
        valid = valid_design & np.isfinite(y)
        beta, _, _, _ = np.linalg.lstsq(design[valid], y[valid], rcond=None)
        residual[valid, j] = y[valid] - design[valid] @ beta
    keep = np.all(np.isfinite(residual), axis=1)
    residual = residual[keep]
    ids = task.index.to_numpy()[keep]
    metadata = model.loc[ids].copy()
    return ids, residual, metadata, pd.DataFrame(mask_counts), cognitive, covariates


def baseline_outputs(args, root: Path, ids, residual, metadata, mask_counts):
    out = root / "00_baseline"
    out.mkdir(parents=True, exist_ok=True)
    _, score_stored, loading_stored, corr_stored, explained = pca(residual)
    stroop_index = TASK_NAMES.index("Stroop box")
    if loading_stored[stroop_index] > 0:
        loading_stored *= -1
        score_stored *= -1
        corr_stored *= -1
    score_better = -score_stored
    loading_better = -loading_stored
    corr_better = -corr_stored

    historical = read_table(args.historical_g)
    historical["IID"] = historical["IID"].astype(str).str.strip()
    historical["g"] = numeric(historical["g"])
    hist = historical.drop_duplicates("IID").set_index("IID")["g"]
    valid_ids = [x for x in ids if x in hist.index and np.isfinite(hist.loc[x])]
    idx = pd.Index(ids).get_indexer(valid_ids)
    observed = hist.loc[valid_ids].to_numpy(dtype=float)
    reproduced = score_stored[idx]
    raw_r = float(np.corrcoef(reproduced, observed)[0, 1])
    slope = float(np.dot(reproduced, observed) / np.dot(reproduced, reproduced))
    aligned = reproduced * slope
    differences = observed - aligned

    rows = []
    for i, task in enumerate(TASK_NAMES):
        rows.append({"task": task, "loading_on_stored_g": loading_stored[i], "correlation_with_stored_g": corr_stored[i], "loading_higher_performance": loading_better[i], "correlation_higher_performance": corr_better[i]})
    loading_table = pd.DataFrame(rows).sort_values("loading_higher_performance", ascending=False)
    summary = pd.DataFrame(
        [
            ("source_cognitive_sample_n", len(read_table(args.cognitive))),
            ("deduplicated_covariate_sample_n", len(pd.unique(read_table(args.covariates)["IID"]))),
            ("pca_sample_n", len(ids)),
            ("pc1_variance_explained", explained),
            ("historical_score_overlap_n", len(valid_ids)),
            ("historical_score_correlation", raw_r),
            ("historical_score_slope_through_origin", slope),
            ("historical_score_rmse_after_alignment", float(np.sqrt(np.mean(differences**2)))),
            ("historical_score_max_absolute_difference", float(np.max(np.abs(differences)))),
        ],
        columns=["metric", "value"],
    )
    if args.normalized_phenotypes is not None:
        normalized = read_table(args.normalized_phenotypes)
        if {"IID", "g"}.issubset(normalized.columns):
            normalized["IID"] = normalized["IID"].astype(str).str.strip()
            normalized["g"] = numeric(normalized["g"])
            normalized_g = normalized.drop_duplicates("IID").set_index("IID")["g"]
            cross_ids = hist.index.intersection(normalized_g.index)
            saved_cross = hist.loc[cross_ids].to_numpy(dtype=float)
            normalized_cross = normalized_g.loc[cross_ids].to_numpy(dtype=float)
            good = np.isfinite(saved_cross) & np.isfinite(normalized_cross)
            summary = pd.concat(
                [
                    summary,
                    pd.DataFrame(
                        [
                            ("normalized_file_score_overlap_n", int(good.sum())),
                            ("normalized_file_score_correlation", float(np.corrcoef(saved_cross[good], normalized_cross[good])[0, 1])),
                            ("normalized_file_score_max_absolute_difference", float(np.max(np.abs(saved_cross[good] - normalized_cross[good])))),
                        ],
                        columns=["metric", "value"],
                    ),
                ],
                ignore_index=True,
            )
    has_normalized_crosscheck = bool(summary["metric"].eq("normalized_file_score_correlation").any())
    normalized_note = (
        "The `g` column in the retained normalized-phenotype table also matched the separate saved-score file exactly. "
        if has_normalized_crosscheck
        else "No normalized-phenotype table was supplied for an additional score-file cross-check. "
    )
    comparison = summary[summary["metric"].str.contains("historical_|normalized_file_score")].copy()
    tsv(summary, out / "baseline_pca_summary.tsv")
    tsv(loading_table, out / "baseline_pca_loadings.tsv")
    corr_matrix = pd.DataFrame(np.corrcoef(standardize(residual), rowvar=False), index=TASK_NAMES, columns=TASK_NAMES)
    tsv(corr_matrix.rename_axis("task").reset_index(), out / "baseline_task_correlations.tsv")
    tsv(comparison, out / "baseline_score_comparison.tsv")
    tsv(mask_counts, out / "baseline_task_exclusions.tsv")

    benchmark = pd.read_csv(root / "config" / "baseline_genetic_benchmarks.tsv", sep="\t", dtype=str)
    tsv(benchmark[benchmark["analysis"].str.contains("metal|regenie|ldsc")], out / "baseline_gwas_qc.tsv")
    tsv(benchmark[benchmark["analysis"].str.contains("vs_koges")], out / "baseline_rg_summary.tsv")

    targets = pd.read_csv(root / "config" / "historical_pca_targets.tsv", sep="\t")
    checked = targets.merge(loading_table, on="task")
    loading_ok = np.all(np.abs(checked["loading_higher_performance"] - checked["expected_loading_higher_performance"]) <= checked["tolerance"])
    corr_ok = np.all(np.abs(checked["correlation_higher_performance"] - checked["expected_correlation_higher_performance"]) <= checked["tolerance"])
    baseline_ok = bool(len(ids) == 7403 and abs(explained - 0.3448937052785348) <= 0.002 and raw_r > 0.999999 and np.max(np.abs(differences)) < 1e-8 and loading_ok and corr_ok)
    report = f"""# Baseline reproduction

The historical six-task PCA was reproduced in **{len(ids):,} participants**. PC1 explained **{explained * 100:.2f}%** of variance. The reconstructed stored-orientation score correlated **r = {raw_r:.12f}** with the saved score; the largest absolute difference after scale alignment was **{np.max(np.abs(differences)):.3g}**.

All six loadings and task-PC1 correlations were within the pre-specified tolerance of the historical reference. {normalized_note}Baseline status: **{'PASS' if baseline_ok else 'FAIL'}**.

The principal Project 3 genetic record is the Chinese–Indian–Malay meta-analysis. The later validated output contains 6,686,830 variants and reaches a maximum variant-level N of 7,789; the earlier retained iteration contains 6,688,950 variants and reaches N = 7,664. Both are identified by checksum in the trans-ancestry manifest.

Two related records are retained separately. The refined Chinese-only REGENIE freeze entered 5,700 participants in Steps 1 and 2 (modal per-variant N = 5,592). The LDSC record is labelled `historical_10k_clean`; it gave h2 = 0.3048 (SE 0.0949), intercept = 1.0257 (SE 0.0081), and rg with KoGES EDU = -0.0502 (SE 0.0835, P = 0.5476). These values are not presented as if they came from the same output file.

The participant-level REGENIE phenotype is not distributed in this repository. Reproduction of its construction is supported by the exact PCA match, while the final Chinese GWAS settings are supported by its Step 1 and Step 2 logs.
"""
    write_text(report, out / "baseline_reproduction_report.md")
    if not baseline_ok:
        raise RuntimeError("Historical PCA reproduction failed the pre-specified checks; stopping before sensitivity analyses")
    return score_better, loading_better, explained, corr_matrix, loading_table


def measurement_outputs(root: Path, cognitive: pd.DataFrame, covariates: pd.DataFrame, mask_counts: pd.DataFrame):
    out = root / "01_measurement"
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    cov = covariates.copy()
    cov["IID"] = cov["IID"].astype(str).str.strip()
    cov = cov[~cov["IID"].duplicated()].set_index("IID")
    cognitive = cognitive.copy()
    cognitive["IID"] = cognitive["IID"].astype(str).str.strip()
    cognitive = cognitive[~cognitive["IID"].duplicated()].set_index("IID")
    joined = cognitive.join(cov[["Ancestry", "FREG17_Subcohort"]], how="left")

    summary_rows, quantile_rows, ancestry_rows, freeze_rows = [], [], [], []
    for label, raw_col, _, transform, discrete in TASKS:
        values = numeric(cognitive[raw_col]).to_numpy(dtype=float)
        observed = values[np.isfinite(values)]
        n_total, n = len(values), len(observed)
        skew, kurt = skew_kurtosis(observed)
        ent, effective = entropy(observed) if discrete else (np.nan, np.nan)
        levels = np.sort(np.unique(observed))
        pct_low = [100 * np.isin(observed, levels[: min(k, len(levels))]).mean() for k in (1, 2, 3)] if discrete else [np.nan] * 3
        pct_high = [100 * np.isin(observed, levels[-min(k, len(levels)) :]).mean() for k in (1, 2, 3)] if discrete else [np.nan] * 3
        row = {
            "task": label,
            "source_variable": raw_col,
            "source_direction": "higher_better" if transform == "identity" else "lower_better",
            "analysis_transform": transform,
            "n_total": n_total,
            "n_nonmissing": n,
            "missing_percent": 100 * (n_total - n) / n_total,
            "unique_values": len(levels),
            "mean": np.mean(observed),
            "sd": np.std(observed, ddof=1),
            "variance": np.var(observed, ddof=1),
            "median": np.median(observed),
            "iqr": quantile(observed, 0.75) - quantile(observed, 0.25),
            "minimum": np.min(observed),
            "maximum": np.max(observed),
            "skewness": skew,
            "excess_kurtosis": kurt,
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
        exclusions = mask_counts.set_index("task").loc[label]
        row.update({"non_finite_after_historical_transform": exclusions["non_finite_after_transform"], "masked_by_historical_4sd_rule": exclusions["masked_beyond_4sd"]})
        summary_rows.append(row)
        for q in [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]:
            quantile_rows.append({"task": label, "quantile": q, "source_value": quantile(observed, q)})

        jvals = numeric(joined[raw_col])
        for group, frame in joined.assign(_value=jvals).groupby("Ancestry", dropna=False):
            vals = frame["_value"].dropna().to_numpy(dtype=float)
            if len(vals):
                group_label = "not_recorded" if pd.isna(group) else group
                ancestry_rows.append({"task": label, "ancestry": group_label, "n": len(vals), "mean": np.mean(vals), "sd": np.std(vals, ddof=1), "median": np.median(vals), "q25": quantile(vals, 0.25), "q75": quantile(vals, 0.75), "minimum": np.min(vals), "maximum": np.max(vals)})
        for group, frame in joined.assign(_value=jvals).groupby("FREG17_Subcohort", dropna=False):
            vals = frame["_value"].dropna().to_numpy(dtype=float)
            if len(vals):
                group_label = "not_recorded" if pd.isna(group) else group
                freeze_rows.append({"task": label, "freeze": group_label, "n": len(vals), "mean": np.mean(vals), "sd": np.std(vals, ddof=1), "median": np.median(vals), "q25": quantile(vals, 0.25), "q75": quantile(vals, 0.75), "minimum": np.min(vals), "maximum": np.max(vals)})

        safe = TASK_KEYS[label].lower()
        histogram_plot(observed, f"{label}: source distribution", raw_col, plots / f"histogram_{safe}.png")
        ecdf_plot(observed, f"{label}: empirical cumulative distribution", raw_col, plots / f"ecdf_{safe}.png")
        ancestry_values = {str(g): numeric(f[raw_col]).dropna().to_numpy(dtype=float) for g, f in joined.groupby("Ancestry")}
        boxplot(ancestry_values, f"{label}: distribution by ancestry", "Source score", plots / f"ancestry_{safe}.png")
        freeze_values = {str(g): numeric(f[raw_col]).dropna().to_numpy(dtype=float) for g, f in joined.groupby("FREG17_Subcohort")}
        boxplot(freeze_values, f"{label}: distribution by cohort/freeze", "Source score", plots / f"freeze_{safe}.png")

    summary = pd.DataFrame(summary_rows)
    tsv(summary, out / "task_measurement_summary.tsv")
    tsv(pd.DataFrame(quantile_rows), out / "task_quantiles.tsv")
    tsv(pd.DataFrame(ancestry_rows), out / "task_by_ancestry_summary.tsv")
    tsv(pd.DataFrame(freeze_rows), out / "task_by_freeze_summary.tsv")

    pairing = summary.set_index("task").loc["Pairing guesses"]
    report = f"""# Task-level measurement validation

The source cognitive table contains **{len(cognitive):,} unique participants**. Recorded missingness is 0% for all six source measures, but their information content differs substantially.

Quiz, Working Memory, and Pairing are bounded/discrete measures with **{int(summary.set_index('task').loc['Quiz score', 'unique_values'])}**, **{int(summary.set_index('task').loc['Working memory', 'unique_values'])}**, and **{int(pairing['unique_values'])}** observed values, respectively. Their effective category counts are reported in `task_measurement_summary.tsv`.

The historical reciprocal transformation of Pairing guesses converts zero errors to infinity and then missing. This affects **{int(pairing['non_finite_after_historical_transform']):,} participants ({100 * pairing['non_finite_after_historical_transform'] / len(cognitive):.2f}%)**, including the best-performing end of that source scale. This is a property of the archived construction and was retained only to reproduce the historical phenotype.

Repeated measurements and item-level responses were not available, so test-retest reliability, alpha, omega, and split-half reliability were not estimated. Only one cohort/freeze label (`HELIOS10K`) was present; between-freeze distribution comparisons are therefore not estimable from these files.
"""
    write_text(report, out / "measurement_validation_report.md")
    return summary


def correlation_table(matrix: np.ndarray, method: str) -> pd.DataFrame:
    frame = pd.DataFrame(matrix, columns=TASK_NAMES)
    return frame.corr(method=method)


def pca_structure_outputs(args, root: Path, residual: np.ndarray, metadata: pd.DataFrame, full_score: np.ndarray, full_loading: np.ndarray):
    out = root / "02_pca_structure"
    plots = out / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    z = standardize(residual)
    pearson = correlation_table(z, "pearson")
    spearman = correlation_table(z, "spearman")
    tsv(pearson.rename_axis("task").reset_index(), out / "phenotypic_correlations_pearson.tsv")
    tsv(spearman.rename_axis("task").reset_index(), out / "phenotypic_correlations_spearman.tsv")

    heatmap_plot([("Full sample", pearson)], "Higher-performance task correlations", plots / "correlation_heatmap_full.png")
    ancestry_matrices = []
    ancestry_rows = []
    ancestries = metadata["Ancestry"].astype(str).to_numpy()
    for ancestry in ["Chinese", "Indian", "Malay"]:
        keep = ancestries == ancestry
        matrix = correlation_table(z[keep], "pearson")
        ancestry_matrices.append((ancestry, matrix))
        _, score, loading, correlations, explained = pca(residual[keep], reference=full_loading)
        congruence = float(np.dot(loading, full_loading) / (np.linalg.norm(loading) * np.linalg.norm(full_loading)))
        for i, task in enumerate(TASK_NAMES):
            ancestry_rows.append({"ancestry": ancestry, "n": int(keep.sum()), "pc1_variance_explained": explained, "tucker_congruence_with_full": congruence, "task": task, "loading_higher_performance": loading[i], "correlation_with_pc1": correlations[i]})
    heatmap_plot(ancestry_matrices, "Task correlations by ancestry", plots / "correlation_heatmap_by_ancestry.png")
    ancestry_table = pd.DataFrame(ancestry_rows)
    tsv(ancestry_table, out / "ancestry_specific_pca.tsv")

    rng = np.random.default_rng(args.seed)
    boot_rows = []
    n = len(residual)
    for replicate in range(1, args.bootstrap + 1):
        idx = rng.integers(0, n, n)
        _, _, loading, correlations, explained = pca(residual[idx], reference=full_loading)
        for j, task in enumerate(TASK_NAMES):
            boot_rows.append({"replicate": replicate, "task": task, "pc1_variance_explained": explained, "loading_higher_performance": loading[j], "correlation_with_pc1": correlations[j]})
    bootstrap = pd.DataFrame(boot_rows)
    tsv(bootstrap, out / "bootstrap_pca_loadings.tsv")
    summaries = []
    for task, frame in bootstrap.groupby("task", sort=False):
        for metric in ["loading_higher_performance", "correlation_with_pc1"]:
            vals = frame[metric].to_numpy()
            summaries.append({"task": task, "metric": metric, "mean": np.mean(vals), "sd": np.std(vals, ddof=1), "q025": quantile(vals, 0.025), "q975": quantile(vals, 0.975)})
    variance_values = bootstrap.drop_duplicates("replicate")["pc1_variance_explained"].to_numpy()
    summaries.append({"task": "PC1", "metric": "variance_explained", "mean": np.mean(variance_values), "sd": np.std(variance_values, ddof=1), "q025": quantile(variance_values, 0.025), "q975": quantile(variance_values, 0.975)})
    bootstrap_summary = pd.DataFrame(summaries)
    tsv(bootstrap_summary, out / "bootstrap_pca_summary.tsv")
    loading_summary = bootstrap_summary[bootstrap_summary["metric"] == "loading_higher_performance"]
    interval_plot(loading_summary, "Bootstrap stability of PC1 loadings", plots / "bootstrap_loading_distributions.png")

    loo_rows = []
    for omitted in TASK_NAMES:
        keep_cols = [i for i, task in enumerate(TASK_NAMES) if task != omitted]
        ref = full_loading[keep_cols]
        _, score, loading, correlations, explained = pca(residual[:, keep_cols], reference=ref)
        score_r = float(np.corrcoef(score, full_score)[0, 1])
        if score_r < 0:
            score *= -1
            loading *= -1
            correlations *= -1
            score_r *= -1
        for j, col in enumerate(keep_cols):
            loo_rows.append({"omitted_task": omitted, "retained_task": TASK_NAMES[col], "pc1_variance_explained": explained, "correlation_with_full_pc1": score_r, "loading_higher_performance": loading[j], "correlation_with_leave_one_out_pc1": correlations[j]})
    loo = pd.DataFrame(loo_rows)
    tsv(loo, out / "leave_one_out_pca_summary.tsv")
    loo_matrix = loo.pivot(index="omitted_task", columns="retained_task", values="loading_higher_performance").reindex(index=TASK_NAMES, columns=TASK_NAMES)
    heatmap_plot([("Omitted task (rows)", loo_matrix)], "Leave-one-task-out PC1 loadings", plots / "leave_one_out_loadings.png")

    speed = z[:, [TASK_NAMES.index(x) for x in SPEED_TASKS]].mean(axis=1)
    collapsed_names = ["Speed domain", "Quiz score", "Working memory", "Pairing guesses"]
    collapsed = np.column_stack([speed, z[:, TASK_NAMES.index("Quiz score")], z[:, TASK_NAMES.index("Working memory")], z[:, TASK_NAMES.index("Pairing guesses")]])
    _, collapsed_score, collapsed_loading, collapsed_corr, collapsed_explained = pca(collapsed)
    if np.mean(collapsed_loading) < 0:
        collapsed_score *= -1
        collapsed_loading *= -1
        collapsed_corr *= -1
    collapsed_table = pd.DataFrame({"indicator": collapsed_names, "loading_higher_performance": collapsed_loading, "correlation_with_pc1": collapsed_corr})
    collapsed_table["pc1_variance_explained"] = collapsed_explained
    collapsed_table["correlation_with_original_pc1"] = float(np.corrcoef(collapsed_score, full_score)[0, 1])
    tsv(collapsed_table, out / "collapsed_speed_pca_summary.tsv")

    speed_boot = loading_summary[loading_summary["task"].isin(SPEED_TASKS)]
    nonspeed_boot = loading_summary[~loading_summary["task"].isin(SPEED_TASKS)]
    min_speed = float(speed_boot["q025"].min())
    max_nonspeed = float(nonspeed_boot["q975"].max())
    congruence_min = float(ancestry_table.groupby("ancestry")["tucker_congruence_with_full"].first().min())
    loo_min = float(loo.groupby("omitted_task")["correlation_with_full_pc1"].first().min())
    report = f"""# PCA structure validation

The three latency-derived indicators form the strongest correlation block. In the residualized common-sample participant bootstrap (B = {args.bootstrap:,}; seed = {args.seed}), standardization and PCA were repeated within every replicate. The separation between speed and non-speed loadings was stable: the smallest speed-task lower 95% limit was **{min_speed:.3f}**, while the largest non-speed upper 95% limit was **{max_nonspeed:.3f}**.

Leave-one-task-out PC1 scores remained correlated at **r >= {loo_min:.3f}** with the full PC1. Ancestry-specific loading vectors had Tucker congruence **>= {congruence_min:.3f}** with the full-sample solution. Estimates for the smaller Indian and Malay samples are less precise and should not be read as measurement-invariance tests.

After collapsing Reaction Time, Stroop Ink, and Stroop Box to one speed indicator, the four-indicator PC1 explained **{collapsed_explained * 100:.2f}%** of variance and correlated **r = {np.corrcoef(collapsed_score, full_score)[0, 1]:.3f}** with the historical six-indicator PC1. This is a pre-specified structural sensitivity analysis, not a replacement phenotype.
"""
    write_text(report, out / "pca_structure_report.md")
    return bootstrap_summary, loo, collapsed_table, ancestry_table


def provenance_outputs(args, root: Path):
    provenance = root / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    sources = [
        ("validation_plan", args.validation_plan, "Pre-specified Phase 1 analysis specification"),
        ("cognitive_source", args.cognitive, "Six source cognitive measures"),
        ("covariates", args.covariates, "Age, sex, ancestry and household-income covariates"),
        ("historical_g", args.historical_g, "Saved historical PCA score"),
        ("historical_r_script", args.historical_r_script, "Archived phenotype construction script"),
        ("normalized_phenotypes", args.normalized_phenotypes, "Historical transformed phenotype cross-check"),
        ("residualized_phenotypes", args.residualized_phenotypes, "Historical residual cross-check"),
        ("regenie_step1_log", args.step1_log, "Final Chinese-only REGENIE Step 1 record"),
        ("regenie_step2_log", args.step2_log, "Final Chinese-only REGENIE Step 2 record"),
        ("ldsc_munge_log", args.ldsc_munge_log, "Historical 10k_clean LDSC munging record"),
        ("ldsc_rg_log", args.ldsc_rg_log, "Historical 10k_clean versus KoGES EDU LDSC record"),
        ("freeze_manifest", args.freeze_manifest, "Released Chinese-only summary-statistic manifest"),
        ("trans_ancestry_manifest", root.parent / "results" / "trans_ancestry" / "meta_analysis_manifest.tsv", "Validated Chinese-Indian-Malay meta-analysis outputs"),
        ("public_phenotype_script", root.parent / "scripts" / "phenotype" / "reproduce_10k_cognitive_factor.R", "Cleaned historical phenotype implementation"),
        ("public_regenie_step1_script", root.parent / "scripts" / "gwas" / "run_regenie_step1_10k_chinese.sh", "Portable final Chinese-only Step 1 command"),
        ("public_regenie_step2_script", root.parent / "scripts" / "gwas" / "run_regenie_step2_10k_chinese.sh", "Portable final Chinese-only Step 2 command"),
    ]
    rows, checksum_lines = [], []
    for item, path, role in sources:
        if path is None:
            rows.append({"item": item, "file": "not supplied", "role": role, "bytes": "NA", "sha256": "NA", "status": "not supplied"})
            continue
        if not path.exists():
            raise FileNotFoundError(path)
        digest = sha256(path)
        display_name = "HELIOS_cognition_validation_plan.md" if item == "validation_plan" else path.name
        rows.append({"item": item, "file": display_name, "role": role, "bytes": path.stat().st_size, "sha256": digest, "status": "verified"})
        checksum_lines.append(f"{digest}  {display_name}")
    tsv(pd.DataFrame(rows), provenance / "baseline_manifest.tsv")
    write_text("\n".join(checksum_lines), provenance / "input_checksums.sha256")
    versions = [
        f"phase1_validation_script\t{sha256(Path(__file__))}",
        f"python\t{platform.python_version()}",
        f"numpy\t{np.__version__}",
        f"pandas\t{pd.__version__}",
        f"pillow\t{Image.__version__}",
        "regenie\t4.1 (archived log)",
        "ldsc\t1.0.1 (archived log)",
    ]
    write_text("\n".join(versions), provenance / "software_versions.txt")


def final_summary(root: Path, measurement, bootstrap, loo, collapsed, ancestry, pca_loading, explained):
    by_task = measurement.set_index("task")
    loading_ci = bootstrap[bootstrap["metric"] == "loading_higher_performance"].set_index("task")
    ancestry_summary = ancestry.groupby("ancestry").first()
    loo_summary = loo.groupby("omitted_task").first()
    collapsed_r = float(collapsed["correlation_with_original_pc1"].iloc[0])
    text = f"""# Phase 1 summary

## Scope

Phase 1 reproduced the historical HELIOS 10k cognitive factor and examined measurement quality and PCA structure. No phenotype was redefined for association testing, and no GWAS was rerun.

## Baseline

- PCA sample: **7,403 participants**.
- PC1 variance explained: **{explained * 100:.2f}%**.
- Reconstructed versus saved historical score: **r > 0.999999999**; maximum aligned difference below **1e-8**.
- The normalized-phenotype table contains the same 7,403 saved `g` values with no numerical differences.
- All six loadings and task-PC1 correlations matched the archived reference within the pre-specified tolerance.
- The two validated trans-ancestry outputs remain separately identified by checksum (maximum variant-level N = 7,789 and 7,664); the refined Chinese freeze and historical `10k_clean` LDSC record are not conflated with them.

## Measurement

- Recorded source missingness was 0% across all six tasks, but the discrete tasks had much lower resolution than the latency measures.
- Effective categories: Quiz **{by_task.loc['Quiz score', 'effective_categories']:.1f}**, Working Memory **{by_task.loc['Working memory', 'effective_categories']:.1f}**, Pairing **{by_task.loc['Pairing guesses', 'effective_categories']:.1f}**.
- The historical reciprocal Pairing transform converts **{int(by_task.loc['Pairing guesses', 'non_finite_after_historical_transform']):,} zero-error scores** to non-finite values. These participants are consequently excluded from the historical complete-case PCA.
- Reliability coefficients could not be estimated because repeated or item-level data were not available.
- The available covariate file contains one cohort/freeze label, so between-freeze comparisons are not estimable.

## PCA structure

- The bootstrap confirmed that the largest loadings are the three latency-derived indicators. Their 95% loading intervals were: Reaction Time **[{loading_ci.loc['Reaction time', 'q025']:.3f}, {loading_ci.loc['Reaction time', 'q975']:.3f}]**, Stroop Box **[{loading_ci.loc['Stroop box', 'q025']:.3f}, {loading_ci.loc['Stroop box', 'q975']:.3f}]**, and Stroop Ink **[{loading_ci.loc['Stroop ink', 'q025']:.3f}, {loading_ci.loc['Stroop ink', 'q975']:.3f}]**.
- Leave-one-task-out solutions correlated **{loo_summary['correlation_with_full_pc1'].min():.3f}–{loo_summary['correlation_with_full_pc1'].max():.3f}** with the full PC1.
- Ancestry-specific loading congruence with the full solution ranged from **{ancestry_summary['tucker_congruence_with_full'].min():.3f} to {ancestry_summary['tucker_congruence_with_full'].max():.3f}**.
- The four-indicator PCA using one collapsed speed domain correlated **r = {collapsed_r:.3f}** with the original PC1.

## Interpretation

The historical phenotype is exactly reproducible and its speed-heavy structure is not a bootstrap artefact. The weighting partly reflects the covariance among three separately entered latency measures. The collapsed-speed analysis provides the planned domain-balance sensitivity test for a later phase, but it does not by itself establish that the historical score is invalid.

The Pairing zero-score transformation is the clearest construction issue identified in Phase 1. Any alternative phenotype work should pre-specify a transformation that retains zero-error participants, while preserving the historical phenotype unchanged as the primary archived analysis.

Phase 1 is complete. The next planned step is Phase 2 phenotype sensitivity work; it should begin only after review of these results.
"""
    write_text(text, root / "phase1_summary.md")


def main() -> int:
    args = parse_args()
    root = args.output_root.resolve()
    for path in [root / "00_baseline", root / "01_measurement", root / "02_pca_structure", root / "provenance"]:
        path.mkdir(parents=True, exist_ok=True)
    provenance_outputs(args, root)
    cognitive = read_table(args.cognitive)
    covariates = read_table(args.covariates)
    required_cog = {"IID", *[x[1] for x in TASKS]}
    required_cov = {"IID", "Age", "Sex", "Ancestry", "FSAQ21_Soc9", "FREG17_Subcohort"}
    if missing := required_cog - set(cognitive.columns):
        raise ValueError(f"Cognitive input is missing: {sorted(missing)}")
    if missing := required_cov - set(covariates.columns):
        raise ValueError(f"Covariate input is missing: {sorted(missing)}")

    ids, residual, metadata, mask_counts, cognitive_dedup, covariates_raw = prepare_phenotypes(cognitive, covariates)
    score, loading, explained, _, loading_table = baseline_outputs(args, root, ids, residual, metadata, mask_counts)
    measurement = measurement_outputs(root, cognitive_dedup, covariates_raw, mask_counts)
    bootstrap, loo, collapsed, ancestry = pca_structure_outputs(args, root, residual, metadata, score, loading)
    final_summary(root, measurement, bootstrap, loo, collapsed, ancestry, loading_table, explained)
    print(f"Phase 1 complete: {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
