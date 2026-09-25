#!/usr/bin/env python3
"""Quantify genome-wide local covariance cancellation from SUPERGNOVA output."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, norm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--orientation-multiplier", type=float, default=-1.0)
    parser.add_argument("--orthant-replicates", type=int, default=100_000)
    parser.add_argument("--interval-replicates", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=20260921)
    return parser.parse_args()


def cancellation(values: np.ndarray) -> float:
    denominator = np.abs(values).sum()
    if denominator == 0:
        return math.nan
    return 1.0 - abs(values.sum()) / denominator


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    n = len(p_values)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def profile_normal_random_effects(y: np.ndarray, variance: np.ndarray) -> dict[str, np.ndarray | float]:
    upper = max(float(np.var(y) * 100.0), float(np.median(variance) * 100.0), 1e-16)

    def fit_at(tau2: float) -> tuple[float, float]:
        total_variance = variance + tau2
        weights = 1.0 / total_variance
        mu = float(np.sum(weights * y) / np.sum(weights))
        nll = 0.5 * float(np.sum(np.log(total_variance) + np.square(y - mu) / total_variance))
        return nll, mu

    result = minimize_scalar(lambda value: fit_at(value)[0], bounds=(0.0, upper), method="bounded")
    candidates = [(0.0, *fit_at(0.0)), (float(result.x), *fit_at(float(result.x)))]
    tau2, nll, mu = min(candidates, key=lambda item: item[1])
    null_nll, null_mu = fit_at(0.0)
    shrinkage = tau2 / (tau2 + variance) if tau2 > 0 else np.zeros_like(variance)
    posterior_mean = mu + shrinkage * (y - mu)
    posterior_variance = tau2 * variance / (tau2 + variance) if tau2 > 0 else np.zeros_like(variance)
    return {
        "mu": mu,
        "tau2": tau2,
        "nll": nll,
        "null_nll": null_nll,
        "null_mu": null_mu,
        "posterior_mean": posterior_mean,
        "posterior_variance": posterior_variance,
    }


def orthant_bootstrap(y: np.ndarray, variance: np.ndarray, replicates: int, rng: np.random.Generator) -> dict[str, float]:
    z = y / np.sqrt(variance)
    observed_against_nonnegative = float(np.square(np.minimum(z, 0.0)).sum())
    observed_against_nonpositive = float(np.square(np.maximum(z, 0.0)).sum())
    exceed_nonnegative = 0
    exceed_nonpositive = 0
    remaining = replicates
    while remaining:
        current = min(500, remaining)
        simulated = rng.standard_normal((current, len(z)))
        against_nonnegative = np.square(np.minimum(simulated, 0.0)).sum(axis=1)
        against_nonpositive = np.square(np.maximum(simulated, 0.0)).sum(axis=1)
        exceed_nonnegative += int((against_nonnegative >= observed_against_nonnegative).sum())
        exceed_nonpositive += int((against_nonpositive >= observed_against_nonpositive).sum())
        remaining -= current
    p_nonnegative = (exceed_nonnegative + 1.0) / (replicates + 1.0)
    p_nonpositive = (exceed_nonpositive + 1.0) / (replicates + 1.0)
    return {
        "statistic_against_all_nonnegative": observed_against_nonnegative,
        "p_against_all_nonnegative": p_nonnegative,
        "statistic_against_all_nonpositive": observed_against_nonpositive,
        "p_against_all_nonpositive": p_nonpositive,
        "intersection_union_p": max(p_nonnegative, p_nonpositive),
        "replicates": replicates,
    }


def bootstrap_cancellation_interval(
    y: np.ndarray,
    variance: np.ndarray,
    fitted_theta: np.ndarray,
    replicates: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    estimates = np.empty(replicates, dtype=float)
    standard_error = np.sqrt(variance)
    for index in range(replicates):
        simulated = fitted_theta + standard_error * rng.standard_normal(len(y))
        fitted = profile_normal_random_effects(simulated, variance)
        estimates[index] = cancellation(np.asarray(fitted["posterior_mean"]))
    return {
        "estimate": cancellation(fitted_theta),
        "lower_95": float(np.quantile(estimates, 0.025)),
        "upper_95": float(np.quantile(estimates, 0.975)),
        "median": float(np.median(estimates)),
        "replicates": replicates,
    }


def summarize(values: np.ndarray, variance: np.ndarray) -> dict[str, float | int]:
    positive = values[values > 0]
    negative = values[values < 0]
    signed = float(values.sum())
    absolute = float(np.abs(values).sum())
    return {
        "n_regions": int(len(values)),
        "n_positive": int((values > 0).sum()),
        "n_negative": int((values < 0).sum()),
        "sum_positive": float(positive.sum()),
        "sum_negative": float(negative.sum()),
        "sum_negative_absolute": float(np.abs(negative).sum()),
        "sum_signed": signed,
        "sum_absolute": absolute,
        "sum_signed_se_independence": float(np.sqrt(variance.sum())),
        "cancellation_fraction": cancellation(values),
    }


def plot_covariance(regions: pd.DataFrame, output: Path) -> None:
    chromosome_max = regions.groupby("chr")["end"].max().sort_index()
    offsets: dict[int, float] = {}
    running = 0.0
    for chromosome, maximum in chromosome_max.items():
        offsets[int(chromosome)] = running
        running += float(maximum)
    x = regions.apply(lambda row: offsets[int(row["chr"])] + row["start"], axis=1)
    colors = np.where(regions["rho_oriented"] >= 0, "#2166ac", "#b2182b")

    fig, axis = plt.subplots(figsize=(12, 4.8))
    axis.axhline(0, color="#333333", linewidth=0.8)
    axis.scatter(x, regions["rho_oriented"], c=colors, s=10, alpha=0.65, linewidths=0)
    for chromosome in range(1, 23):
        if chromosome in offsets:
            axis.axvline(offsets[chromosome], color="#e6e6e6", linewidth=0.5, zorder=0)
    ticks = [offsets[c] + chromosome_max.loc[c] / 2 for c in chromosome_max.index]
    axis.set_xticks(ticks)
    axis.set_xticklabels([str(int(c)) for c in chromosome_max.index], fontsize=8)
    axis.set_xlabel("Chromosome")
    axis.set_ylabel("Local genetic covariance\n(higher HELIOS g = better performance)")
    axis.set_title("HELIOS general cognitive ability and KoGES educational attainment")
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def plot_cumulative(regions: pd.DataFrame, output: Path) -> None:
    ordered = regions.sort_values(["chr", "start"]).reset_index(drop=True)
    values = ordered["rho_oriented"].to_numpy()
    cumulative_positive = np.cumsum(np.maximum(values, 0.0))
    cumulative_negative = np.cumsum(np.minimum(values, 0.0))
    cumulative_net = np.cumsum(values)

    fig, axis = plt.subplots(figsize=(11, 4.8))
    axis.plot(cumulative_positive, color="#2166ac", label="Positive covariance", linewidth=1.5)
    axis.plot(cumulative_negative, color="#b2182b", label="Negative covariance", linewidth=1.5)
    axis.plot(cumulative_net, color="#222222", label="Net covariance", linewidth=2.0)
    axis.axhline(0, color="#777777", linewidth=0.7)
    chromosome_ends = ordered.groupby("chr").size().cumsum()
    chromosome_starts = chromosome_ends.shift(fill_value=0)
    ticks = ((chromosome_starts + chromosome_ends) / 2).to_numpy()
    axis.set_xticks(ticks)
    axis.set_xticklabels([str(int(c)) for c in chromosome_ends.index], fontsize=8)
    axis.set_xlabel("Chromosome")
    axis.set_ylabel("Cumulative local genetic covariance")
    axis.set_title("Opposing regional covariance contributions")
    axis.legend(frameon=False, ncol=3, loc="upper left")
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def plot_decomposition(raw: dict[str, float | int], shrunk: dict[str, float | int], output: Path) -> None:
    labels = ["Raw", "Noise-aware"]
    positive = [raw["sum_positive"], shrunk["sum_positive"]]
    negative = [raw["sum_negative"], shrunk["sum_negative"]]
    net = [raw["sum_signed"], shrunk["sum_signed"]]
    x = np.arange(2)
    width = 0.24
    fig, axis = plt.subplots(figsize=(7.2, 4.8))
    axis.bar(x - width, positive, width, color="#2166ac", label="Positive")
    axis.bar(x, negative, width, color="#b2182b", label="Negative")
    axis.bar(x + width, net, width, color="#333333", label="Net")
    axis.axhline(0, color="#777777", linewidth=0.7)
    axis.set_xticks(x)
    axis.set_xticklabels(labels)
    axis.set_ylabel("Summed local genetic covariance")
    axis.set_title("Regional covariance decomposition")
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    tables = args.outdir / "tables"
    figures = args.outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    regions = pd.read_csv(args.input, sep=r"\s+")
    required = {"chr", "start", "end", "rho", "var", "p", "m"}
    missing = required - set(regions.columns)
    if missing:
        raise ValueError(f"Missing SUPERGNOVA columns: {sorted(missing)}")
    regions = regions.loc[np.isfinite(regions["rho"]) & np.isfinite(regions["var"]) & (regions["var"] > 0)].copy()
    regions["rho_oriented"] = args.orientation_multiplier * regions["rho"]
    regions["se"] = np.sqrt(regions["var"])
    regions["z_oriented"] = regions["rho_oriented"] / regions["se"]
    regions["p_recalculated"] = 2 * norm.sf(np.abs(regions["z_oriented"]))
    regions["p_bh"] = bh_adjust(regions["p_recalculated"].to_numpy())

    y = regions["rho_oriented"].to_numpy()
    variance = regions["var"].to_numpy()
    raw = summarize(y, variance)
    random_effects = profile_normal_random_effects(y, variance)
    posterior_mean = np.asarray(random_effects["posterior_mean"])
    posterior_variance = np.asarray(random_effects["posterior_variance"])
    regions["rho_posterior_mean"] = posterior_mean
    regions["rho_posterior_se"] = np.sqrt(posterior_variance)
    regions["posterior_probability_positive"] = norm.cdf(
        posterior_mean / np.where(posterior_variance > 0, np.sqrt(posterior_variance), np.inf)
    )
    noise_aware = summarize(posterior_mean, posterior_variance)
    inverse_variance_weights = 1.0 / variance
    fixed_effect_mean = float(np.sum(inverse_variance_weights * y) / np.sum(inverse_variance_weights))
    heterogeneity_q = float(np.sum(inverse_variance_weights * np.square(y - fixed_effect_mean)))
    heterogeneity_p = float(chi2.sf(heterogeneity_q, len(y) - 1))
    random_effects_lrt = max(0.0, 2.0 * (float(random_effects["null_nll"]) - float(random_effects["nll"])))
    random_effects_lrt_p = float(0.5 * chi2.sf(random_effects_lrt, 1)) if random_effects_lrt > 0 else 1.0

    rng = np.random.default_rng(args.seed)
    orthant = orthant_bootstrap(y, variance, args.orthant_replicates, rng)
    interval = bootstrap_cancellation_interval(
        y,
        variance,
        posterior_mean,
        args.interval_replicates,
        rng,
    )

    robustness_rows = []
    masks: dict[str, np.ndarray] = {
        "all_regions": np.ones(len(regions), dtype=bool),
        "exclude_mhc_chr6_25_34mb": ~((regions["chr"] == 6) & (regions["start"] <= 34_000_000) & (regions["end"] >= 25_000_000)).to_numpy(),
    }
    absolute_order = np.argsort(np.abs(y))[::-1]
    for count in (1, 5, 10):
        mask = np.ones(len(regions), dtype=bool)
        mask[absolute_order[:count]] = False
        masks[f"exclude_top_{count}_absolute"] = mask
    for label, mask in masks.items():
        subset_y = y[mask]
        subset_variance = variance[mask]
        subset_fit = profile_normal_random_effects(subset_y, subset_variance)
        row = {"analysis": label, **summarize(subset_y, subset_variance)}
        row["noise_aware_cancellation_fraction"] = cancellation(np.asarray(subset_fit["posterior_mean"]))
        robustness_rows.append(row)
    for chromosome in range(1, 23):
        mask = regions["chr"].to_numpy() != chromosome
        subset_y = y[mask]
        subset_variance = variance[mask]
        subset_fit = profile_normal_random_effects(subset_y, subset_variance)
        row = {"analysis": f"leave_chr{chromosome}_out", **summarize(subset_y, subset_variance)}
        row["noise_aware_cancellation_fraction"] = cancellation(np.asarray(subset_fit["posterior_mean"]))
        robustness_rows.append(row)
    robustness = pd.DataFrame(robustness_rows)

    regions.to_csv(tables / "regional_covariance.tsv", sep="\t", index=False)
    robustness.to_csv(tables / "robustness.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"estimator": "raw", **raw},
            {"estimator": "normal_random_effects_posterior_mean", **noise_aware},
        ]
    ).to_csv(tables / "covariance_decomposition.tsv", sep="\t", index=False)
    pd.DataFrame([orthant]).to_csv(tables / "orthant_bootstrap.tsv", sep="\t", index=False)
    pd.DataFrame([interval]).to_csv(tables / "cancellation_interval.tsv", sep="\t", index=False)

    summary = {
        "orientation": "HELIOS g multiplied by -1 so positive values indicate better cognitive performance; EDU positive indicates higher attainment",
        "input_regions_with_estimates": len(regions),
        "nominal_regions": int((regions["p_recalculated"] < 0.05).sum()),
        "fdr_regions": int((regions["p_bh"] < 0.05).sum()),
        "raw": raw,
        "normal_random_effects": {
            "mu": float(random_effects["mu"]),
            "tau2": float(random_effects["tau2"]),
            "likelihood_ratio_statistic_vs_tau2_zero": random_effects_lrt,
            "mixture_chi_square_p_vs_tau2_zero": random_effects_lrt_p,
            "decomposition": noise_aware,
            "regions_posterior_probability_positive_gt_0_95": int((regions["posterior_probability_positive"] > 0.95).sum()),
            "regions_posterior_probability_negative_gt_0_95": int((regions["posterior_probability_positive"] < 0.05).sum()),
        },
        "fixed_effect_heterogeneity": {
            "inverse_variance_weighted_mean": fixed_effect_mean,
            "cochran_q": heterogeneity_q,
            "degrees_of_freedom": len(y) - 1,
            "p": heterogeneity_p,
        },
        "orthant_bootstrap": orthant,
        "cancellation_interval": interval,
        "seed": args.seed,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    plot_covariance(regions, figures / "local_covariance_genome.png")
    plot_cumulative(regions, figures / "cumulative_covariance.png")
    plot_decomposition(raw, noise_aware, figures / "covariance_decomposition.png")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
