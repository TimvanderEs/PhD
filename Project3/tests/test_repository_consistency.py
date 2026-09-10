#!/usr/bin/env python3
"""Check publication terminology, workflow names, and local documentation links."""

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def read_tsv(relative_path: str) -> list[dict[str, str]]:
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_local_markdown_links() -> None:
    link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    for markdown_file in ROOT.rglob("*.md"):
        for destination in link_pattern.findall(markdown_file.read_text(encoding="utf-8")):
            destination = destination.strip().strip("<>")
            if destination.startswith(("http://", "https://", "mailto:", "#")):
                continue
            path_text = destination.split("#", 1)[0]
            if not path_text:
                continue
            target = (markdown_file.parent / path_text).resolve()
            assert target.exists(), f"Broken link in {markdown_file}: {destination}"


def main() -> int:
    required_files = (
        "scripts/classify_pleio_locus_direction.R",
        "scripts/run_pleio_locus_classification.sh",
        "scripts/gprofiler/plot_gprofiler_enrichment.R",
        "scripts/pleio/isf/EAS_MDD_EA_CF.isf",
        "scripts/pleio/isf/EAS_SCZ_EA_CF.isf",
        "scripts/lava/postprocess_lava.py",
        "scripts/sumstats_qc/README.md",
        "provenance/sumstatsqc_runs.tsv",
        "inputs/smr_heidi/run_manifest.tsv",
        "CODE_AVAILABILITY.md",
    )
    for relative_path in required_files:
        assert (ROOT / relative_path).is_file(), relative_path

    retired_files = (
        "scripts/classify_PLEIO_FUMA_loci_v2.R",
        "scripts/run_PLEIO_reclassification_audit.sh",
        "scripts/gprofiler/plot_gprofiler_s8.R",
        "scripts/pleio/isf/EAS_MDD_EA_G.isf",
        "scripts/pleio/isf/EAS_SCZ_EA_G.isf",
    )
    for relative_path in retired_files:
        assert not (ROOT / relative_path).exists(), relative_path

    classifier = read("scripts/classify_pleio_locus_direction.R")
    for retired_term in (
        "legacy_",
        "corrected_",
        "trait_dominant_or_unclear",
        "mixed_dual",
    ):
        assert retired_term not in classifier, retired_term
    for required_term in (
        '"EA_Z"',
        '"CF_Z"',
        '"Concordant"',
        '"Discordant"',
        '"Dual"',
        '"Mixed"',
        '"Unassigned"',
        "pleio_prioritised",
    ):
        assert required_term in classifier, required_term

    synthetic_header = read("examples/synthetic/raw_merged_PLEIO.tsv").splitlines()[0]
    assert synthetic_header == "SNP\tCHR\tBP\tEA_Z\tCF_Z\tMDD_Z\tSCZ_Z"

    readme = read("README.md")
    assert "Supplementary Table S8" in readme
    assert "Supplementary Figure S7" in readme
    assert "cognitive function" in readme
    assert "Supplementary Figure S8" not in readme

    reader_docs = "\n".join(
        markdown_file.read_text(encoding="utf-8")
        for markdown_file in ROOT.rglob("*.md")
    )
    for retired_text in (
        "cognitive performance",
        "classify_PLEIO_FUMA_loci_v2.R",
        "run_PLEIO_reclassification_audit.sh",
        "plot_gprofiler_s8.R",
        "Supplementary Table 8",
    ):
        assert retired_text not in reader_docs, retired_text

    gprofiler_plot = read("scripts/gprofiler/plot_gprofiler_enrichment.R")
    assert "Supplementary_Figure_S7" in gprofiler_plot
    assert "Supplementary_Figure_S8" not in gprofiler_plot

    for manifest in (
        "inputs/fuma_magma/single_trait_run_manifest.tsv",
        "inputs/fuma_magma/pleio_run_manifest.tsv",
    ):
        assert read(manifest).splitlines()[0].split("\t")[1] == "analysis_label"

    lava_counts = read_tsv("provenance/lava_pair_test_counts.tsv")
    totals: dict[tuple[str, str], int] = {}
    for row in lava_counts:
        key = (row["analysis_set"], row["ancestry"])
        totals[key] = totals.get(key, 0) + int(row["n_bivariate_tests"])
    assert totals[("archived_runner_output", "EAS")] == 726
    assert totals[("exact_threshold_recalculation", "EAS")] == 698
    assert totals[("publication_workbook", "EAS")] == 697
    assert totals[("archived_runner_output", "EUR")] == 3142

    exact_ea_scz = next(
        row for row in lava_counts
        if row["analysis_set"] == "exact_threshold_recalculation"
        and row["ancestry"] == "EAS"
        and row["trait_pair"] == "EA--SCZ"
    )
    workbook_ea_scz = next(
        row for row in lava_counts
        if row["analysis_set"] == "publication_workbook"
        and row["ancestry"] == "EAS"
        and row["trait_pair"] == "EA--SCZ"
    )
    assert int(exact_ea_scz["n_bivariate_tests"]) == 329
    assert int(workbook_ea_scz["n_bivariate_tests"]) == 328

    input_expectations = {
        "scripts/lava/config/EUR_manuscript_target_input.info.tsv": {
            "EA": ("1", "0"), "MDD": ("525197", "3362335"),
            "SCZ": ("53386", "77258"), "CF": ("1", "0"),
        },
        "scripts/lava/config/EAS_manuscript_target_input.info.tsv": {
            "MDD": ("15771", "178777"), "EA": ("1", "0"),
            "SCZ": ("22778", "35362"), "CF": ("1", "0"),
        },
    }
    for path, expected in input_expectations.items():
        observed = {row["phenotype"]: (row["cases"], row["controls"]) for row in read_tsv(path)}
        assert observed == expected

    recovered_eas = {
        row["phenotype"]: (row["cases"], row["controls"])
        for row in read_tsv("scripts/lava/config/EAS_recovered_run_input.info.tsv")
    }
    recovered_eur = {
        row["phenotype"]: (row["cases"], row["controls"])
        for row in read_tsv("scripts/lava/config/EUR_recovered_run_input.info.tsv")
    }
    assert recovered_eas["EA"] == ("15771", "178777")
    assert recovered_eas["SCZ"] == ("13305", "16244")
    assert recovered_eur["MDD"] == ("15771", "178777")
    assert recovered_eur["SCZ"] == ("13305", "16244")

    sumstats_runs = read_tsv("provenance/sumstatsqc_runs.tsv")
    assert len(sumstats_runs) == 8
    assert sum(row["log_status"] == "recovered" for row in sumstats_runs) == 5
    eas_scz = next(row for row in sumstats_runs if row["ancestry"] == "EAS" and row["trait"] == "SCZ")
    assert (eas_scz["log_cases"], eas_scz["manuscript_cases"]) == ("27888", "22778")

    smr_runs = read_tsv("inputs/smr_heidi/run_manifest.tsv")
    assert len(smr_runs) == 7
    assert all(row["portal_job_id"] == "not_recovered" for row in smr_runs)

    fingerprints = read_tsv("provenance/final_output_fingerprints.tsv")
    lava_fingerprints = [row for row in fingerprints if row["workflow"] == "LAVA"]
    assert len(lava_fingerprints) == 4
    assert all(row["provenance_status"] == "archived_requires_rerun" for row in lava_fingerprints)

    ldsc_runner = read("scripts/ldsc/run_genomicsem_ldsc.R")
    for explicit_name in (
        "_genetic_covariance_raw.csv",
        "_genetic_covariance_nearPD.csv",
        "_genetic_correlations_raw.csv",
        "_genetic_correlations_nearPD.csv",
        "_matrix_adjustment_QC.tsv",
    ):
        assert explicit_name in ldsc_runner
    assert 'paste0(out_prefix, "_genetic_covariance.csv")' not in ldsc_runner
    assert 'paste0(out_prefix, "_genetic_correlations.csv")' not in ldsc_runner

    code_availability = read("CODE_AVAILABILITY.md")
    assert "https://github.com/TimvanderEs/PhD/tree/main/Project3" in code_availability
    assert "exact reviewed version is identified by its Git commit" in code_availability

    check_local_markdown_links()

    print("PASS: repository terminology and supplementary references are consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
