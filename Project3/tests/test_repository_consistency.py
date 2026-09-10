#!/usr/bin/env python3
"""Check publication terminology, workflow names, and local documentation links."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


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

    check_local_markdown_links()

    print("PASS: repository terminology and supplementary references are consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
