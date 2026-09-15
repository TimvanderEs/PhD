#!/usr/bin/env python3
"""
Reanalyse merged webSMR/HEIDI outputs using standard multiple-testing correction
and rebuild reviewer-facing sensitivity and PLEIO/LAVA convergence tables.

Expected inputs:
  --smr-root points to a directory containing EU_webSMR and EAS_SMR_web.
  --supp-tsv optionally points to canonical PLEIO/LAVA integration tables.

No third-party Python packages are required.

Primary significance definition
-------------------------------
For each ancestry-trait analysis:
    p_SMR < 0.05 / N_valid_SMR_tests
    p_HEIDI > 0.01

Here N_valid_SMR_tests is the number of valid rows in the merged 15-resource
webSMR output. Because each row is a probe-by-QTL-resource test, this is a
conservative, reviewer-defensible correction across the full merged analysis.

Additional sensitivity tiers
----------------------------
1. Study-wide Bonferroni: p_SMR < 0.05 / total valid tests across all analyses.
2. Ultra-stringent legacy: p_SMR < 5e-8.
3. HEIDI p > 0.05.
4. At least 10 HEIDI SNPs.
5. SMR-multi also passes the trait-specific Bonferroni threshold.
6. Unique-probe Bonferroni (reported as sensitivity only).

Outputs
-------
00_SMR_reanalysis_QC.tsv
20_SMR_HEIDI_all_results_reanalysed.tsv.gz
21_SMR_threshold_comparison.tsv
22A_SMR_bonf_hits.tsv
22B_SMR_ultra_hits.tsv
22C_SMR_studywide_hits.tsv
22D_SMR_bonf_unique_probe_sensitivity.tsv
23A_SMR_bonf_unique_genes.tsv
23B_SMR_ultra_unique_genes.tsv
23C_SMR_cross_ancestry_bonf.tsv
24A_PLEIO_SMR_threshold_summary.tsv
24B_PLEIO_SMR_threshold_detail.tsv
24C_PLEIO_multilayer_thresholds.tsv
99_SMR_reanalysis_checks.tsv

The PLEIO/LAVA files are optional. If the canonical supplementary TSVs are
present, the integration tables are generated. Otherwise the SMR reanalysis
still completes and the missing integration files are listed in QC.
"""

from __future__ import print_function

import argparse
import csv
import gzip
import hashlib
import math
import os
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path


FILE_SPECS = [
    ("EUR", "EA",  "EU_webSMR/EU_EDU_SMR_merged_raw.zip", "cognitive_educational"),
    ("EUR", "CF",  "EU_webSMR/EU_COG_SMR_merged_raw.zip", "cognitive_educational"),
    ("EUR", "MDD", "EU_webSMR/EU_MDD_SMR_merged_raw.zip", "psychiatric"),
    ("EUR", "SCZ", "EU_webSMR/EU_SCZ_SMR_merged_raw.zip", "psychiatric"),
    ("EAS", "EA",  "EAS_SMR_web/EDU_EAS_eSMR.merged.tsv", "cognitive_educational"),
    ("EAS", "MDD", "EAS_SMR_web/EAS_MDD_eSMR.merged.tsv", "psychiatric"),
    ("EAS", "SCZ", "EAS_SMR_web/SCZ_EAS_eSMR.merged.tsv", "psychiatric"),
]

MODEL_TRAITS = {
    "FOURTRAIT": {"EA", "CF", "MDD", "SCZ"},
    "MDD_3TRAIT": {"EA", "CF", "MDD"},
    "SCZ_3TRAIT": {"EA", "CF", "SCZ"},
}

DISPLAY_MODEL = {
    ("EUR", "FOURTRAIT"): "EUR_FOURTRAIT",
    ("EUR", "MDD_3TRAIT"): "EUR_MDD_EA_CF",
    ("EUR", "SCZ_3TRAIT"): "EUR_SCZ_EA_CF",
    ("EAS", "FOURTRAIT"): "EAS_4TRAIT",
    ("EAS", "MDD_3TRAIT"): "EAS_MDD_EA_CF",
    ("EAS", "SCZ_3TRAIT"): "EAS_SCZ_EA_CF",
}

EXPECTED_PLEIO_TOTALS = {
    ("EUR", "FOURTRAIT"): 998,
    ("EUR", "MDD_3TRAIT"): 887,
    ("EUR", "SCZ_3TRAIT"): 826,
    ("EAS", "FOURTRAIT"): 23,
    ("EAS", "MDD_3TRAIT"): 8,
    ("EAS", "SCZ_3TRAIT"): 27,
}

BASE_OUTPUT_COLUMNS = [
    "ancestry", "trait", "trait_domain", "qtl_name", "gene", "probe_id", "ensg",
    "probe_chr", "probe_bp", "top_snp", "top_snp_chr", "top_snp_bp",
    "A1", "A2", "Freq",
    "b_GWAS", "se_GWAS", "p_GWAS",
    "b_eQTL", "se_eQTL", "p_eQTL",
    "b_SMR", "se_SMR", "p_SMR", "p_SMR_multi", "p_HEIDI", "n_HEIDI_SNPs",
    "source_file", "source_inner_file"
]

FLAG_COLUMNS = [
    "n_valid_tests_trait",
    "bonf_threshold_all_tests",
    "n_unique_probes_trait",
    "bonf_threshold_unique_probes",
    "studywide_bonf_threshold",
    "passes_SMR_bonf_all_tests",
    "passes_SMR_bonf_unique_probes",
    "passes_SMR_studywide",
    "passes_SMR_5e8",
    "passes_HEIDI_001",
    "passes_HEIDI_005",
    "HEIDI_n_ge_10",
    "passes_SMR_multi_bonf",
    "primary_bonf_HEIDI_hit",
    "studywide_HEIDI_hit",
    "ultra_5e8_HEIDI_hit",
    "bonf_HEIDI005_hit",
    "bonf_HEIDI_n10_hit",
    "bonf_SMRmulti_supported_hit",
]


def eprint(*args):
    print(*args, file=sys.stderr)


def clean_text(x):
    if x is None:
        return ""
    return str(x).strip()


def parse_float(x):
    s = clean_text(x)
    if not s or s.lower() in {"nan", "na", "none", ".", "null"}:
        return None
    try:
        v = float(s)
    except Exception:
        return None
    if not math.isfinite(v):
        return None
    return v


def parse_int(x):
    v = parse_float(x)
    if v is None:
        return None
    return int(round(v))


def bool_text(v):
    return "TRUE" if bool(v) else "FALSE"


def norm_ensg(x):
    s = clean_text(x)
    if not s:
        return ""
    m = re.search(r"(ENSG\d+)", s, re.I)
    return m.group(1).upper() if m else s.upper()


def norm_gene(x):
    return clean_text(x).upper()


def norm_ancestry(x):
    s = clean_text(x).upper()
    if s in {"EU", "EUR", "EUROPEAN"}:
        return "EUR"
    if s in {"EA", "EAS", "EAST_ASIAN", "EAST ASIAN"}:
        return "EAS"
    return s


def norm_trait(x):
    s = clean_text(x).upper()
    mapping = {
        "EDU": "EA", "EDUCATION": "EA", "EDUCATIONAL_ATTAINMENT": "EA",
        "G": "CF", "COG": "CF", "COGNITIVE_FUNCTION": "CF", "HEL10K": "CF",
        "DEPRESSION": "MDD", "SCHIZOPHRENIA": "SCZ",
    }
    return mapping.get(s, s)


def norm_model(x):
    s = clean_text(x).upper().replace("-", "_").replace(" ", "_")
    if "FOUR" in s or "4TRAIT" in s or "4_TRAIT" in s:
        return "FOURTRAIT"
    if "MDD" in s and ("3TRAIT" in s or "3_TRAIT" in s or "EA_CF" in s or "EA_G" in s):
        return "MDD_3TRAIT"
    if "SCZ" in s and ("3TRAIT" in s or "3_TRAIT" in s or "EA_CF" in s or "EA_G" in s):
        return "SCZ_3TRAIT"
    return s


def pick(row, names, default=""):
    lower = {str(k).lower(): k for k in row.keys()}
    for name in names:
        if name in row:
            return row.get(name, default)
        k = lower.get(name.lower())
        if k is not None:
            return row.get(k, default)
    return default


def open_tabular(path):
    """
    Return (text_handle, inner_name, closer).
    Caller must call closer().
    """
    path = Path(path)
    if path.suffix.lower() == ".zip":
        zf = zipfile.ZipFile(str(path), "r")
        candidates = [
            n for n in zf.namelist()
            if not n.endswith("/") and n.lower().endswith((".tsv", ".txt", ".csv"))
        ]
        if not candidates:
            zf.close()
            raise RuntimeError("No tabular file found inside {0}".format(path))
        preferred = [n for n in candidates if "merged" in n.lower()]
        inner = preferred[0] if preferred else candidates[0]
        raw = zf.open(inner, "r")
        import io
        txt = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
        def closer():
            try:
                txt.close()
            finally:
                zf.close()
        return txt, inner, closer

    fh = path.open("r", encoding="utf-8", errors="replace", newline="")
    return fh, path.name, fh.close


def sniff_delimiter(handle):
    pos = handle.tell()
    sample = handle.read(8192)
    handle.seek(pos)
    if sample.count("\t") >= sample.count(","):
        return "\t"
    return ","


def iterate_raw(path):
    handle, inner, closer = open_tabular(path)
    try:
        delim = sniff_delimiter(handle)
        reader = csv.DictReader(handle, delimiter=delim)
        for row in reader:
            yield row, inner
    finally:
        closer()


def canonicalise_raw(row, ancestry, trait, domain, source_file, inner):
    gene = pick(row, ["Gene", "gene", "symbol"])
    if not clean_text(gene):
        # In EUR files, pandas-exported gene symbol is usually in the "index" column.
        candidate = pick(row, ["index"])
        if candidate and not str(candidate).upper().startswith("ENSG"):
            gene = candidate

    probe = pick(row, ["probeID", "probe_id", "ProbeID"])
    ensg = pick(row, ["gene_id", "ensg", "ensembl_gene_id"])
    if not clean_text(ensg):
        ensg = probe

    out = {
        "ancestry": ancestry,
        "trait": trait,
        "trait_domain": domain,
        "qtl_name": clean_text(pick(row, ["qtl_name", "QTL", "qtl"])),
        "gene": clean_text(gene),
        "probe_id": clean_text(probe),
        "ensg": norm_ensg(ensg),
        "probe_chr": clean_text(pick(row, ["ProbeChr", "probe_chr"])),
        "probe_bp": clean_text(pick(row, ["Probe_bp", "probe_bp"])),
        "top_snp": clean_text(pick(row, ["topSNP", "top_snp"])),
        "top_snp_chr": clean_text(pick(row, ["topSNP_chr", "top_snp_chr"])),
        "top_snp_bp": clean_text(pick(row, ["topSNP_bp", "top_snp_bp"])),
        "A1": clean_text(pick(row, ["A1", "a1"])),
        "A2": clean_text(pick(row, ["A2", "a2"])),
        "Freq": clean_text(pick(row, ["Freq", "freq", "EAF"])),
        "b_GWAS": clean_text(pick(row, ["b_GWAS", "b_gwas"])),
        "se_GWAS": clean_text(pick(row, ["se_GWAS", "se_gwas"])),
        "p_GWAS": clean_text(pick(row, ["p_GWAS", "p_gwas"])),
        "b_eQTL": clean_text(pick(row, ["b_eQTL", "b_eqtl"])),
        "se_eQTL": clean_text(pick(row, ["se_eQTL", "se_eqtl"])),
        "p_eQTL": clean_text(pick(row, ["p_eQTL", "p_eqtl"])),
        "b_SMR": clean_text(pick(row, ["b_SMR", "b_smr"])),
        "se_SMR": clean_text(pick(row, ["se_SMR", "se_smr"])),
        "p_SMR": clean_text(pick(row, ["p_SMR", "p_smr"])),
        "p_SMR_multi": clean_text(pick(row, ["p_SMR_multi", "p_smr_multi"])),
        "p_HEIDI": clean_text(pick(row, ["p_HEIDI", "p_heidi"])),
        "n_HEIDI_SNPs": clean_text(pick(row, ["nsnp_HEIDI", "n_HEIDI_SNPs", "n_heidi_snps"])),
        "source_file": source_file,
        "source_inner_file": inner,
    }
    return out


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_gzip_tsv_header(path, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = gzip.open(str(path), "wt", encoding="utf-8", newline="")
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    return fh, writer


def md5(path):
    h = hashlib.md5()
    with Path(path).open("rb") as fh:
        while True:
            b = fh.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def first_pass(smr_root):
    stats = {}
    qc = []
    for ancestry, trait, rel, domain in FILE_SPECS:
        path = Path(smr_root) / rel
        if not path.exists():
            raise FileNotFoundError("Missing raw SMR file: {0}".format(path))
        n_rows = 0
        n_valid = 0
        probes = set()
        genes = set()
        qtls = set()
        for raw, inner in iterate_raw(path):
            n_rows += 1
            c = canonicalise_raw(raw, ancestry, trait, domain, path.name, inner)
            p = parse_float(c["p_SMR"])
            if p is not None:
                n_valid += 1
            if c["probe_id"]:
                probes.add(c["probe_id"])
            if c["ensg"]:
                genes.add(c["ensg"])
            elif c["gene"]:
                genes.add(norm_gene(c["gene"]))
            if c["qtl_name"]:
                qtls.add(c["qtl_name"])
        key = (ancestry, trait)
        stats[key] = {
            "ancestry": ancestry,
            "trait": trait,
            "trait_domain": domain,
            "source_path": str(path),
            "source_file": path.name,
            "n_rows": n_rows,
            "n_valid_p_SMR_tests": n_valid,
            "n_unique_probe_ids": len(probes),
            "n_unique_genes_all_rows": len(genes),
            "n_qtl_resources": len(qtls),
            "bonf_threshold_all_tests": (0.05 / n_valid) if n_valid else None,
            "bonf_threshold_unique_probes": (0.05 / len(probes)) if probes else None,
        }
        qc.append({
            "check": "RAW_FILE_PRESENT_AND_READABLE",
            "ancestry": ancestry,
            "trait": trait,
            "status": "PASS",
            "observed": n_rows,
            "expected_or_note": str(path),
        })
    total_valid = sum(x["n_valid_p_SMR_tests"] for x in stats.values())
    studywide = 0.05 / total_valid if total_valid else None
    return stats, studywide, qc


def derive_flags(c, st, studywide):
    p = parse_float(c["p_SMR"])
    pm = parse_float(c["p_SMR_multi"])
    ph = parse_float(c["p_HEIDI"])
    nh = parse_int(c["n_HEIDI_SNPs"])
    t_all = st["bonf_threshold_all_tests"]
    t_probe = st["bonf_threshold_unique_probes"]

    pass_bonf = p is not None and t_all is not None and p < t_all
    pass_probe = p is not None and t_probe is not None and p < t_probe
    pass_study = p is not None and studywide is not None and p < studywide
    pass_ultra = p is not None and p < 5e-8
    h001 = ph is not None and ph > 0.01
    h005 = ph is not None and ph >= 0.05
    hn10 = nh is not None and nh >= 10
    pm_bonf = pm is not None and t_all is not None and pm < t_all

    flags = {
        "n_valid_tests_trait": st["n_valid_p_SMR_tests"],
        "bonf_threshold_all_tests": t_all,
        "n_unique_probes_trait": st["n_unique_probe_ids"],
        "bonf_threshold_unique_probes": t_probe,
        "studywide_bonf_threshold": studywide,
        "passes_SMR_bonf_all_tests": bool_text(pass_bonf),
        "passes_SMR_bonf_unique_probes": bool_text(pass_probe),
        "passes_SMR_studywide": bool_text(pass_study),
        "passes_SMR_5e8": bool_text(pass_ultra),
        "passes_HEIDI_001": bool_text(h001),
        "passes_HEIDI_005": bool_text(h005),
        "HEIDI_n_ge_10": bool_text(hn10),
        "passes_SMR_multi_bonf": bool_text(pm_bonf),
        "primary_bonf_HEIDI_hit": bool_text(pass_bonf and h001),
        "studywide_HEIDI_hit": bool_text(pass_study and h001),
        "ultra_5e8_HEIDI_hit": bool_text(pass_ultra and h001),
        "bonf_HEIDI005_hit": bool_text(pass_bonf and h005),
        "bonf_HEIDI_n10_hit": bool_text(pass_bonf and h001 and hn10),
        "bonf_SMRmulti_supported_hit": bool_text(pass_bonf and h001 and pm_bonf),
    }
    return flags


def second_pass(smr_root, output_dir, stats, studywide):
    all_path = Path(output_dir) / "20_SMR_HEIDI_all_results_reanalysed.tsv.gz"
    all_fh, all_writer = write_gzip_tsv_header(all_path, BASE_OUTPUT_COLUMNS + FLAG_COLUMNS)

    tier_rows = {
        "bonf": [],
        "ultra": [],
        "studywide": [],
        "probe": [],
    }
    try:
        for ancestry, trait, rel, domain in FILE_SPECS:
            path = Path(smr_root) / rel
            st = stats[(ancestry, trait)]
            for raw, inner in iterate_raw(path):
                c = canonicalise_raw(raw, ancestry, trait, domain, path.name, inner)
                flags = derive_flags(c, st, studywide)
                row = dict(c)
                row.update(flags)
                all_writer.writerow(row)
                if flags["primary_bonf_HEIDI_hit"] == "TRUE":
                    tier_rows["bonf"].append(row)
                if flags["ultra_5e8_HEIDI_hit"] == "TRUE":
                    tier_rows["ultra"].append(row)
                if flags["studywide_HEIDI_hit"] == "TRUE":
                    tier_rows["studywide"].append(row)
                if flags["passes_SMR_bonf_unique_probes"] == "TRUE" and flags["passes_HEIDI_001"] == "TRUE":
                    tier_rows["probe"].append(row)
    finally:
        all_fh.close()

    fields = BASE_OUTPUT_COLUMNS + FLAG_COLUMNS
    write_tsv(Path(output_dir) / "22A_SMR_bonf_hits.tsv", tier_rows["bonf"], fields)
    write_tsv(Path(output_dir) / "22B_SMR_ultra_hits.tsv", tier_rows["ultra"], fields)
    write_tsv(Path(output_dir) / "22C_SMR_studywide_hits.tsv", tier_rows["studywide"], fields)
    write_tsv(Path(output_dir) / "22D_SMR_bonf_unique_probe_sensitivity.tsv", tier_rows["probe"], fields)
    return tier_rows


def gene_key(row):
    ensg = norm_ensg(row.get("ensg", ""))
    if ensg:
        return ("ENSG", ensg)
    return ("SYMBOL", norm_gene(row.get("gene", "")))


def unique_gene_rows(rows):
    best = {}
    for row in rows:
        key = (row["ancestry"], row["trait"], gene_key(row))
        p = parse_float(row.get("p_SMR"))
        old = best.get(key)
        if old is None or (p is not None and (parse_float(old.get("p_SMR")) is None or p < parse_float(old.get("p_SMR")))):
            best[key] = row
    output = []
    rank_groups = defaultdict(list)
    for row in best.values():
        rank_groups[(row["ancestry"], row["trait"])].append(row)
    for grp, arr in sorted(rank_groups.items()):
        arr.sort(key=lambda r: parse_float(r.get("p_SMR")) if parse_float(r.get("p_SMR")) is not None else 1.0)
        for rank, row in enumerate(arr, start=1):
            r = dict(row)
            r["unique_gene_rank_within_trait"] = rank
            output.append(r)
    return output


def threshold_summary(stats, studywide, tier_rows):
    result = []
    for key in sorted(stats):
        st = stats[key]
        ancestry, trait = key
        subsets = {}
        for tier, rows in tier_rows.items():
            subsets[tier] = [r for r in rows if r["ancestry"] == ancestry and r["trait"] == trait]
        bonf = subsets["bonf"]
        ultra = subsets["ultra"]
        study = subsets["studywide"]
        probe = subsets["probe"]
        p_smr_bonf_before_heidi = 0
        # Derive from all-results flags later would require reread; use raw second pass counts stored below.
        result.append({
            "ancestry": ancestry,
            "trait": trait,
            "trait_domain": st["trait_domain"],
            "n_rows": st["n_rows"],
            "n_valid_p_SMR_tests": st["n_valid_p_SMR_tests"],
            "n_unique_probe_ids": st["n_unique_probe_ids"],
            "n_unique_genes_all_rows": st["n_unique_genes_all_rows"],
            "n_qtl_resources": st["n_qtl_resources"],
            "bonf_threshold_all_tests": st["bonf_threshold_all_tests"],
            "bonf_threshold_unique_probes": st["bonf_threshold_unique_probes"],
            "studywide_bonf_threshold": studywide,
            "legacy_ultra_threshold": 5e-8,
            "bonf_HEIDI001_rows": len(bonf),
            "bonf_HEIDI001_unique_genes": len({gene_key(r) for r in bonf}),
            "bonf_HEIDI005_rows": sum(1 for r in bonf if r["passes_HEIDI_005"] == "TRUE"),
            "bonf_HEIDI_n10_rows": sum(1 for r in bonf if r["HEIDI_n_ge_10"] == "TRUE"),
            "bonf_SMRmulti_supported_rows": sum(1 for r in bonf if r["passes_SMR_multi_bonf"] == "TRUE"),
            "studywide_HEIDI001_rows": len(study),
            "studywide_HEIDI001_unique_genes": len({gene_key(r) for r in study}),
            "ultra_HEIDI001_rows": len(ultra),
            "ultra_HEIDI001_unique_genes": len({gene_key(r) for r in ultra}),
            "unique_probe_bonf_HEIDI001_rows": len(probe),
            "unique_probe_bonf_HEIDI001_unique_genes": len({gene_key(r) for r in probe}),
            "additional_bonf_rows_vs_ultra": len(bonf) - len(ultra),
            "additional_bonf_unique_genes_vs_ultra":
                len({gene_key(r) for r in bonf}) - len({gene_key(r) for r in ultra}),
        })
    return result


def cross_ancestry_table(unique_rows):
    by_gene = defaultdict(lambda: defaultdict(list))
    for r in unique_rows:
        key = gene_key(r)
        by_gene[key][r["ancestry"]].append(r)
    out = []
    for key, anc in by_gene.items():
        if "EUR" not in anc or "EAS" not in anc:
            continue
        eur_traits = sorted({x["trait"] for x in anc["EUR"]})
        eas_traits = sorted({x["trait"] for x in anc["EAS"]})
        eur_best = min(anc["EUR"], key=lambda x: parse_float(x["p_SMR"]) or 1.0)
        eas_best = min(anc["EAS"], key=lambda x: parse_float(x["p_SMR"]) or 1.0)
        out.append({
            "gene": eur_best["gene"] or eas_best["gene"],
            "ensg": eur_best["ensg"] or eas_best["ensg"],
            "EUR_traits": ";".join(eur_traits),
            "EAS_traits": ";".join(eas_traits),
            "EUR_best_trait": eur_best["trait"],
            "EUR_best_p_SMR": eur_best["p_SMR"],
            "EUR_best_p_HEIDI": eur_best["p_HEIDI"],
            "EUR_best_b_SMR": eur_best["b_SMR"],
            "EAS_best_trait": eas_best["trait"],
            "EAS_best_p_SMR": eas_best["p_SMR"],
            "EAS_best_p_HEIDI": eas_best["p_HEIDI"],
            "EAS_best_b_SMR": eas_best["b_SMR"],
            "same_trait_present": bool_text(bool(set(eur_traits) & set(eas_traits))),
            "effect_direction_concordant_best": bool_text(
                (parse_float(eur_best["b_SMR"]) or 0) * (parse_float(eas_best["b_SMR"]) or 0) > 0
            ),
        })
    out.sort(key=lambda r: min(parse_float(r["EUR_best_p_SMR"]) or 1, parse_float(r["EAS_best_p_SMR"]) or 1))
    return out


def read_tsv_rows(path):
    path = Path(path)
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [dict(r) for r in reader]


def detect_col(rows, candidates):
    if not rows:
        return None
    cols = list(rows[0].keys())
    lower = {c.lower(): c for c in cols}
    for c in candidates:
        if c in cols:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def to_bool(x):
    s = clean_text(x).upper()
    return s in {"TRUE", "T", "YES", "Y", "1"}


def integration_inputs(supp_tsv):
    supp = Path(supp_tsv)
    locus = supp / "12_PLEIO_locus_master.tsv"
    gene = supp / "14_PLEIO_gene_direction_map.tsv"
    lava = supp / "11_LAVA_PLEIO_overlap.tsv"
    return locus, gene, lava


def normalise_locus_master(rows):
    if not rows:
        return []
    c_anc = detect_col(rows, ["ancestry", "ANCESTRY"])
    c_mod = detect_col(rows, ["model", "pleio_model", "MODEL", "analysis_model"])
    c_disp = detect_col(rows, ["display_model", "model_display", "trait", "analysis"])
    c_locus = detect_col(rows, ["pleio_locus_id", "locus_id", "canonical_locus_id"])
    c_fuma = detect_col(rows, ["fuma_locus_id", "GenomicLocus", "genomic_locus", "pleio_fuma_locus"])
    c_chr = detect_col(rows, ["chr", "CHR", "chromosome"])
    c_start = detect_col(rows, ["start", "start_bp", "START"])
    c_end = detect_col(rows, ["end", "end_bp", "END"])
    c_dir = detect_col(rows, ["locus_direction_class", "direction_class", "pleio_direction_class"])
    c_prio = detect_col(rows, ["pleio_prioritised", "novel_vs_single_trait_FUMA", "PLEIO_prioritised"])
    c_status = detect_col(rows, ["prioritisation_status", "single_trait_overlap_status"])
    out = []
    counters = defaultdict(int)
    for r in rows:
        anc = norm_ancestry(r.get(c_anc, ""))
        model_raw = r.get(c_mod, "") if c_mod else (r.get(c_disp, "") if c_disp else "")
        model = norm_model(model_raw)
        if (anc, model) not in EXPECTED_PLEIO_TOTALS:
            continue
        counters[(anc, model)] += 1
        locus_id = clean_text(r.get(c_locus, "")) if c_locus else ""
        if not locus_id:
            locus_id = "{0}_{1}_L{2:04d}".format(anc, model, counters[(anc, model)])
        out.append({
            "ancestry": anc,
            "pleio_model": model,
            "display_model": clean_text(r.get(c_disp, "")) if c_disp else DISPLAY_MODEL.get((anc, model), ""),
            "pleio_locus_id": locus_id,
            "pleio_fuma_locus": clean_text(r.get(c_fuma, "")) if c_fuma else "",
            "chr": clean_text(r.get(c_chr, "")) if c_chr else "",
            "start": clean_text(r.get(c_start, "")) if c_start else "",
            "end": clean_text(r.get(c_end, "")) if c_end else "",
            "pleio_direction_class": clean_text(r.get(c_dir, "")) if c_dir else "",
            "pleio_prioritised": to_bool(r.get(c_prio, "")) if c_prio else False,
            "prioritisation_status": clean_text(r.get(c_status, "")) if c_status else "",
        })
    return out


def normalise_gene_map(rows, loci):
    if not rows:
        return []
    c_anc = detect_col(rows, ["ancestry", "ANCESTRY"])
    c_mod = detect_col(rows, ["model", "pleio_model", "MODEL", "analysis_model"])
    c_disp = detect_col(rows, ["display_model", "model_display", "trait", "analysis"])
    c_locus = detect_col(rows, ["pleio_locus_id", "locus_id", "canonical_locus_id"])
    c_fuma = detect_col(rows, ["fuma_locus_id", "GenomicLocus", "genomic_locus", "pleio_fuma_locus"])
    c_chr = detect_col(rows, ["chr", "CHR", "chromosome"])
    c_start = detect_col(rows, ["start", "start_bp", "START"])
    c_end = detect_col(rows, ["end", "end_bp", "END"])
    c_gene = detect_col(rows, ["gene", "symbol", "gene_symbol", "mapped_gene", "Gene"])
    c_ensg = detect_col(rows, ["ensg", "gene_id", "ensembl_gene_id", "Ensembl_ID"])

    locus_lookup = {}
    for l in loci:
        keys = [
            (l["ancestry"], l["pleio_model"], l["pleio_locus_id"]),
            (l["ancestry"], l["pleio_model"], l["pleio_fuma_locus"]),
            (l["ancestry"], l["pleio_model"], l["chr"], l["start"], l["end"]),
        ]
        for k in keys:
            locus_lookup[k] = l

    out = []
    seen = set()
    for r in rows:
        anc = norm_ancestry(r.get(c_anc, ""))
        model_raw = r.get(c_mod, "") if c_mod else (r.get(c_disp, "") if c_disp else "")
        model = norm_model(model_raw)
        if (anc, model) not in EXPECTED_PLEIO_TOTALS:
            continue
        loc = None
        if c_locus:
            loc = locus_lookup.get((anc, model, clean_text(r.get(c_locus, ""))))
        if loc is None and c_fuma:
            loc = locus_lookup.get((anc, model, clean_text(r.get(c_fuma, ""))))
        if loc is None and c_chr and c_start and c_end:
            loc = locus_lookup.get((
                anc, model, clean_text(r.get(c_chr, "")),
                clean_text(r.get(c_start, "")), clean_text(r.get(c_end, ""))
            ))
        if loc is None:
            continue
        gene = clean_text(r.get(c_gene, "")) if c_gene else ""
        ensg = norm_ensg(r.get(c_ensg, "")) if c_ensg else ""
        if not gene and not ensg:
            continue
        key = (loc["pleio_locus_id"], ensg, norm_gene(gene))
        if key in seen:
            continue
        seen.add(key)
        x = dict(loc)
        x.update({"gene": gene, "ensg": ensg})
        out.append(x)
    return out


def lava_supported_loci(rows, loci):
    if not rows:
        return set()
    c_anc = detect_col(rows, ["ancestry", "pleio_ancestry", "ANCESTRY"])
    c_mod = detect_col(rows, ["model", "pleio_model", "MODEL"])
    c_locus = detect_col(rows, ["pleio_locus_id", "locus_id", "canonical_locus_id"])
    c_fuma = detect_col(rows, ["fuma_locus_id", "GenomicLocus", "pleio_fuma_locus"])
    c_chr = detect_col(rows, ["pleio_chr", "chr", "CHR"])
    c_start = detect_col(rows, ["pleio_start", "start", "start_bp"])
    c_end = detect_col(rows, ["pleio_end", "end", "end_bp"])

    by_id = {(l["ancestry"], l["pleio_model"], l["pleio_locus_id"]): l["pleio_locus_id"] for l in loci}
    by_fuma = {(l["ancestry"], l["pleio_model"], l["pleio_fuma_locus"]): l["pleio_locus_id"] for l in loci if l["pleio_fuma_locus"]}
    by_coord = {(l["ancestry"], l["pleio_model"], l["chr"], l["start"], l["end"]): l["pleio_locus_id"] for l in loci}
    supported = set()
    for r in rows:
        anc = norm_ancestry(r.get(c_anc, "")) if c_anc else ""
        model = norm_model(r.get(c_mod, "")) if c_mod else ""
        loc_id = None
        if c_locus:
            loc_id = by_id.get((anc, model, clean_text(r.get(c_locus, ""))))
        if loc_id is None and c_fuma:
            loc_id = by_fuma.get((anc, model, clean_text(r.get(c_fuma, ""))))
        if loc_id is None and c_chr and c_start and c_end:
            loc_id = by_coord.get((
                anc, model, clean_text(r.get(c_chr, "")),
                clean_text(r.get(c_start, "")), clean_text(r.get(c_end, ""))
            ))
        if loc_id:
            supported.add(loc_id)
    return supported


def pleio_integration(output_dir, supp_tsv, tier_rows, qc):
    if not supp_tsv:
        qc.append({
            "check": "PLEIO_INTEGRATION_INPUTS",
            "ancestry": "",
            "trait": "",
            "status": "SKIPPED",
            "observed": "",
            "expected_or_note": "Set --supp-tsv to build optional PLEIO/LAVA integration tables.",
        })
        return

    locus_path, gene_path, lava_path = integration_inputs(supp_tsv)
    missing = [str(p) for p in (locus_path, gene_path) if not p.exists()]
    if missing:
        qc.append({
            "check": "PLEIO_INTEGRATION_INPUTS",
            "ancestry": "",
            "trait": "",
            "status": "SKIPPED",
            "observed": ";".join(missing),
            "expected_or_note": "Upload/preserve canonical PLEIO TSVs to build overlap tables.",
        })
        return

    loci = normalise_locus_master(read_tsv_rows(locus_path))
    gene_map = normalise_gene_map(read_tsv_rows(gene_path), loci)
    lava_support = lava_supported_loci(read_tsv_rows(lava_path), loci) if lava_path.exists() else set()

    observed = defaultdict(int)
    for l in loci:
        observed[(l["ancestry"], l["pleio_model"])] += 1
    for key, exp in EXPECTED_PLEIO_TOTALS.items():
        got = observed.get(key, 0)
        qc.append({
            "check": "CANONICAL_PLEIO_TOTAL",
            "ancestry": key[0],
            "trait": key[1],
            "status": "PASS" if got == exp else "FAIL",
            "observed": got,
            "expected_or_note": exp,
        })

    # Fast match dictionaries.
    by_ensg = defaultdict(list)
    by_symbol = defaultdict(list)
    for g in gene_map:
        if g["ensg"]:
            by_ensg[(g["ancestry"], g["ensg"])].append(g)
        if g["gene"]:
            by_symbol[(g["ancestry"], norm_gene(g["gene"]))].append(g)

    summary_rows = []
    detail_rows = []
    multilayer_rows = []

    tier_defs = [
        ("BONF_PRIMARY", tier_rows["bonf"]),
        ("STUDYWIDE_BONF", tier_rows["studywide"]),
        ("ULTRA_5E8", tier_rows["ultra"]),
        ("UNIQUE_PROBE_BONF_SENSITIVITY", tier_rows["probe"]),
    ]

    for tier_name, rows in tier_defs:
        unique = unique_gene_rows(rows)
        for anc in ("EAS", "EUR"):
            for model in ("FOURTRAIT", "MDD_3TRAIT", "SCZ_3TRAIT"):
                relevant = MODEL_TRAITS[model]
                eligible = [r for r in unique if r["ancestry"] == anc and r["trait"] in relevant]
                trait_gene_pairs = {(r["trait"], gene_key(r)) for r in eligible}
                unique_genes = {gene_key(r) for r in eligible}
                mapped_genes = set()
                mapped_loci = set()
                prioritised_loci = set()
                both_lava = set()

                for r in eligible:
                    matches = []
                    if r["ensg"]:
                        matches.extend(by_ensg.get((anc, norm_ensg(r["ensg"])), []))
                    if not matches and r["gene"]:
                        matches.extend(by_symbol.get((anc, norm_gene(r["gene"])), []))
                    matches = [m for m in matches if m["pleio_model"] == model]
                    if matches:
                        mapped_genes.add(gene_key(r))
                    for m in matches:
                        mapped_loci.add(m["pleio_locus_id"])
                        if m["pleio_prioritised"]:
                            prioritised_loci.add(m["pleio_locus_id"])
                        if m["pleio_locus_id"] in lava_support:
                            both_lava.add(m["pleio_locus_id"])
                        d = dict(m)
                        d.update({
                            "analysis_tier": tier_name,
                            "SMR_trait": r["trait"],
                            "SMR_gene": r["gene"],
                            "SMR_ensg": r["ensg"],
                            "SMR_probe_id": r["probe_id"],
                            "SMR_qtl_name": r["qtl_name"],
                            "SMR_top_snp": r["top_snp"],
                            "b_SMR": r["b_SMR"],
                            "se_SMR": r["se_SMR"],
                            "p_SMR": r["p_SMR"],
                            "p_SMR_multi": r["p_SMR_multi"],
                            "p_HEIDI": r["p_HEIDI"],
                            "n_HEIDI_SNPs": r["n_HEIDI_SNPs"],
                            "LAVA_overlap_any": bool_text(m["pleio_locus_id"] in lava_support),
                            "interpretation": "DESCRIPTIVE_REGULATORY_CONVERGENCE_NOT_CAUSAL_MEDIATION",
                        })
                        detail_rows.append(d)

                summary_rows.append({
                    "analysis_tier": tier_name,
                    "ancestry": anc,
                    "pleio_model": model,
                    "display_model": DISPLAY_MODEL[(anc, model)],
                    "relevant_SMR_traits": ";".join(sorted(relevant)),
                    "n_unique_SMR_trait_gene_pairs": len(trait_gene_pairs),
                    "n_unique_SMR_genes": len(unique_genes),
                    "n_SMR_genes_mapped_to_PLEIO": len(mapped_genes),
                    "pct_SMR_genes_mapped_to_PLEIO":
                        (100.0 * len(mapped_genes) / len(unique_genes)) if unique_genes else "",
                    "n_PLEIO_loci_with_SMR_gene": len(mapped_loci),
                    "n_PLEIO_prioritised_loci_with_SMR_gene": len(prioritised_loci),
                    "n_loci_with_both_LAVA_and_SMR_support": len(both_lava),
                    "interpretation": "DESCRIPTIVE_REGULATORY_CONVERGENCE_NOT_CAUSAL_MEDIATION",
                })

                for l in [x for x in loci if x["ancestry"] == anc and x["pleio_model"] == model]:
                    has_smr = l["pleio_locus_id"] in mapped_loci
                    has_lava = l["pleio_locus_id"] in lava_support
                    if has_smr or has_lava:
                        z = dict(l)
                        z.update({
                            "analysis_tier": tier_name,
                            "SMR_support_any": bool_text(has_smr),
                            "LAVA_support_any": bool_text(has_lava),
                            "both_SMR_and_LAVA": bool_text(has_smr and has_lava),
                            "interpretation": "DESCRIPTIVE_CONVERGENCE_NOT_INDEPENDENT_REPLICATION",
                        })
                        multilayer_rows.append(z)

    summary_fields = [
        "analysis_tier", "ancestry", "pleio_model", "display_model",
        "relevant_SMR_traits", "n_unique_SMR_trait_gene_pairs", "n_unique_SMR_genes",
        "n_SMR_genes_mapped_to_PLEIO", "pct_SMR_genes_mapped_to_PLEIO",
        "n_PLEIO_loci_with_SMR_gene", "n_PLEIO_prioritised_loci_with_SMR_gene",
        "n_loci_with_both_LAVA_and_SMR_support", "interpretation"
    ]
    detail_fields = [
        "analysis_tier", "ancestry", "pleio_model", "display_model",
        "pleio_locus_id", "pleio_fuma_locus", "chr", "start", "end",
        "pleio_direction_class", "pleio_prioritised", "prioritisation_status",
        "gene", "ensg", "SMR_trait", "SMR_gene", "SMR_ensg", "SMR_probe_id",
        "SMR_qtl_name", "SMR_top_snp", "b_SMR", "se_SMR", "p_SMR",
        "p_SMR_multi", "p_HEIDI", "n_HEIDI_SNPs", "LAVA_overlap_any",
        "interpretation"
    ]
    multilayer_fields = [
        "analysis_tier", "ancestry", "pleio_model", "display_model",
        "pleio_locus_id", "pleio_fuma_locus", "chr", "start", "end",
        "pleio_direction_class", "pleio_prioritised", "prioritisation_status",
        "SMR_support_any", "LAVA_support_any", "both_SMR_and_LAVA",
        "interpretation"
    ]
    write_tsv(Path(output_dir) / "24A_PLEIO_SMR_threshold_summary.tsv", summary_rows, summary_fields)
    write_tsv(Path(output_dir) / "24B_PLEIO_SMR_threshold_detail.tsv", detail_rows, detail_fields)
    write_tsv(Path(output_dir) / "24C_PLEIO_multilayer_thresholds.tsv", multilayer_rows, multilayer_fields)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--smr-root",
        required=True,
        help="Directory containing EAS_SMR_web and EU_webSMR."
    )
    ap.add_argument(
        "--supp-tsv",
        default="",
        help="Optional canonical supplementary TSV directory for PLEIO/LAVA integration."
    )
    ap.add_argument(
        "--output",
        default="SMR_REANALYSIS",
        help="Output directory (default: ./SMR_REANALYSIS)."
    )
    args = ap.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    stats, studywide, qc = first_pass(args.smr_root)
    tier_rows = second_pass(args.smr_root, out, stats, studywide)

    summary = threshold_summary(stats, studywide, tier_rows)
    write_tsv(out / "21_SMR_threshold_comparison.tsv", summary, list(summary[0].keys()))

    bonf_unique = unique_gene_rows(tier_rows["bonf"])
    ultra_unique = unique_gene_rows(tier_rows["ultra"])
    unique_fields = BASE_OUTPUT_COLUMNS + FLAG_COLUMNS + ["unique_gene_rank_within_trait"]
    write_tsv(out / "23A_SMR_bonf_unique_genes.tsv", bonf_unique, unique_fields)
    write_tsv(out / "23B_SMR_ultra_unique_genes.tsv", ultra_unique, unique_fields)

    cross = cross_ancestry_table(bonf_unique)
    cross_fields = [
        "gene", "ensg", "EUR_traits", "EAS_traits",
        "EUR_best_trait", "EUR_best_p_SMR", "EUR_best_p_HEIDI", "EUR_best_b_SMR",
        "EAS_best_trait", "EAS_best_p_SMR", "EAS_best_p_HEIDI", "EAS_best_b_SMR",
        "same_trait_present", "effect_direction_concordant_best"
    ]
    write_tsv(out / "23C_SMR_cross_ancestry_bonf.tsv", cross, cross_fields)

    # Raw-file integrity checks, including the historical EAS MDD/SCZ duplication concern.
    eas_mdd = Path(args.smr_root) / "EAS_SMR_web/EAS_MDD_eSMR.merged.tsv"
    eas_scz = Path(args.smr_root) / "EAS_SMR_web/SCZ_EAS_eSMR.merged.tsv"
    if eas_mdd.exists() and eas_scz.exists():
        same = md5(eas_mdd) == md5(eas_scz)
        qc.append({
            "check": "EAS_MDD_NOT_DUPLICATE_OF_EAS_SCZ",
            "ancestry": "EAS",
            "trait": "MDD",
            "status": "FAIL" if same else "PASS",
            "observed": "same_md5" if same else "different_md5",
            "expected_or_note": "Different hashes required.",
        })

    for s in summary:
        qc.append({
            "check": "BONFERRONI_GAIN_VS_ULTRA",
            "ancestry": s["ancestry"],
            "trait": s["trait"],
            "status": "PASS",
            "observed": "rows:+{0};genes:+{1}".format(
                s["additional_bonf_rows_vs_ultra"],
                s["additional_bonf_unique_genes_vs_ultra"]
            ),
            "expected_or_note": "Descriptive increase after replacing 5e-8 with 0.05/N.",
        })

    pleio_integration(out, args.supp_tsv, tier_rows, qc)

    qc_fields = ["check", "ancestry", "trait", "status", "observed", "expected_or_note"]
    write_tsv(out / "00_SMR_reanalysis_QC.tsv", qc, qc_fields)
    write_tsv(out / "99_SMR_reanalysis_checks.tsv", qc, qc_fields)

    print("\nSMR/HEIDI reanalysis complete.")
    print("Output:", out)
    print("Study-wide threshold: {0:.6g}".format(studywide))
    print("\nPrimary Bonferroni + HEIDI p>0.01:")
    for s in summary:
        print(
            "{0} {1}: threshold={2:.6g}, rows={3}, unique_genes={4}, "
            "additional_vs_5e-8={5} rows/{6} genes".format(
                s["ancestry"], s["trait"], s["bonf_threshold_all_tests"],
                s["bonf_HEIDI001_rows"], s["bonf_HEIDI001_unique_genes"],
                s["additional_bonf_rows_vs_ultra"],
                s["additional_bonf_unique_genes_vs_ultra"],
            )
        )
    print("\nKey file:", out / "21_SMR_threshold_comparison.tsv")


if __name__ == "__main__":
    main()
