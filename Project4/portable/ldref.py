#!/usr/bin/env python3
"""Portable LDSC runner. Python 3 standard library only; engines are external."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

VERSION = "0.1.0"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = next(handle).split()
        require(len(set(header)) == len(header), f"Duplicate headers: {path}")
        for line_no, line in enumerate(handle, 2):
            fields = line.split()
            require(len(fields) == len(header), f"Wrong field count: {path}:{line_no}")
            yield dict(zip(header, fields))


def number(value):
    result = float(value)
    require(math.isfinite(result), f"Non-finite numeric value: {value}")
    return result


def load(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    require(config.get("schema_version") == 1, "Expected schema_version 1")
    require(config.get("build") in ("GRCh37", "GRCh38"), "Declare GRCh37 or GRCh38")
    require(config.get("m_mode") in ("common", "all"), "m_mode must be common or all")
    require(type(config.get("n_blocks")) is int and config["n_blocks"] >= 2, "n_blocks must be an integer >= 2")
    require(config.get("panels") and config.get("traits"), "At least one panel and trait required")
    require(config.get("analyses") and set(config["analyses"]) <= {"h2", "rg"}, "analyses must contain h2 and/or rg")
    require(len(config["analyses"]) == len(set(config["analyses"])), "Duplicate analyses")
    require("weights" not in config["panels"], "Panel name weights is reserved")
    require("rg" not in config["analyses"] or len(config["traits"]) >= 2, "rg requires at least two traits")
    for group in ("panels", "traits"):
        for name, entry in config[group].items():
            require(re.fullmatch(r"[A-Za-z0-9_-]+", name), f"Unsafe identifier: {name}")
            require(entry.get("build") == config["build"], f"Build mismatch: {name}")
            require(entry.get("ancestry") and entry.get("source"), f"Declare ancestry and source: {name}")
    for entry in [*config["panels"].values(), config["weights"]]:
        require(entry.get("build") == config["build"], "Weight/reference build mismatch")
        require("{chr}" in entry["prefix"], "Panel prefixes require {chr}")
        require(entry["prefix"].count("{chr}") == 1, "Exactly one chromosome placeholder required")
        entry["prefix"] = str((path.parent / entry["prefix"]).resolve())
    for entry in config["traits"].values():
        entry["sumstats"] = str((path.parent / entry["sumstats"]).resolve())
        require(entry.get("scale") == "observed", "v0.1 reports observed-scale h2 only")
    config["engine"]["script"] = str((path.parent / config["engine"]["script"]).resolve())
    engine_python = config["engine"]["python"]
    if "/" in engine_python or "\\" in engine_python:
        config["engine"]["python"] = str((path.parent / engine_python).resolve())
    require(config["engine"].get("version"), "Declare engine version/commit")
    return config


def validate(config):
    manifest, warnings = {}, []
    reference_ids = None
    weight_ids = set()
    for name, panel in [*config["panels"].items(), ("weights", config["weights"])]:
        ids, seen = [], set()
        for chrom in range(1, 23):
            prefix = panel["prefix"].replace("{chr}", str(chrom))
            file = prefix + ".l2.ldscore.gz"
            count = 0
            last_bp = 0
            for row in rows(file):
                require(set(row) == {"CHR", "SNP", "BP", "L2"}, f"Expected unpartitioned CHR SNP BP L2: {file}")
                require(int(row["CHR"]) == chrom, f"Wrong chromosome: {file}")
                bp = int(row["BP"])
                require(bp >= last_bp and bp > 0, f"Coordinates not ordered: {file}")
                last_bp = bp
                number(row["L2"])
                require(row["SNP"] not in seen, f"Duplicate panel SNP: {row['SNP']}")
                seen.add(row["SNP"])
                ids.append((chrom, bp, row["SNP"]))
                count += 1
            require(count, f"Empty chromosome file: {file}")
            manifest[file] = sha(file)
            if name != "weights":
                suffix = ".l2.M_5_50" if config["m_mode"] == "common" else ".l2.M"
                mfile = prefix + suffix
                values = Path(mfile).read_text().split()
                require(len(values) == 1 and number(values[0]) > 0, f"Invalid SNP-count file: {mfile}")
                manifest[mfile] = sha(mfile)
        if name == "weights":
            weight_ids = seen
        elif reference_ids is None:
            reference_ids = ids
        else:
            require(ids == reference_ids, f"Panel {name} has a different SNP universe/order/coordinates; harmonise explicitly first")
    common = {item[2] for item in reference_ids} & weight_ids
    require(common, "No reference/weight overlap")
    summaries = {}
    for name, trait in config["traits"].items():
        seen, overlap, n_min, n_max, chi = set(), 0, math.inf, 0, 0
        for row in rows(trait["sumstats"]):
            require({"SNP", "A1", "A2", "N", "Z"} <= set(row), "Munged sumstats require SNP A1 A2 N Z")
            require(row["SNP"] not in seen, f"Duplicate sumstats SNP: {name}/{row['SNP']}")
            seen.add(row["SNP"])
            require(row["A1"] in "ACGT" and len(row["A1"]) == 1 and row["A2"] in "ACGT" and len(row["A2"]) == 1 and row["A1"] != row["A2"], f"Invalid alleles: {name}/{row['SNP']}")
            n, z = number(row["N"]), number(row["Z"])
            require(n > 0, f"Nonpositive N: {name}")
            n_min, n_max = min(n_min, n), max(n_max, n)
            if row["SNP"] in common:
                overlap += 1
                chi += z*z
        require(overlap >= config["n_blocks"], f"Insufficient overlapping SNPs: {name}")
        require(n_max <= number(trait["max_n"]), f"N exceeds declared study maximum: {name}")
        mean_chi = chi / overlap
        summaries[name] = dict(rows=len(seen), overlap=overlap, n_min=n_min, n_max=n_max, mean_chi2_before_engine_filter=mean_chi)
        if overlap < 200000:
            warnings.append(f"{name}: fewer than 200,000 overlapping SNPs")
        if mean_chi < 1.02:
            warnings.append(f"{name}: weak mean chi-square; h2 and rg may be imprecise")
        for panel_name, panel in config["panels"].items():
            if trait["ancestry"] != panel["ancestry"]:
                warnings.append(f"{name}/{panel_name}: declared ancestries differ; justify this comparison")
        manifest[trait["sumstats"]] = sha(trait["sumstats"])
    return dict(status="passed", files=manifest, traits=summaries, warnings=warnings,
                caveat="Build, ancestry and phenotype identity are declarations; LD-score files cannot verify allele alignment or source identity.")


def jobs(config, output):
    engine = [config["engine"]["python"], config["engine"]["script"]]
    for panel, entry in config["panels"].items():
        base = ["--ref-ld-chr", entry["prefix"].replace("{chr}", "@"),
                "--w-ld-chr", config["weights"]["prefix"].replace("{chr}", "@"),
                "--n-blocks", str(config["n_blocks"]), "--print-delete-vals"]
        if config["m_mode"] == "all":
            base += ["--not-M-5-50"]
        for analysis in config["analyses"]:
            combinations = [(t,) for t in config["traits"]] if analysis == "h2" else itertools.combinations(config["traits"], 2)
            for traits in combinations:
                prefix = output / analysis / panel / "__".join(traits)
                command = engine + ["--" + analysis, ",".join(config["traits"][t]["sumstats"] for t in traits)] + base + ["--out", str(prefix)]
                yield dict(analysis=analysis, panel=panel, traits=list(traits), prefix=str(prefix), argv=command)


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def write_table(path, records):
    require(records, "No records to write")
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def summarise(output):
    state = json.loads((output / "run.json").read_text())
    require(state["status"] == "complete", "Run did not complete; inspect run.json and logs")
    records, correlations = [], []
    for job in state["jobs"]:
        log = Path(job["prefix"] + ".log").read_text()
        if job["analysis"] != "h2":
            match = re.search(r"Genetic Correlation:\s+([-+\d.eE]+)\s+\(([-+\d.eE]+)\)", log)
            if match:
                rg, se = map(number, match.groups())
                require(se > 0, "Nonpositive rg SE")
                correlations.append(dict(panel=job["panel"], trait_a=job["traits"][0], trait_b=job["traits"][1],
                                         rg=rg, se=se, status="estimated", p_normal=math.erfc(abs(rg/se)/math.sqrt(2))))
            else:
                require("Summary of Genetic Correlation Results" in log, f"Incomplete rg output: {job['prefix']}")
                correlations.append(dict(panel=job["panel"], trait_a=job["traits"][0], trait_b=job["traits"][1],
                                         rg="NA", se="NA", status="not_estimable_see_log", p_normal="NA"))
            continue
        match = re.search(r"Total Observed scale h2:\s+([-+\d.eE]+)\s+\(([-+\d.eE]+)\)", log)
        require(match, f"No observed-scale h2 result: {job['prefix']}")
        h2, se = map(number, match.groups())
        require(se > 0, "Nonpositive h2 SE")
        records.append(dict(panel=job["panel"], trait=job["traits"][0], h2=h2, se=se,
                            ci_lower=h2-1.96*se, ci_upper=h2+1.96*se, h2_over_se=h2/se))
    if records:
        write_table(output / "h2_summary.tsv", records)
    if correlations:
        write_table(output / "rg_summary.tsv", correlations)


def paired(path_a, path_b, full_a, full_b):
    a, b = list(rows(path_a)), list(rows(path_b))
    require(len(a) == len(b) and len(a) >= 2, "Unequal or insufficient blocks")
    keys = [r["block"] for r in a]
    require(len(set(keys)) == len(keys), "Duplicate block identifiers")
    require(keys == [r["block"] for r in b], "Block order differs")
    require(all(r["omitted_snps_sha256"] == s["omitted_snps_sha256"] and re.fullmatch(r"[0-9a-f]{64}", r["omitted_snps_sha256"]) for r,s in zip(a,b)), "Block SNP identities differ or are invalid")
    diffs = [number(r["delete_h2"]) - number(s["delete_h2"]) for r,s in zip(a,b)]
    mean = sum(diffs) / len(diffs)
    se = math.sqrt((len(diffs)-1)/len(diffs)*sum((d-mean)**2 for d in diffs))
    delta = number(full_a) - number(full_b)
    return dict(difference=delta, paired_se=se, ci_lower=delta-1.96*se, ci_upper=delta+1.96*se,
                p=math.erfc(abs(delta/se)/math.sqrt(2)) if se else None, blocks=len(diffs),
                caution="Block membership must come from the engine; matching block counts alone is insufficient. No multiplicity correction applied.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("validate", "plan", "run"):
        p = sub.add_parser(action)
        p.add_argument("config", type=Path)
        p.add_argument("--output", type=Path, required=action != "validate", default=Path("ldref-output"))
    p = sub.add_parser("summarise")
    p.add_argument("output", type=Path)
    p = sub.add_parser("compare")
    p.add_argument("a", type=Path)
    p.add_argument("b", type=Path)
    p.add_argument("--h2-a", type=float, required=True)
    p.add_argument("--h2-b", type=float, required=True)
    args = parser.parse_args()
    if args.action == "compare":
        print(json.dumps(paired(args.a, args.b, args.h2_a, args.h2_b), indent=2))
        return
    if args.action == "summarise":
        summarise(args.output.resolve())
        return
    config = load(args.config)
    report = validate(config)
    if args.action == "validate":
        print(json.dumps(report, indent=2))
        return
    output = args.output.resolve()
    require(not output.exists(), "Output exists; use a new directory to preserve prior runs")
    output.mkdir(parents=True)
    plan = list(jobs(config, output))
    require(plan, "No jobs requested")
    write_json(output / "config.resolved.json", config)
    write_json(output / "validation.json", report)
    write_json(output / "commands.json", plan)
    if args.action == "plan":
        print(f"Validated plan: {output / 'commands.json'}")
        return
    require(shutil.which(config["engine"]["python"]), "Engine Python not found")
    engine_path = config["engine"]["script"]
    require(Path(engine_path).is_file(), "Engine script not found")
    state = dict(status="running", workflow_version=VERSION, workflow_sha256=sha(__file__),
                 engine_sha256=sha(engine_path), jobs=[])
    write_json(output / "run.json", state)
    try:
        for job in plan:
            prefix = Path(job["prefix"])
            prefix.parent.mkdir(parents=True, exist_ok=True)
            with Path(str(prefix) + ".stdout.log").open("w") as log:
                result = subprocess.run(job["argv"], stdout=log, stderr=subprocess.STDOUT, check=False)
            job["returncode"] = result.returncode
            state["jobs"].append(job)
            require(result.returncode == 0, f"Engine failed: {prefix}.stdout.log")
            require(Path(str(prefix)+".log").is_file(), f"Engine did not create result log: {prefix}")
            write_json(output / "run.json", state)
        state["status"] = "complete"
    except BaseException:
        state["status"] = "failed"
        raise
    finally:
        write_json(output / "run.json", state)
    try:
        summarise(output)
    except BaseException:
        state["status"] = "failed_result_validation"
        write_json(output / "run.json", state)
        raise
    print(f"Completed: {output}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, StopIteration, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
