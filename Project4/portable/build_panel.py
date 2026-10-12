#!/usr/bin/env python3
"""Build standard or PC-adjusted LD scores from QC'd chromosome PLINK files."""
import argparse
import json
import subprocess
from pathlib import Path

from ldref import number, require, sha, write_json


def prepare(config_path, output, execute=False):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    require(config.get("schema_version") == 1, "Expected schema_version 1")
    require(config.get("build") in ("GRCh37", "GRCh38"), "Declare build")
    require(config.get("ancestry") and config.get("source"), "Declare ancestry and source")
    require(config.get("qc_complete") is True, "Inputs must have completed genotype/sample QC")
    require(config.get("engine", {}).get("version"), "Declare engine commit/version")
    window = number(config["window_cm"])
    require(window > 0, "Positive window_cm required")
    prefix_template = str((config_path.parent / config["bfile"]).resolve())
    require(prefix_template.count("{chr}") == 1, "bfile requires one {chr} placeholder")
    engine = config["engine"]
    interpreter = engine["python"]
    if "/" in interpreter or "\\" in interpreter:
        interpreter = str((config_path.parent / interpreter).resolve())
    script = str((config_path.parent / engine["script"]).resolve())
    cov_template = config.get("covariates")
    if cov_template:
        cov_template = str((config_path.parent / cov_template).resolve())
    manifest, commands = {}, []
    previous_ids = None
    output = Path(output).resolve()
    require(not output.exists(), "Use a fresh output directory")
    for chrom in range(1, 23):
        prefix = prefix_template.replace("{chr}", str(chrom))
        ids = []
        for line in Path(prefix+".fam").read_text().splitlines():
            fields = line.split()
            require(len(fields) == 6, "Malformed FAM")
            ids.append(tuple(fields[:2]))
        require(len(ids) >= 3 and len(set(ids)) == len(ids), "Empty/small sample or duplicate FAM IDs")
        require(previous_ids is None or ids == previous_ids, "Sample order differs across chromosomes")
        previous_ids = ids
        last_bp, last_cm, markers, cm_values = 0, -1, set(), set()
        for line in Path(prefix+".bim").read_text().splitlines():
            fields = line.split()
            require(len(fields) == 6 and int(fields[0]) == chrom, "Wrong chromosome/malformed BIM")
            cm, bp = number(fields[2]), int(fields[3])
            require(bp > 0 and bp >= last_bp and cm >= 0 and cm >= last_cm, "Invalid or unordered map")
            require(fields[1] not in markers, "Duplicate BIM SNP")
            require(fields[4] in ("A","C","G","T") and fields[5] in ("A","C","G","T") and fields[4] != fields[5], "Non-SNP BIM alleles")
            markers.add(fields[1]); cm_values.add(cm)
            last_bp, last_cm = bp, cm
        require(len(cm_values) > 1, "Genetic map is constant/missing; interpolate a build-matched map first")
        with Path(prefix+".bed").open("rb") as handle:
            require(handle.read(3) == b"\x6c\x1b\x01", "Expected SNP-major PLINK BED")
        require(Path(prefix+".bed").stat().st_size == 3 + ((len(ids)+3)//4)*len(markers), "BED size incompatible with BIM/FAM")
        for suffix in (".bed", ".bim", ".fam"):
            manifest[prefix+suffix] = sha(prefix+suffix)
        command = [interpreter, script, "--bfile", prefix, "--l2", "--ld-wind-cm", str(window)]
        if cov_template:
            covfile = cov_template.replace("{chr}", str(chrom))
            covrows = [line.split() for line in Path(covfile).read_text().splitlines()]
            require([tuple(row[:2]) for row in covrows] == ids, "Covariate IDs/order do not match FAM; headerless FID IID PCs required")
            k = len(covrows[0])-2
            require(0 < k < len(ids)-1 and all(len(row) == k+2 for row in covrows), "Invalid PC dimensions")
            columns = [[number(row[j+2]) for row in covrows] for j in range(k)]
            for j, col in enumerate(columns):
                require(abs(sum(col)) < 1e-5, "PC columns must be centred")
                for i in range(j+1):
                    dot = sum(a*b for a,b in zip(col, columns[i]))
                    require(abs(dot - (1 if i == j else 0)) < 1e-5, "PC columns must be orthonormal; recompute PCs after sample selection")
            manifest[covfile] = sha(covfile)
            command += ["--cov", covfile]
        command += ["--out", str(output / str(chrom))]
        commands.append(command)
    output.mkdir(parents=True)
    state = dict(status="planned", config=config, files=manifest, commands=commands, samples=len(previous_ids),
                 builder_sha256=sha(__file__))
    write_json(output / "build.json", state)
    if execute:
        state["status"] = "running"
        try:
            state["engine_sha256"] = sha(script)
            for chrom, command in enumerate(commands, 1):
                with (output / f"{chrom}.stdout.log").open("w") as log:
                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
                require(result.returncode == 0, f"Chromosome {chrom} failed; inspect stdout log")
                for suffix in (".l2.ldscore.gz", ".l2.M", ".l2.M_5_50"):
                    require((output / f"{chrom}{suffix}").is_file(), f"Missing chromosome {chrom} output {suffix}")
                state["completed_chromosomes"] = chrom
                write_json(output / "build.json", state)
            state["status"] = "complete"
        except BaseException:
            state["status"] = "failed"
            raise
        finally:
            write_json(output / "build.json", state)
    return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", help="Run engine after checks; default only writes a plan")
    args = parser.parse_args()
    try:
        state = prepare(args.config, args.output, args.execute)
        print(state["status"])
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, f"ERROR: {error}\n")
