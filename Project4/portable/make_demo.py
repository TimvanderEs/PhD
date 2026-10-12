#!/usr/bin/env python3
"""Generate synthetic input-format fixtures. These are NOT scientific results."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def create(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    for panel in ("panel_a", "panel_b", "weights"):
        directory = root / panel
        directory.mkdir()
        for chrom in range(1, 23):
            with gzip.open(directory / f"{chrom}.l2.ldscore.gz", "wt") as handle:
                handle.write("CHR SNP BP L2\n")
                for i in range(1, 11):
                    handle.write(f"{chrom} rs{chrom}_{i} {i*1000} {1+i/10+(0.1 if panel == 'panel_b' else 0)}\n")
            for suffix in ("M", "M_5_50"):
                (directory / f"{chrom}.l2.{suffix}").write_text("10\n")
    with gzip.open(root / "trait.sumstats.gz", "wt") as handle:
        handle.write("SNP A1 A2 N Z\n")
        for chrom in range(1, 23):
            for i in range(1, 11):
                handle.write(f"rs{chrom}_{i} A C 1000 {i/10}\n")
    config = dict(schema_version=1, build="GRCh38", m_mode="all", n_blocks=4,
                  analyses=["h2"], engine=dict(python="python2", script="external/ldsc.py", version="demo-no-engine"),
                  weights=dict(prefix="weights/{chr}", build="GRCh38"),
                  panels={p: dict(prefix=p+"/{chr}", build="GRCh38", ancestry="synthetic", source="fabricated format fixture") for p in ("panel_a", "panel_b")},
                  traits=dict(toy=dict(sumstats="trait.sumstats.gz", build="GRCh38", ancestry="synthetic", source="fabricated format fixture", scale="observed", max_n=1000)))
    (root / "config.json").write_text(json.dumps(config, indent=2)+"\n")
    for name, values in (("a", [0.19, 0.20, 0.21, 0.20]), ("b", [0.18, 0.19, 0.19, 0.20])):
        with (root / f"blocks_{name}.tsv").open("w") as handle:
            handle.write("block omitted_snps_sha256 delete_h2\n")
            for i, value in enumerate(values):
                digest = hashlib.sha256(f"synthetic-block-{i}".encode()).hexdigest()
                handle.write(f"{i} {digest} {value}\n")
    return root / "config.json"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(create(parser.parse_args().output))
