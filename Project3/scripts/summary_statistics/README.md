# Summary-statistic preparation

- `convert_regenie_to_ldsc.sh` converts one REGENIE association file to the
  column format used for LDSC preparation. It is a path-independent version of
  the recovered EC2 script.
- `convert_helios_b38_to_rsid37.py` maps HELIOS build 38 variant identifiers to
  rsIDs and build 37 positions using a project-specific Rosetta table. The copy
  here is byte-identical to the executed local and EC2 copies.

The final build-37 freeze and downstream exports are fingerprinted in
[`../../results/summary_statistics/freeze_manifest.tsv`](../../results/summary_statistics/freeze_manifest.tsv).
The Rosetta table and full summary statistics are not distributed in this code
repository.

```bash
./convert_regenie_to_ldsc.sh input.regenie output.sumstats.tsv.gz
python convert_helios_b38_to_rsid37.py \
  input.tsv.gz rosetta.tsv.gz output.tsv.gz
```
