# Portable LD-reference workflow · v0.1

Validate harmonised, munged GWAS inputs and unpartitioned LD-score panels, inspect
commands, run h²/rg analyses and export observed-scale h² summaries. The wrapper
uses Python 3.9+ and its standard library. LDSC/cov-LDSC is a separate engine with
its own environment and licence. No controlled data or upstream engine code is
included.

The analysis runner starts from existing LD-score panels. `build_panel.py` builds
new panels from prepared chromosome PLINK files, with optional PC adjustment.
This release is a local prototype, not a claim of a new statistical method.

## Try it without controlled data

Run from this directory, using fresh output paths:

```sh
python3 make_demo.py /tmp/ldref-demo
python3 ldref.py validate /tmp/ldref-demo/config.json
python3 ldref.py plan /tmp/ldref-demo/config.json --output /tmp/ldref-plan
python3 ldref.py compare /tmp/ldref-demo/blocks_a.tsv /tmp/ldref-demo/blocks_b.tsv --h2-a 0.20 --h2-b 0.19
python3 -m unittest -v test_ldref
```

The example fabricates 220 SNPs across 22 chromosomes. It tests input validation,
planning and paired variance arithmetic. It is too small for scientific LDSC
inference and does not run a real engine. Expected comparison: difference 0.01,
paired SE approximately 0.012247. Tests use a clearly labelled fake backend to
check process failures, result parsing and output preservation. A real-engine
public example and a pinned container are still required for publication release.

## Use your inputs

Copy `example.config.json` outside the repository and replace all placeholders.
Paths resolve relative to the configuration file; absolute paths also work.
Panel prefixes contain `{chr}`, for example `/data/panel/{chr}` or
`/data/panel/chr{chr}.reference`. The workflow passes the equivalent `@`
chromosome placeholder to LDSC. All 22 autosomes are required.

```sh
python3 ldref.py validate /data/my-analysis.json
python3 ldref.py plan /data/my-analysis.json --output /data/planned-analysis
python3 ldref.py run /data/my-analysis.json --output /data/new-analysis
python3 ldref.py summarise /data/new-analysis
```

Use separate fresh directories for `plan` and `run`. Neither overwrites earlier
results. Jobs execute sequentially without shell interpolation. Failures stop the
run and are recorded with the process exit code and captured output.

Configuration fields:

| Field | Meaning |
|---|---|
| `engine.python`, `engine.script` | Interpreter and upstream LDSC/cov-LDSC entry point |
| `engine.version` | Verified upstream release or commit; script hash is captured too |
| `build` | Declared build shared by GWAS, panels and weights |
| `panels` | Named reference prefixes, ancestry and source release |
| `weights` | Fixed regression-weight LD-score prefix used for all comparisons |
| `traits` | Munged files, ancestry, source accession/phenotype, observed scale and maximum plausible N |
| `analyses` | `h2`, `rg` or both; rg generates all trait pairs within each reference |
| `n_blocks` | Requested number of LDSC jackknife blocks |
| `m_mode` | `common` uses `.l2.M_5_50`; `all` uses `.l2.M` and `--not-M-5-50` |

Do not substitute a nominal cohort size for variant-level N. Binary traits are
reported on the observed scale; liability conversion is outside v0.1.

## Build a panel

Copy `build.config.json` and supply genotype paths, source, build and engine.
Complete sample selection, relatedness filtering, variant QC and interpolation
of a build-matched genetic map first. BIM column 3 must contain genetic positions
in centimorgans. The [historical construction workflow](../WORKFLOW.md) provides
the Singapore-specific preprocessing record. The portable builder consumes its
prepared inputs; it does not infer QC choices or regenerate PCs automatically.

```sh
python3 build_panel.py /data/build.json --output /data/build-plan
python3 build_panel.py /data/build.json --output /data/new-panel --execute
```

Set `covariates` to null for standard LD scores or to a headerless `FID IID PC1
PC2 ...` file for cov-LDSC. A `{chr}` placeholder is optional in the covariate
path. The builder checks BED dimensions, FAM IDs, map order, SNP alleles and
sample consistency across chromosomes. Covariates must match FAM order exactly
and be centred and orthonormal. Unit length alone is insufficient if PC columns
are correlated after sample restriction. Recompute PCs on the chosen reference
sample when needed. A scientific rationale for window size and PC count remains
the researcher's responsibility.

Build outputs include the native per-chromosome panel files and `build.json`
with commands, hashes and completion status. Run `ldref.py validate` with an
analysis configuration pointing at the new panel before association-summary
analyses. Construction completion alone does not establish calibration.

## Input checks and interpretation

The validator reads every row, checks finite LD scores/Z/N, valid allele letters,
duplicates, ordered coordinates, complete chromosome bundles, positive SNP-count
files, reference-to-weight overlap and N against the declared maximum. Compared
reference panels must have identical SNP IDs, coordinates and order. Low overlap,
weak mean chi-square and declared ancestry differences generate warnings.

The validator cannot establish phenotype identity or genome build from rsIDs,
verify allele alignment against genotypes, or determine whether a reference's
covariates match the GWAS. Source declarations require independent verification.
Negative LD scores are allowed because finite-sample corrected estimates may be
negative. SNP-count files describe the variants contributing to the heritability
model and must not be regenerated from the number of regression rows. Validation
checks their format, not the scientific correctness of those counts.

v0.1 supports a single `L2` column only. It does not support partitioned LDSC,
multiple weight panels, custom intercepts or automatic build conversion. rg
estimates appear in `rg_summary.tsv`; non-estimable results retain an explicit
status and link by job name to the native logs. P values in that summary are
normal approximations recomputed from the printed estimate and SE, and may
differ slightly from upstream P values due to rounding.
An ancestry mismatch warning may be expected in a benchmark, but is not a
justification for using a pooled panel for an ancestry-specific GWAS.

## Outputs

Each run stores `config.resolved.json`, `validation.json` (input SHA-256 hashes
and warnings), `commands.json`, `run.json`, per-job stdout and native engine logs.
Successful h² jobs produce `h2_summary.tsv` with h², SE, 95% normal-approximation
intervals and h²/SE. Negative estimates and interval limits are retained.
Resolved paths may identify controlled infrastructure: review before publishing.

## Paired comparisons

`compare` accepts two whitespace-delimited tables with columns:

```
block omitted_snps_sha256 delete_h2
```

The block hash must represent the **actual omitted regression SNPs** from the
engine, with an identical serialization convention in both runs. Do not invent
block hashes from ordinal block numbers. The command checks equal block IDs,
order and SNP-membership hashes. Raw `.delete` files do not contain this evidence,
so v0.1 intentionally does not automatically convert them into paired tables.
The caller must also establish matching model scale, SNP universe and settings.

The full-sample contrast is h²(A) − h²(B). Paired SE uses full-precision delete
differences and the `(B−1)/B` jackknife factor. It never bias-corrects a rounded log
estimate by multiplying it by B. The P value is an unadjusted normal approximation;
prespecify and correct for multiple contrasts/traits before making claims.
Zero SE returns no P value. Non-significance is not an equivalence test.

## Next release gates

- Pin and exercise an upstream engine environment/container on public data.
- Exercise the builder on real public genotypes and add numerical tests for the
  local vectorisation patch. Format/order checks are already included.
- Add native block-membership capture before automatic paired inference.
- Add explicit multiple-testing families for inference across traits/panels.
- Validate approved reference bundles and their separate distribution terms.

Credit the [original LDSC](https://github.com/bulik/ldsc) and
[cov-LDSC](https://github.com/immunogenomics/cov-ldsc) implementations. The wrapper
uses the repository's MIT licence; upstream engines retain their own licences.
