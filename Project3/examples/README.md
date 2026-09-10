# Synthetic validation example

The files under `synthetic/` mimic the minimum FUMA and merged-PLEIO inputs required by the directional-classification workflow. All loci, identifiers, and Z scores are artificial.

Run the example from the `Project3` directory:

```bash
bash examples/run_synthetic_example.sh
```

Expected locus labels are:

| Locus | Expected locus class | Reason |
| --- | --- | --- |
| 1 | `Concordant` | EA and CF Z scores are positive and the psychiatric Z scores are negative. |
| 2 | `Mixed` | The psychiatric domain contains both positive and negative qualifying Z scores. |
| 3 | `Discordant` | The cognitive/educational and psychiatric domains are both negative. |
| 4 | `Unassigned` | Neither EA nor CF reaches the default directional-information threshold. |
| 5 | `Dual` | One SNP is Concordant and the other is Discordant, with no Mixed SNP. |

The example also exercises the single-trait overlap audit: loci 1 and 2 overlap
single-trait FUMA loci, while loci 3, 4, and 5 do not.

The `synthetic/lava/` files separately exercise LAVA consolidation. Run:

```bash
python3 tests/test_lava_postprocessing.py
```

That test verifies the exact `0.05 / number of regions` eligibility rule,
retention of a valid zero P value, canonical EA/CF labels when input trait order
is reversed, and BH correction within trait-pair families.
