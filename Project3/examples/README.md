# Synthetic validation example

The files under `synthetic/` mimic the minimum FUMA and merged-PLEIO inputs required by the directional-classification workflow. All loci, identifiers, and Z scores are artificial.

Run the example from the `Project3` directory:

```bash
bash examples/run_synthetic_example.sh
```

Expected locus labels are:

| Locus | Expected corrected label | Reason |
| --- | --- | --- |
| 1 | `concordant` | Cognitive Z scores are positive and psychiatric Z scores are negative. |
| 2 | `mixed` | The psychiatric domain contains both positive and negative qualifying Z scores. |
| 3 | `discordant` | Cognitive and psychiatric domains are both negative. |
| 4 | `trait_dominant_or_unclear` | No cognitive Z score reaches the default threshold. |

The example also exercises the single-trait overlap audit: loci 1 and 2 overlap single-trait FUMA loci, while loci 3 and 4 do not.
