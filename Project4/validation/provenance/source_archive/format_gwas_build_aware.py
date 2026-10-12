import polars as pl
import argparse
import re
import sys

def find_col(patterns, columns):
    for p in patterns:
        for c in columns:
            if re.search(p, c.lower()):
                return c
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--universe", choices=["singapore", "hm3"], required=True)
    args = parser.parse_args()

    # 1. Load GWAS with broad null handling
    try:
        gwas_lazy = pl.scan_csv(args.gwas, separator="\t", null_values=["NA", "nan", "NaN", "."], ignore_errors=True)
        key_lazy = pl.scan_csv(args.key, separator="\t")
        gwas_cols = gwas_lazy.collect_schema().names()
        key_cols = key_lazy.collect_schema().names()
    except Exception as e:
        print(f"❌ Error reading {args.gwas}: {e}")
        return

    # 2. Broad Header Mapping
    chr_col = find_col([r'^chr$', r'^chromosome$', r'^#chrom$', r'^chrom$'], gwas_cols)
    pos_col = find_col([r'^pos$', r'^bp$', r'^position$', r'^base_pair_location$'], gwas_cols)
    a1_col  = find_col([r'^ref$', r'^a1$', r'^effect_allele$', r'^inc_allele$'], gwas_cols)
    a2_col  = find_col([r'^alt$', r'^a2$', r'^other_allele$', r'^dec_allele$'], gwas_cols)
    beta_col= find_col([r'^beta_meta$', r'^beta$', r'^b$', r'^effect$'], gwas_cols)
    se_col  = find_col([r'^se_meta$', r'^se$', r'^standard_error$'], gwas_cols)
    p_col   = find_col([r'^p_value_meta$', r'^p_value$', r'^p$', r'^pval$', r'^p-value$'], gwas_cols)
    neg_p   = find_col([r'^neglog10_pval_meta$', r'^neglog10_p_meta$', r'^lp$'], gwas_cols)
    n_col   = find_col([r'^n_meta$', r'^n$', r'^sample_size$', r'^total_n$'], gwas_cols)

    # 3. Handle UKB 'variant' column (Splitting into components)
    variant_col = find_col([r'^variant$'], gwas_cols)
    if not chr_col and variant_col:
        gwas_lazy = gwas_lazy.with_columns([
            pl.col(variant_col).str.split(":").list.get(0).alias("CHR_derived"),
            pl.col(variant_col).str.split(":").list.get(1).cast(pl.Int64).alias("POS_derived"),
            pl.col(variant_col).str.split(":").list.get(2).alias("A1_derived"),
            pl.col(variant_col).str.split(":").list.get(3).alias("A2_derived")
        ])
        chr_col, pos_col = "CHR_derived", "POS_derived"
        if not a1_col: a1_col = "A1_derived"
        if not a2_col: a2_col = "A2_derived"

    # 4. Critical Safety Check: Skip if columns are None
    required = {"CHR": chr_col, "POS": pos_col, "A1": a1_col, "A2": a2_col, "BETA": beta_col, "SE": se_col}
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"❌ Skipping {args.gwas}: Missing headers {missing}")
        return

    # 5. P-value Recovery
    if not p_col and neg_p:
        gwas_lazy = gwas_lazy.with_columns(p_derived = 10**(-pl.col(neg_p).cast(pl.Float64)))
        p_col = "p_derived"
    elif not p_col:
        print(f"❌ Skipping {args.gwas}: No P-value column found.")
        return

    # 6. Build Detection & Join
    gwas_lazy = gwas_lazy.with_columns(gwas_key = pl.col(chr_col).cast(pl.Utf8) + ":" + pl.col(pos_col).cast(pl.Utf8))
    best_df, max_snps, det_build = None, 0, "None"

    for build in ["BP_37", "BP_38"]:
        temp_key = key_lazy.with_columns(k_key = pl.col("CHR").cast(pl.Utf8) + ":" + pl.col(build).cast(pl.Utf8))
        joined = gwas_lazy.join(temp_key, left_on="gwas_key", right_on="k_key", how="inner")
        count = joined.select(pl.len()).collect().item()
        if count > max_snps:
            max_snps, best_df, det_build = count, joined, build

    if max_snps < 1000:
        print(f"❌ Skipping {args.gwas}: Low SNP match count ({max_snps})")
        return

    # 7. Universe-Specific ID Selection
    if args.universe == "singapore":
        id_expr = ("chr" + pl.col("CHR").cast(pl.Utf8) + ":" + 
                   pl.col(det_build).cast(pl.Utf8) + ":" + 
                   pl.col("A1") + ":" + pl.col("A2"))
    else:
        rsid_name = "RSID" if "RSID" in key_cols else "SNP"
        id_expr = pl.col(rsid_name)

    n_val = pl.col(n_col) if n_col else pl.lit(400000)
    
    # 8. Final Select and Write
    (best_df
     .with_columns(Z_raw = pl.col(beta_col).cast(pl.Float64) / pl.col(se_col).cast(pl.Float64))
     .with_columns(
         Z = pl.when((pl.col(a1_col).str.to_uppercase() == pl.col("A2")) & 
                     (pl.col(a2_col).str.to_uppercase() == pl.col("A1")))
               .then(pl.col("Z_raw") * -1).otherwise(pl.col("Z_raw")),
         SNP_ID = id_expr
     )
     .select([
         pl.col("SNP_ID").alias("SNP"),
         pl.col("A1"), pl.col("A2"),
         pl.col("Z"), pl.col(p_col).alias("P"), n_val.alias("N")
     ])
     .collect(engine="streaming")
     .write_csv(args.out, separator="\t")
    )
    print(f"✅ {args.universe.upper()} [{det_build}] Matched: {max_snps} SNPs")

if __name__ == "__main__":
    main()
