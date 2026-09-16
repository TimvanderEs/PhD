#!/usr/bin/env python3
import sys
import csv
import gzip
import math

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: python convert_helios_b38_to_rsid37.py <input> <rosetta.tsv.gz> <output.tsv.gz>\n"
    )
    sys.exit(1)

infile = sys.argv[1]
rosetta = sys.argv[2]
outfile = sys.argv[3]

def open_auto(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)

def z_to_p(z):
    return math.erfc(abs(float(z)) / math.sqrt(2.0))

def first_present(fieldnames, candidates):
    for c in candidates:
        if c in fieldnames:
            return c
    return None

# ----------------------------
# Load Rosetta mapping
# ----------------------------
lookup = {}

with open_auto(rosetta, 'rt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        rsid = row['RSID']
        chr_ = row['CHR']
        bp37 = row['BP_37']
        bp38 = row['BP_38']
        a1 = row['A1'].upper()
        a2 = row['A2'].upper()

        key_fwd = "chr{0}:{1}:{2}:{3}".format(chr_, bp38, a1, a2)
        key_rev = "chr{0}:{1}:{2}:{3}".format(chr_, bp38, a2, a1)

        lookup[key_fwd] = (rsid, chr_, bp37)
        lookup[key_rev] = (rsid, chr_, bp37)

# ----------------------------
# Parse input
# ----------------------------
n_in = 0
n_mapped = 0
n_dropped = 0
n_missing_marker = 0
n_missing_stat = 0
n_missing_n = 0
n_bad = 0

seen = set()

with open_auto(infile, 'rt') as fin, open_auto(outfile, 'wt') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    fields = reader.fieldnames

    marker_col = first_present(fields, ['SNP', 'ID', 'MarkerName'])
    a1_col = first_present(fields, ['A1', 'Allele1'])
    a2_col = first_present(fields, ['A2', 'Allele2'])
    n_col = first_present(fields, ['N', 'N_trans'])

    z_col = first_present(fields, ['Z'])
    beta_col = first_present(fields, ['BETA', 'Effect', 'BETA_trans'])
    se_col = first_present(fields, ['SE', 'StdErr', 'SE_trans'])
    p_col = first_present(fields, ['P', 'P.value', 'P-value', 'P_trans'])

    if marker_col is None:
        sys.stderr.write("Missing required marker column. Tried: SNP, ID, MarkerName\n")
        sys.exit(1)

    if a1_col is None or a2_col is None:
        sys.stderr.write("Missing required allele columns. Need A1/A2 or Allele1/Allele2\n")
        sys.exit(1)

    if n_col is None:
        sys.stderr.write("Missing sample size column. Tried: N, N_trans\n")
        sys.exit(1)

    if z_col is None and not (beta_col and se_col):
        sys.stderr.write("Need either Z or a BETA/SE pair.\n")
        sys.exit(1)

    writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
    writer.writerow(['SNP', 'CHR', 'BP', 'A1', 'A2', 'N', 'Z', 'P'])

    for row in reader:
        n_in += 1

        marker = row.get(marker_col, '')
        if marker == '':
            n_missing_marker += 1
            n_dropped += 1
            continue

        a1 = row[a1_col].upper()
        a2 = row[a2_col].upper()

        if marker not in lookup:
            n_dropped += 1
            continue

        rsid, chr_, bp37 = lookup[marker]

        # deduplicate by rsid + alleles
        dedup_key = (rsid, a1, a2)
        if dedup_key in seen:
            continue
        seen.add(dedup_key)

        nval = row.get(n_col, '')
        if nval == '' or nval == 'NA':
            n_missing_n += 1
            n_dropped += 1
            continue

        try:
            n_float = float(nval)
        except Exception:
            n_bad += 1
            n_dropped += 1
            continue

        z = None
        p = None

        # Prefer explicit Z if available
        if z_col and row.get(z_col, '') not in ('', 'NA'):
            try:
                z = float(row[z_col])
            except Exception:
                z = None

        # Otherwise compute Z from beta/se
        if z is None and beta_col and se_col:
            beta_val = row.get(beta_col, '')
            se_val = row.get(se_col, '')
            if beta_val not in ('', 'NA') and se_val not in ('', 'NA'):
                try:
                    beta = float(beta_val)
                    se = float(se_val)
                    if se == 0:
                        n_missing_stat += 1
                        n_dropped += 1
                        continue
                    z = beta / se
                except Exception:
                    z = None

        if z is None:
            n_missing_stat += 1
            n_dropped += 1
            continue

        # Use explicit p if present, else derive from z
        if p_col and row.get(p_col, '') not in ('', 'NA'):
            try:
                p = float(row[p_col])
            except Exception:
                p = None

        if p is None:
            p = z_to_p(z)

        if p == 0.0:
            p = 1e-300

        writer.writerow([
            rsid,
            chr_,
            bp37,
            a1,
            a2,
            int(n_float) if n_float.is_integer() else n_float,
            z,
            p
        ])
        n_mapped += 1

sys.stderr.write("Input rows:          {0}\n".format(n_in))
sys.stderr.write("Mapped rows:         {0}\n".format(n_mapped))
sys.stderr.write("Dropped rows:        {0}\n".format(n_dropped))
sys.stderr.write("Missing marker rows: {0}\n".format(n_missing_marker))
sys.stderr.write("Missing N rows:      {0}\n".format(n_missing_n))
sys.stderr.write("Missing stat rows:   {0}\n".format(n_missing_stat))
sys.stderr.write("Bad parse rows:      {0}\n".format(n_bad))
