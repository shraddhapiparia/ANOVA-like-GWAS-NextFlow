"""
Generates the tiny synthetic PLINK + phenotype dataset used for smoke-testing.
Outputs: example.bed, example.bim, example.fam, example_pheno.csv
Run from this directory: python3 make_example_data.py

Design constraints:
  - 25 samples (> 15 model params needed by snp_anova.R's ANOVA model)
  - 10 SNPs (fits in one batch with start=1, end=10)
  - 3 Cluster levels (~8-9 samples each) for the factor ANOVA
  - Phenotype rows are in the same order as the FAM file (no merge in snp_anova.R)
  - Seed=42 for reproducibility
"""

import random, struct, csv, os

random.seed(42)
OUT = os.path.dirname(os.path.abspath(__file__))

N_SAMPLES = 25
N_SNPS    = 10

# .fam  (FID IID PaternalID MaternalID Sex Phenotype)
fam_ids = [f"S{i}" for i in range(1, N_SAMPLES+1)]
with open(f"{OUT}/example.fam", "w") as f:
    for sid in fam_ids:
        sex = random.choice([1, 2])
        f.write(f"1 {sid} 0 0 {sex} -9\n")

# .bim  (chr snp_id cM position A1 A2)
snp_ids = [f"rs{i}" for i in range(1, N_SNPS+1)]
alleles = [("A","G"),("A","C"),("T","C"),("G","T"),("A","T"),
           ("C","G"),("A","G"),("C","T"),("A","C"),("G","T")]
with open(f"{OUT}/example.bim", "w") as f:
    for idx, sid in enumerate(snp_ids):
        pos = (idx+1) * 100000
        a1, a2 = alleles[idx]
        f.write(f"1\t{sid}\t0\t{pos}\t{a1}\t{a2}\n")

# .bed  (PLINK SNP-major binary format)
# Magic bytes: 0x6c 0x1b 0x01 (SNP-major mode)
# Genotype bit pairs per sample, LSB-first in each byte:
#   00 = hom A1/A1  01 = missing  10 = het  11 = hom A2/A2
BYTES_PER_SNP = (N_SAMPLES + 3) // 4
geno_map = {0: 0b00, 1: 0b10, 2: 0b11}

bed_bytes = bytearray([0x6c, 0x1b, 0x01])
for snp_idx in range(N_SNPS):
    snp_genos = []
    for s in range(N_SAMPLES):
        r = random.random()
        snp_genos.append(0 if r < 0.50 else (1 if r < 0.85 else 2))
    for b in range(BYTES_PER_SNP):
        byte_val = 0
        for bit in range(4):
            s = b*4 + bit
            g = snp_genos[s] if s < N_SAMPLES else 0
            byte_val |= geno_map[g] << (bit * 2)
        bed_bytes.append(byte_val)

with open(f"{OUT}/example.bed", "wb") as f:
    f.write(bed_bytes)

# phenotype CSV  (row order matches .fam)
clusters = [1]*8 + [2]*9 + [3]*8
random.shuffle(clusters)

cols = ["IID","Cluster","Age","Sex"] + [f"PC{i}" for i in range(1,11)] + ["pc1_clinical"]
with open(f"{OUT}/example_pheno.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=cols)
    writer.writeheader()
    for i, sid in enumerate(fam_ids):
        row = {"IID": sid, "Cluster": clusters[i],
               "Age": random.randint(30, 70), "Sex": random.choice([0, 1])}
        for pc in range(1, 11):
            row[f"PC{pc}"] = round(random.gauss(0, 0.05), 6)
        row["pc1_clinical"] = round(clusters[i] * 0.3 + random.gauss(0, 0.5), 6)
        writer.writerow(row)

print(f"Generated {N_SAMPLES} samples, {N_SNPS} SNPs -> {os.listdir(OUT)}")
