"""
Custom code 3: known-gene, protein-coding subset of the SQANTI3-filtered isoforms
=================================================================================
Usage : python 03_coding_subset.py <filter_result_RulesFilter_result_classification.txt>
Input : classification table written by the SQANTI3 rules filter
Output (written next to the input):
  1. known_gene_associated_isoforms.tsv  - FSM/ISM/NIC/NNC/Fusion isoforms of annotated (known) genes
  2. known_gene_coding_isoforms.tsv      - of those, the isoforms with coding == "coding"
"""

import pandas as pd
import os
import sys

# ──────────────────────────────────────────────
# 0. Paths
# ──────────────────────────────────────────────
INPUT_FILE = sys.argv[1] if len(sys.argv) > 1 else (
    "03_sqanti3/filter_out/"
    "filter_result_RulesFilter_result_classification.txt"
)
OUTPUT_DIR = os.path.dirname(INPUT_FILE)
OUT_ALL    = os.path.join(OUTPUT_DIR, "known_gene_associated_isoforms.tsv")
OUT_CODING = os.path.join(OUTPUT_DIR, "known_gene_coding_isoforms.tsv")

# ──────────────────────────────────────────────
# 1. Read the table
# ──────────────────────────────────────────────
print(">> Reading classification file ...")
df = pd.read_csv(INPUT_FILE, sep="\t", low_memory=False)
print(f"   Total rows       : {len(df):,}")
print(f"   Total columns    : {len(df.columns)}")

# strip spaces from column names
df.columns = df.columns.str.strip()

# the "coding" column must be column 30 (index 29)
print(f"\n>> Column [29] = '{df.columns[29]}'")
if df.columns[29] != "coding":
    print(f"[WARNING] Expected 'coding' at col[29], got '{df.columns[29]}'")
    print("          check the column order.")
    sys.exit(1)

# ──────────────────────────────────────────────
# 2. Keep filter_result == "Isoform"
# ──────────────────────────────────────────────
df_isoform = df[df["filter_result"] == "Isoform"].copy()
print(f"\n>> filter_result == 'Isoform' : {len(df_isoform):,} rows")
print(f"   (removed, e.g. Artifact: {len(df) - len(df_isoform):,} rows)")

# ──────────────────────────────────────────────
# 3. Structural categories kept (FSM/ISM/NIC/NNC/Fusion)
#    compared in lower case with the SQANTI3 category names
# ──────────────────────────────────────────────
TARGET_CATS = {
    "full-splice_match",
    "incomplete-splice_match",
    "novel_in_catalog",
    "novel_not_in_catalog",
    "fusion",
}

df_isoform["_cat_lower"] = df_isoform["structural_category"].str.strip().str.lower()
df_cat = df_isoform[df_isoform["_cat_lower"].isin(TARGET_CATS)].copy()

# category names -> short labels (for plots and output)
cat_map = {
    "full-splice_match"      : "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog"       : "NIC",
    "novel_not_in_catalog"   : "NNC",
    "fusion"                 : "Fusion",
}
df_cat["structural_category"] = df_cat["_cat_lower"].map(cat_map)

print(f"\n>> After category filter (FSM/ISM/NIC/NNC/Fusion): {len(df_cat):,} rows")
print(df_cat["structural_category"].value_counts().to_string())

# ──────────────────────────────────────────────
# 4. Known / novel gene by associated_gene
#    novel genes carry the "novelGene_" prefix
#    e.g. novelGene_FAM43B_AS, novelGene_51394
# ──────────────────────────────────────────────
def is_novel_gene(gene_name):
    """novel gene if the name starts with novelGene_ or is missing"""
    if pd.isna(gene_name):
        return True
    return str(gene_name).strip().lower().startswith("novelgene")

df_cat["_is_novel_gene"] = df_cat["associated_gene"].apply(is_novel_gene)

n_known = (~df_cat["_is_novel_gene"]).sum()
n_novel = df_cat["_is_novel_gene"].sum()
print(f"\n>> Gene origin:")
print(f"   Known gene-associated : {n_known:,} ({n_known/len(df_cat)*100:.1f}%)")
print(f"   Novel gene-associated : {n_novel:,} ({n_novel/len(df_cat)*100:.1f}%)")

# ──────────────────────────────────────────────
# 5. Counts by category x gene origin x coding
# ──────────────────────────────────────────────
df_cat["_gene_origin"] = df_cat["_is_novel_gene"].map(
    {True: "Novel gene", False: "Known gene"}
)
df_cat["_coding_clean"] = df_cat["coding"].str.strip()

summary = (
    df_cat.groupby(["structural_category", "_gene_origin", "_coding_clean"])
    .size()
    .reset_index(name="count")
)
print("\n>> Detailed summary (category x gene origin x coding):")
print(summary.to_string(index=False))

# ──────────────────────────────────────────────
# 6. Write tables
# ──────────────────────────────────────────────
drop_cols = ["_cat_lower", "_is_novel_gene", "_gene_origin", "_coding_clean"]
df_cat_clean = df_cat.drop(columns=[c for c in drop_cols if c in df_cat.columns])

# table 1: known gene-associated isoforms
df_known = df_cat_clean[~df_cat["_is_novel_gene"]].copy()
df_known.to_csv(OUT_ALL, index=False, sep="\t")
print(f"\n>> [1] Known gene-associated isoforms")
print(f"        file : {OUT_ALL}")
print(f"        rows : {len(df_known):,}")
print(df_known["structural_category"].value_counts().to_string())

# table 2: known gene and coding == "coding"
df_known_coding = df_known[df_known["coding"].str.strip() == "coding"].copy()
df_known_coding.to_csv(OUT_CODING, index=False, sep="\t")
print(f"\n>> [2] Known gene + coding isoforms")
print(f"        file : {OUT_CODING}")
print(f"        rows : {len(df_known_coding):,}")
print(df_known_coding["structural_category"].value_counts().to_string())

print("\n>> All done.")