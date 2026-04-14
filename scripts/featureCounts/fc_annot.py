import pandas as pd
import re
import os
import glob

# ── CONFIG ────────────────────────────────────────────────────────────────────
GTF = "TE_hg_clean.gtf"

FC_DIRS = {
    "DNA":    "DNA",
    "RNA_s0": "RNA_s0",
}
# ──────────────────────────────────────────────────────────────────────────────

# 1. Build gene_id -> label mapping from GTF
print("Parsing GTF...")
gene_label = {}
with open(GTF) as f:
    for line in f:
        if line.startswith("#") or "\texon\t" not in line:
            continue
        
        chrom = line.split("\t")[0]
        gid = re.search(r'gene_id "([^"]+)"', line)
        cid = re.search(r'class_id "([^"]+)"', line)
        gt  = re.search(r'gene_type "([^"]+)"', line)

        if gid:
            if cid:
                label = cid.group(1)
            elif gt:
                label = gt.group(1)
                # Reclassify chrM protein-coding genes
                if chrom == "chrM" and label == "protein_coding":
                    label = "Mt_protein_coding"
            else:
                label = "unknown"
            gene_label[gid.group(1)] = label

print(f"  {len(gene_label)} unique gene_ids mapped\n")

# 2. Process each fc.txt file
records = []

for condition, directory in FC_DIRS.items():
    fc_files = glob.glob(os.path.join(directory, "*.txt"))
    fc_files = [f for f in fc_files if not f.endswith(".summary")]

    for fc_path in sorted(fc_files):
        sample = os.path.basename(fc_path).replace("_fc.txt", "")
        summary_path = fc_path + ".summary"

        # ── Parse .summary for alignment stats ────────────────────────────
        total_alignments = None
        assigned = None

        if os.path.exists(summary_path):
            status_counts = {}
            with open(summary_path) as f:
                next(f)  # skip header line ("Status\t/path/to.bam")
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        status_counts[parts[0]] = int(parts[1])

            assigned         = status_counts.get("Assigned", 0)
            total_alignments = sum(status_counts.values())  # sum all rows

        # ── Parse fc.txt counts ────────────────────────────────────────────
        fc = pd.read_csv(fc_path, sep="\t", comment="#")
        fc.columns = [c.strip() for c in fc.columns]

        gene_col  = fc.columns[0]   # Geneid
        count_col = fc.columns[-1]  # last column = BAM counts

        fc = fc[[gene_col, count_col]].copy()
        fc.columns = ["gene_id", "count"]

        fc["label"] = fc["gene_id"].map(gene_label).fillna("unknown")

        label_counts = fc.groupby("label")["count"].sum()

        row = {
            "sample":              sample,
            "condition":           condition,
            "total_alignments":    total_alignments,
            "assigned_alignments": assigned,
            "pct_assigned":        round(100 * assigned / total_alignments, 2)
                                   if total_alignments and assigned else None,
        }
        row.update(label_counts.to_dict())
        records.append(row)

# 3. Build final table
df = pd.DataFrame(records)

meta_cols  = ["sample", "condition", "total_alignments", "assigned_alignments", "pct_assigned"]
label_cols = sorted([c for c in df.columns if c not in meta_cols])
df = df[meta_cols + label_cols].fillna(0)

for col in label_cols:
    df[col] = df[col].astype(int)

df.to_csv("featurecounts_summary.csv", index=False)
print("Saved: featurecounts_summary.csv")
print(df[meta_cols + label_cols[:5]].to_string(index=False))
