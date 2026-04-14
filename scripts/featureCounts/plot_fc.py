import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

df = pd.read_csv("featurecounts_summary.csv")

# ── CONFIG ────────────────────────────────────────────────────────────────────
meta_cols = ["sample", "condition", "total_alignments", "assigned_alignments", "pct_assigned"]
label_cols = [c for c in df.columns if c not in meta_cols]

GROUP_MAP = {
    "protein_coding":                     "protein_coding",
    "Mt_protein_coding":                  "Mt_protein_coding",
    "Mt_rRNA":                            "Mt_noncoding",
    "Mt_tRNA":                            "Mt_noncoding",
    "lncRNA":                             "noncoding_RNA",
    "miRNA":                              "noncoding_RNA",
    "misc_RNA":                           "noncoding_RNA",
    "rRNA":                               "noncoding_RNA",
    "rRNA_pseudogene":                    "noncoding_RNA",
    "scaRNA":                             "noncoding_RNA",
    "scRNA":                              "noncoding_RNA",
    "snoRNA":                             "noncoding_RNA",
    "snRNA":                              "noncoding_RNA",
    "sRNA":                               "noncoding_RNA",
    "vault_RNA":                          "noncoding_RNA",
    "ribozyme":                           "noncoding_RNA",
    "SINE":                               "SINEs",
    "SINE?":                              "SINEs",
    "LINE":                               "LINEs",
    "LTR":                                "LTRs",
    "LTR?":                               "LTRs",
    "DNA":                                "Other_repetitive",
    "DNA?":                               "Other_repetitive",
    "RC":                                 "Other_repetitive",
    "RC?":                                "Other_repetitive",
    "Retroposon":                         "Other_repetitive",
    "RNA":                                "Other_repetitive",
    "Satellite":                          "Other_repetitive",
    "Unknown":                            "Other_repetitive",
    "processed_pseudogene":               "Other_assigned",
    "unprocessed_pseudogene":             "Other_assigned",
    "transcribed_processed_pseudogene":   "Other_assigned",
    "transcribed_unprocessed_pseudogene": "Other_assigned",
    "transcribed_unitary_pseudogene":     "Other_assigned",
    "translated_processed_pseudogene":    "Other_assigned",
    "translated_unprocessed_pseudogene":  "Other_assigned",
    "unitary_pseudogene":                 "Other_assigned",
    "polymorphic_pseudogene":             "Other_assigned",
    "pseudogene":                         "Other_assigned",
    "TEC":                                "Other_assigned",
    "IG_C_gene":                          "Other_assigned",
    "IG_C_pseudogene":                    "Other_assigned",
    "IG_D_gene":                          "Other_assigned",
    "IG_J_gene":                          "Other_assigned",
    "IG_J_pseudogene":                    "Other_assigned",
    "IG_pseudogene":                      "Other_assigned",
    "IG_V_gene":                          "Other_assigned",
    "IG_V_pseudogene":                    "Other_assigned",
    "TR_C_gene":                          "Other_assigned",
    "TR_D_gene":                          "Other_assigned",
    "TR_J_gene":                          "Other_assigned",
    "TR_J_pseudogene":                    "Other_assigned",
    "TR_V_gene":                          "Other_assigned",
    "TR_V_pseudogene":                    "Other_assigned",
}

COLORS = {
    "protein_coding":      "#6BB5D6",
    "Mt_protein_coding":   "#1A7A2E",
    "Mt_noncoding":        "#4ABFBF",
    "noncoding_RNA":       "#ADD8E6",
    "SINEs":               "#8B9A2D",
    "LINEs":               "#D4C97A",
    "LTRs":                "#E8897A",
    "Other_repetitive":    "#7B2D5E",
    "Other_assigned":      "#D46BB5",
    "intronic_intergenic": "#2D2D8B",
}

STACK_ORDER = [
    "protein_coding", "Mt_protein_coding", "Mt_noncoding", "noncoding_RNA",
    "SINEs", "LINEs", "LTRs", "Other_repetitive", "Other_assigned",
    "intronic_intergenic",
]

legend_labels = {
    "protein_coding":      "Protein Coding",
    "Mt_protein_coding":   "Mitochondrial Protein Coding",
    "Mt_noncoding":        "Mitochondrial Noncoding",
    "noncoding_RNA":       "Noncoding RNA",
    "SINEs":               "SINEs",
    "LINEs":               "LINEs",
    "LTRs":                "LTRs",
    "Other_repetitive":    "Other Repetitive",
    "Other_assigned":      "Other Assigned",
    "intronic_intergenic": "Intergenic/Intronic",
}

# ── PREPARE DATA ──────────────────────────────────────────────────────────────
plot_df = df.copy()

group_data = {}
for col in label_cols:
    group = GROUP_MAP.get(col, col)
    group_data[group] = group_data.get(group, 0) + plot_df[col].values

group_df = pd.DataFrame(group_data, index=plot_df.index)
group_df["intronic_intergenic"] = plot_df["total_alignments"] - plot_df["assigned_alignments"]
group_df["sample"]    = plot_df["sample"]
group_df["condition"] = plot_df["condition"]
group_df["total"]     = plot_df["total_alignments"]

biotype_cols = [c for c in group_df.columns if c not in ["sample", "condition", "total"]]
for col in biotype_cols:
    group_df[col] = group_df[col] / group_df["total"] * 100

STACK_ORDER = [c for c in STACK_ORDER if c in group_df.columns]

gdf = group_df.set_index("sample")

def avg_samples(names):
    return gdf.loc[names, biotype_cols].mean()

# ── PLOT FUNCTION ─────────────────────────────────────────────────────────────
def make_plot(rna_condition, dna_condition, output_prefix):

    # Averaged pairs only: RNA first then DNA, labeled by individual
    bars = [
        (avg_samples([f"120916_R_{rna_condition}", f"40TL_R_{rna_condition}"]), "GLAHM:120916", "RNA"),
        (avg_samples([f"120916_D_{dna_condition}", f"44TL_D_{dna_condition}"]), "GLAHM:120916", "DNA"),
        (avg_samples([f"120917_R_{rna_condition}", f"48TL_R_{rna_condition}"]), "GLAHM:120917", "RNA"),
        (avg_samples([f"120917_D_{dna_condition}", f"45TL_D_{dna_condition}"]), "GLAHM:120917", "DNA"),
    ]

    x_positions = np.arange(len(bars))

    fig, ax = plt.subplots(figsize=(6, 4))

    for i, (props, label, condition) in enumerate(bars):
        x = x_positions[i]
        bottom = 0
        width = 0.7

        for biotype in STACK_ORDER:
            val = props.get(biotype, 0)
            ax.bar(x, val, bottom=bottom,
                   color=COLORS.get(biotype, "#AAAAAA"),
                   width=width,
                   edgecolor="none")
            bottom += val

        # Condition label above bar — black, bold, no color coding
        ax.text(x, 101.5, condition, ha="center", va="bottom",
                fontsize=8, fontweight="bold", color="black")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(["", "", "", ""], fontsize=9)

    pairs = [("GLAHM:120916", 0, 1), ("GLAHM:120917", 2, 3)]
    for label, i1, i2 in pairs:
        center = (x_positions[i1] + x_positions[i2]) / 2
        ax.text(center, -6, label, ha="center", va="top",
                fontsize=9, fontweight="bold", transform=ax.transData)
    ax.set_ylim(0, 108)
    ax.set_ylabel("Percentage (%) of Mapped Reads", fontsize=11)
    ax.set_xlim(x_positions[0] - 0.6, x_positions[-1] + 0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        mpatches.Patch(color=COLORS[b], label=legend_labels[b])
        for b in STACK_ORDER if b in COLORS
    ]
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=9, frameon=False, title_fontsize=10)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_prefix}.pdf / .png")
    plt.show()


# ── GENERATE BOTH PLOTS ───────────────────────────────────────────────────────
# Plot 1: RNA (s0-stranded) + DNA (unstranded)
make_plot(rna_condition="RNA_s0", dna_condition="DNA",
          output_prefix="biotype_composition_RNA_s0_DNA")