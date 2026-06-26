import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

DESEQ_FILE = "Dk N_HF.tabular"
RAW_FILE = "RNA-seq_raw.csv"
OUTPUT_FILE = "Sodium_channel_profiling_HF_vs_Normal.png"

TITLE = "Sodium channel profiling: Heart Failure vs Normal"

REFERENCE_LABEL = "Normal"
COMPARISON_LABEL = "Heart Failure"

RESULTS_HAS_HEADER = False
FLIP_LOG2FC = False

COLOR_MIN = -3
COLOR_MAX = 3

#input relevant genes 
GENE_GROUPS = {
    "Voltage-gated Na+\nchannels": [
        "SCN11A",
        "SCN8A",
        "SCN4A",
        "SCN10A",
        "SCN5A",
        "SCN1A",
        "SCN3A",
        "SCN9A",
        "SCN2A"
    ],
    "Nav beta\nsubunits": [
        "SCN2B",
        "SCN3B",
        "SCN4B",
        "SCN1B"
    ],
    "Na+/H+\nexchangers": [
        "SLC9A1",
        "SLC9A8",
        "SLC9A7",
        "SLC9A6"
    ],
    "Na+/Ca2+\nexchangers": [
        "SLC8A3",
        "SLC8A2",
        "SLC8A1",
        "SLC8B1"
    ],
    "Na+/HCO3-\nexchangers": [
        "SLC4A4",
        "SLC4A5",
        "SLC4A7"
    ],
    "Na+/K+\npump": [
        "ATP1A3",
        "ATP1A2",
        "ATP1A1",
        "ATP1B1",
        "FXYD1"
    ]
}

if RESULTS_HAS_HEADER:
    deseq = pd.read_csv(DESEQ_FILE, sep="\t")

    if "gene_id" not in deseq.columns:
        deseq = deseq.rename(columns={deseq.columns[0]: "gene_id"})

else:
    deseq = pd.read_csv(DESEQ_FILE, sep="\t", header=None)

    deseq.columns = [
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    ]

deseq["gene_id"] = (
    deseq["gene_id"]
    .astype(str)
    .str.replace(r"\..*$", "", regex=True)
)

deseq["log2FoldChange"] = pd.to_numeric(
    deseq["log2FoldChange"],
    errors="coerce"
)

deseq["padj"] = pd.to_numeric(
    deseq["padj"],
    errors="coerce"
)

deseq = deseq.dropna(subset=["log2FoldChange"])

if FLIP_LOG2FC:
    deseq["log2FoldChange"] = -deseq["log2FoldChange"]

annotation = pd.read_csv(
    RAW_FILE,
    usecols=["gene_id", "gene_name"]
)

annotation["gene_id"] = (
    annotation["gene_id"]
    .astype(str)
    .str.replace(r"\..*$", "", regex=True)
)

annotation = annotation.drop_duplicates(subset="gene_id")

deseq = deseq.merge(
    annotation,
    on="gene_id",
    how="left"
)

deseq = deseq.dropna(subset=["gene_name"])

deseq = deseq.sort_values(
    "padj",
    na_position="last"
)

deseq = deseq.drop_duplicates(
    subset="gene_name"
)

deseq_by_gene = deseq.set_index("gene_name")

def get_stars(padj):
    if pd.isna(padj) or padj >= 0.05:
        return ""
    elif padj < 0.0001:
        return "****"
    elif padj < 0.001:
        return "***"
    elif padj < 0.01:
        return "**"
    else:
        return "*"

rows = []
group_positions = []

group_items = list(GENE_GROUPS.items())

for group_number, (group_name, genes) in enumerate(group_items):
    start_row = len(rows)

    for gene in genes:
        if gene in deseq_by_gene.index:
            row = deseq_by_gene.loc[gene]

            rows.append({
                "gene": gene,
                "log2FoldChange": row["log2FoldChange"],
                "padj": row["padj"],
                "blank": False
            })

    end_row = len(rows) - 1

    if end_row >= start_row:
        middle_row = (start_row + end_row) / 2

        group_positions.append(
            (group_name, start_row, end_row, middle_row)
        )

        if group_number != len(group_items) - 1:
            rows.append({
                "gene": "",
                "log2FoldChange": np.nan,
                "padj": np.nan,
                "blank": True
            })

heatmap_values = np.array([
    [-row["log2FoldChange"], row["log2FoldChange"]]
    if not row["blank"]
    else [np.nan, np.nan]
    for row in rows
])

cmap = LinearSegmentedColormap.from_list(
    "blue_white_red",
    ["blue", "white", "red"]
)

cmap.set_bad("white")

fig_height = max(7, len(rows) * 0.23)

fig = plt.figure(
    figsize=(7.2, fig_height),
    facecolor="white"
)

gs = fig.add_gridspec(
    1,
    3,
    width_ratios=[1.45, 2.8, 0.15],
    wspace=0.04
)

ax_group = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[0, 1])
cax = fig.add_subplot(gs[0, 2])

ax_group.set_xlim(0, 1)
ax_group.set_ylim(len(rows) - 0.5, -0.5)

for group_name, start, end, middle in group_positions:
    ax_group.vlines(
        0.90,
        start - 0.5,
        end + 0.5,
        color="black",
        linewidth=0.8
    )

    ax_group.text(
        0.35,
        middle,
        group_name,
        ha="center",
        va="center",
        fontsize=7,
        fontweight="bold"
    )

ax_group.axis("off")

im = ax_heat.imshow(
    heatmap_values,
    aspect="auto",
    cmap=cmap,
    vmin=COLOR_MIN,
    vmax=COLOR_MAX,
    interpolation="none"
)

ax_heat.set_xticks([0, 1])

ax_heat.set_xticklabels(
    [REFERENCE_LABEL, COMPARISON_LABEL],
    fontsize=11
)

ax_heat.set_yticks(np.arange(len(rows)))

ax_heat.set_yticklabels(
    [row["gene"] for row in rows],
    fontsize=7,
    fontweight="bold"
)

ax_heat.tick_params(
    axis="y",
    length=0,
    pad=2
)

ax_heat.tick_params(
    axis="x",
    length=0,
    pad=6
)

for spine in ax_heat.spines.values():
    spine.set_visible(False)

for y, row in enumerate(rows):
    stars = get_stars(row["padj"])

    if stars != "":
        ax_heat.text(
            1,
            y,
            stars,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold"
        )

cbar = plt.colorbar(
    im,
    cax=cax,
    ticks=np.arange(COLOR_MIN, COLOR_MAX + 1, 1)
)

cbar.set_label(
    "log2FC",
    fontsize=11
)

cbar.ax.tick_params(
    labelsize=8,
    length=0
)

fig.suptitle(
    TITLE,
    x=0.05,
    y=0.985,
    ha="left",
    fontsize=12
)

fig.subplots_adjust(
    left=0.05,
    right=0.94,
    top=0.94,
    bottom=0.06
)

plt.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.show()