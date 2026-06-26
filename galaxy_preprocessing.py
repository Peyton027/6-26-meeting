import pandas as pd

raw = pd.read_csv("RNA-seq_raw.csv")

normal_samples = ["N1", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"]
hf_samples = ["HF1", "HF2", "HF3", "HF4", "HF5", "HF6", "HF7", "HF8", "HF9", "HF10"]

samples = normal_samples + hf_samples

target_genes = [
    "CACNA1C", "CACNA1D", "CACNA1S", "CACNA1F", "CACNA1A",
    "CACNA1B", "CACNA1E", "CACNA1G", "CACNA1H", "CACNA1I",
    "RYR1", "RYR2", "RYR3", "CASQ2", "TRDN", "ASPH", "FKBP1B",
    "ATP2A1", "ATP2A2", "ATP2A3", "PLN", "SLN",
    "ITPR1", "ITPR2", "ITPR3",
    "SLC8A1", "SLC8A2", "SLC8A3", "ATP2B1", "ATP2B4",
    "MCU", "MCUB", "MICU1", "MICU2", "MICU3", "SMDT1", "MCUR1",
    "SLC8B1", "LETM1",
    "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
    "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2",
    "TRPC1", "TRPC3", "TRPC4", "TRPC6",
    "TRPM3", "TRPM4", "TRPM7", "TRPV1", "TRPV2", "PKD2",
    "CALM1", "CALM2", "CALM3",
    "KCNA4", "KCNA5", "KCNB1", "KCND2", "KCND3",
    "KCNQ1", "KCNH2", "KCNE1", "KCNE2", "KCNE3",
    "KCNAB1", "KCNAB2", "KCNIP2", "DPP6", "DPP10",
    "KCNJ2", "KCNJ4", "KCNJ12", "KCNJ14", "KCNJ3",
    "KCNJ5", "KCNJ8", "KCNJ11", "ABCC8", "ABCC9",
    "KCNK1", "KCNK2", "KCNK3", "KCNK6",
    "KCNN1", "KCNN2", "KCNN3", "KCNN4", "KCNMA1", "KCNMB1",
    "HCN1", "HCN2", "HCN4",
    "ATP1A1", "ATP1A2", "ATP1A3", "ATP1B1", "FXYD1",
    "TRPM6", "MAGT1", "CNNM1", "CNNM2", "CNNM3", "CNNM4",
    "NIPA1", "NIPA2", "NIPAL1", "NIPAL2", "NIPAL3", "NIPAL4",
    "SLC41A1", "SLC41A2", "MRS2", "SLC41A3", "MFN2",
    "CLDN16", "CLDN19", "TMEM94",
    "SCN5A", "SCN10A", "SCN8A", "SCN1A", "SCN2A",
    "SCN3A", "SCN4A", "SCN9A", "SCN11A",
    "SCN1B", "SCN2B", "SCN3B", "SCN4B",
    "SLC9A1", "SLC9A6", "SLC9A7", "SLC9A8",
    "SLC4A4", "SLC4A5", "SLC4A7",
    "CLCN1", "CLCN2", "CLCN3", "CLCN4", "CLCN5",
    "CLCN6", "CLCN7", "CLCNKA", "CLCNKB",
    "ANO1", "ANO2", "ANO6", "ANO10",
    "CFTR",
    "SLC4A1", "SLC4A2", "SLC4A3",
    "SLC12A4", "SLC12A5", "SLC12A6", "SLC12A7",
    "SLC12A1", "SLC12A2"
]

selected = raw[raw["gene_name"].isin(target_genes)].copy()

selected["gene_name"] = pd.Categorical(
    selected["gene_name"],
    categories=target_genes,
    ordered=True
)

selected = selected.sort_values("gene_name")

counts_matrix = selected[["gene_id"] + samples]

counts_matrix.to_csv(
    "Galaxy_DESeq2_Normal_vs_HF_counts_matrix.tsv",
    sep="\t",
    index=False
)

metadata = pd.DataFrame({
    "sample_id": samples,
    "condition": ["Normal"] * 9 + ["HF"] * 10
})

metadata.to_csv(
    "Galaxy_DESeq2_Normal_vs_HF_sample_metadata.tsv",
    sep="\t",
    index=False
)

annotation = selected[[
    "gene_id",
    "gene_name",
    "gene_description",
    "gene_biotype"
]]

annotation.to_csv(
    "Normal_vs_HF_target_gene_annotation.tsv",
    sep="\t",
    index=False
)