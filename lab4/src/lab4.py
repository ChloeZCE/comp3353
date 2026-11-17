"""COMP3353 Lab Session 4 - scRNA-seq analysis of the 10x PBMC 1k v3 dataset.

Run from the lab4/ folder:  python3 src/lab4.py
Plots are written to results/, text output to stdout.
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)

OUT = "results"
os.makedirs(OUT, exist_ok=True)
sc.settings.figdir = OUT
sc.settings.verbosity = 0


def save(name):
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, name), dpi=200, bbox_inches="tight")
    plt.close("all")


# ---------------------------------------------------------------- Question 1
adata = sc.read_10x_mtx("data/raw_feature_bc_matrix/")
print("== Q1 ==")
print(adata)
print("X:", type(adata.X), adata.X.shape)
print("obs (barcodes):", adata.obs_names[:3].tolist(), "...")
print(adata.var.head())
print(f"Unique barcodes: {adata.obs_names.nunique()}  Genes: {adata.n_vars}")

# ---------------------------------------------------------------- Question 2
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
print("\n== Q2 ==")
print(adata.obs.head(10))
tc = adata.obs["total_counts"]
print("Mitochondrial genes:", int(adata.var["mt"].sum()))
print("Barcodes with total_counts > 0  :", int((tc > 0).sum()))
print("Barcodes with total_counts > 100:", int((tc > 100).sum()))

# ---------------------------------------------------------------- Question 3
UMI_THR = 500
counts = sorted(adata.obs["total_counts"].values)
ranks = list(reversed(np.arange(1, len(counts) + 1)))
plt.figure(figsize=(4, 4))
plt.loglog(ranks, counts, color="#2b5c8f", lw=2)
plt.axhline(UMI_THR, color="red", ls="--", lw=1.2, label=f"UMI threshold = {UMI_THR}")
plt.xlabel("Barcode rank (descending UMI count)")
plt.ylabel("Total UMI counts per barcode")
plt.title("Barcode rank plot (knee plot)")
plt.legend(fontsize=8)
save("q3_barcode_rank.png")
print("\n== Q3 ==")
print(f"UMI_THR = {UMI_THR}; barcodes kept: {int((tc >= UMI_THR).sum())}")

# ---------------------------------------------------------------- Question 4
MT_THR = 20
filtered_adata = adata[adata.obs["total_counts"] >= UMI_THR].copy()
mt_ratio = filtered_adata.obs["pct_counts_mt"].values
plt.figure(figsize=(4, 4))
plt.hist(mt_ratio, bins=np.arange(0, 100, 2), color="#2b5c8f")
plt.axvline(MT_THR, color="red", ls="--", lw=1.2, label=f"MT threshold = {MT_THR}%")
plt.xlabel("Mitochondrial counts (% of total UMIs)")
plt.ylabel("Number of cells")
plt.title("Mitochondrial read fraction")
plt.legend(fontsize=8)
save("q4_mt_hist.png")
print("\n== Q4 ==")
print(f"MT_THR = {MT_THR}%; cells kept: {int((mt_ratio < MT_THR).sum())} / {len(mt_ratio)}")

# ---------------------------------------------------------------- Question 5
print("\n== Q5 ==")
filtered_adata = adata[adata.obs["total_counts"] >= UMI_THR].copy()
filtered_adata = filtered_adata[filtered_adata.obs["pct_counts_mt"] < MT_THR].copy()
sc.pp.filter_genes(filtered_adata, min_cells=10)
print("Filtered dimensions:", filtered_adata.shape)

filtered_adata.layers["counts"] = filtered_adata.X.copy()
sc.pp.highly_variable_genes(
    filtered_adata, layer="counts", flavor="seurat_v3", n_top_genes=2000
)
sc.pp.normalize_total(filtered_adata, target_sum=1e4)
sc.pp.log1p(filtered_adata)
filtered_adata.layers["log"] = filtered_adata.X.copy()
filtered_adata.raw = filtered_adata

sc.pp.pca(filtered_adata, n_comps=40, mask_var="highly_variable")
sc.pp.neighbors(filtered_adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(filtered_adata, min_dist=0.5, spread=1.0)
pcs = filtered_adata.obsm["X_pca"]
umap_embs = filtered_adata.obsm["X_umap"]
vr = filtered_adata.uns["pca"]["variance_ratio"]

plt.figure(figsize=(4.5, 4.5))
plt.scatter(pcs[:, 0], pcs[:, 1], s=4, c="#2b5c8f", alpha=0.7, linewidths=0)
plt.xlabel(f"PC1 ({vr[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({vr[1]*100:.1f}% variance)")
plt.title(f"PCA of PBMC cells (n={filtered_adata.n_obs})")
save("q5_pca.png")

plt.figure(figsize=(4.5, 4.5))
plt.scatter(umap_embs[:, 0], umap_embs[:, 1], s=4, c="#2b5c8f", alpha=0.7, linewidths=0)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title(f"UMAP of PBMC cells (n={filtered_adata.n_obs})")
save("q5_umap.png")

# ---------------------------------------------------------------- Question 6
print("\n== Q6 ==")
for res in [0.2, 0.5, 0.8, 1.2]:
    sc.tl.leiden(filtered_adata, resolution=res, key_added=f"leiden_res_{res}")
    print(f"resolution {res}: {filtered_adata.obs[f'leiden_res_{res}'].nunique()} clusters")
sc.pl.umap(
    filtered_adata,
    color=["leiden_res_0.2", "leiden_res_0.5", "leiden_res_0.8", "leiden_res_1.2"],
    legend_loc="on data", legend_fontsize=10, frameon=False, ncols=2,
    show=False,
)
save("q6_leiden_resolutions.png")

LEIDEN_RES = 0.8
sc.tl.leiden(filtered_adata, resolution=LEIDEN_RES, key_added="leiden_clusters")
print(f"Chosen LEIDEN_RES = {LEIDEN_RES}")
print(filtered_adata.obs["leiden_clusters"].value_counts().sort_index().to_string())
sc.pl.umap(
    filtered_adata, color="leiden_clusters", legend_loc="on data",
    legend_fontsize=10, frameon=False,
    title=f"Leiden clusters (resolution = {LEIDEN_RES})", show=False,
)
save("q6_leiden_chosen.png")

# ---------------------------------------------------------------- Question 7
print("\n== Q7 ==")
sc.tl.rank_genes_groups(
    filtered_adata, groupby="leiden_clusters", reference="rest",
    method="wilcoxon", corr_method="benjamini-hochberg",
    layer="log", use_raw=False, key_added="all_clusters",
)
markers = sc.get.rank_genes_groups_df(filtered_adata, group=None, key="all_clusters")
noise_prefixes = ("RPS", "RPL", "MT-", "RP")
markers_filtered = markers[~markers["names"].str.startswith(noise_prefixes)].copy()
print(markers_filtered.head(10).to_string())

# Top 8 up-regulated markers per cluster (padj < 0.05, ranked by score)
print("\nTop up-regulated markers per cluster:")
for g, df in markers_filtered.groupby("group", observed=True):
    top = df[(df.pvals_adj < 0.05) & (df.logfoldchanges > 0)].head(8)
    print(f"  cluster {g}: {', '.join(top['names'])}")

CLUSTER_OF_INTEREST = os.environ.get("COI", "3")
LFC_THR, PADJ_THR = 1.0, 0.05
de = markers_filtered[markers_filtered["group"] == CLUSTER_OF_INTEREST].copy()
de["neglog10_padj"] = -np.log10(de["pvals_adj"].clip(lower=1e-300))
# |LFC| > 15 are artefacts of genes with ~zero expression in one group
# (scanpy's log fold change explodes when a group mean is 0) -> not plotted.
n_art = int((de.logfoldchanges.abs() > 15).sum())
de = de[de.logfoldchanges.abs() <= 15]
print(f"Dropped {n_art} genes with |LFC| > 15 (zero-expression artefacts) from volcano plot")
up = (de.logfoldchanges > LFC_THR) & (de.pvals_adj < PADJ_THR)
down = (de.logfoldchanges < -LFC_THR) & (de.pvals_adj < PADJ_THR)

plt.figure(figsize=(5.5, 5))
plt.scatter(de.loc[~(up | down), "logfoldchanges"], de.loc[~(up | down), "neglog10_padj"],
            s=4, c="lightgrey", label="Not significant")
plt.scatter(de.loc[up, "logfoldchanges"], de.loc[up, "neglog10_padj"],
            s=6, c="#d62728", label=f"Up ({int(up.sum())})")
plt.scatter(de.loc[down, "logfoldchanges"], de.loc[down, "neglog10_padj"],
            s=6, c="#1f77b4", label=f"Down ({int(down.sum())})")
plt.axhline(-np.log10(PADJ_THR), color="k", ls="--", lw=0.8)
plt.axvline(LFC_THR, color="k", ls="--", lw=0.8)
plt.axvline(-LFC_THR, color="k", ls="--", lw=0.8)
for _, r in de[up].sort_values("neglog10_padj", ascending=False).head(10).iterrows():
    plt.annotate(r["names"], (r.logfoldchanges, r.neglog10_padj), fontsize=7,
                 xytext=(3, 2), textcoords="offset points")
plt.xlabel("log2 fold change (cluster vs rest)")
plt.ylabel("-log10(adjusted p-value)")
plt.title(f"Volcano plot: cluster {CLUSTER_OF_INTEREST} vs rest")
plt.legend(fontsize=7, loc="upper left")
save(f"q7_volcano_cluster{CLUSTER_OF_INTEREST}.png")
print(f"\nCluster {CLUSTER_OF_INTEREST} top up-regulated genes (|LFC|>{LFC_THR}, padj<{PADJ_THR}):")
print(de[up].sort_values("scores", ascending=False).head(15)[
    ["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False))

# ---------------------------------------------------------------- Question 8
print("\n== Q8 ==")
marker_genes = {
    "CD4+ T cells": ["CD3D", "IL7R", "TCF7"],
    "CD8+ T cells": ["CD3D", "GZMA", "GZMK", "CD8A"],
    "NK cells": ["NKG7", "GNLY", "PRF1", "KLRD1"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "Monocytes": ["CD14", "S100A9", "LYZ"],
    "FCGR3A+ Mono": ["FCGR3A", "MS4A7"],
    "Dendritic cells": ["HLA-DRA", "HLA-DPA1", "CST3"],
    "pDC": ["LILRA4", "IL3RA"],
    "Platelets": ["PPBP", "TUBB1", "CAVIN2"],
}
sc.pl.dotplot(filtered_adata, var_names=marker_genes, groupby="leiden_clusters",
              standard_scale="var", show=False)
save("q8_dotplot.png")
sc.pl.umap(filtered_adata, color=["CD3D", "MS4A1", "CD14", "NKG7", "PPBP", "HLA-DPA1"],
           ncols=3, frameon=False, layer="log", show=False)
save("q8_marker_umaps.png")

# mean log expression of canonical markers per cluster (to support annotation)
genes = sorted({g for v in marker_genes.values() for g in v} | {"CD8B", "CCR7", "FCER1A"})
genes = [g for g in genes if g in filtered_adata.var_names]
expr = sc.get.obs_df(filtered_adata, keys=genes + ["leiden_clusters"], layer="log")
print(expr.groupby("leiden_clusters", observed=True).mean().round(2).T.to_string())

cluster_labels = {
    "0": "CD14+ Monocytes",
    "1": "Naive T cells (CD4/CD8 mixed)",
    "2": "CD4+ Memory T cells",
    "3": "Naive B cells",
    "4": "GZMK+ CD8+ / MAIT T cells",
    "5": "NK cells",
    "6": "CD8+ Cytotoxic T cells",
    "7": "Memory B cells",
    "8": "FCGR3A+ Monocytes",
    "9": "Conventional DCs",
    "10": "Platelets",
    "11": "Plasmacytoid DCs",
    "12": "Plasma cells",
}
if True:
    filtered_adata.obs["cell_type"] = (
        filtered_adata.obs["leiden_clusters"].map(cluster_labels).astype("category")
    )
    sc.pl.umap(filtered_adata, color="cell_type", legend_loc="on data",
               legend_fontsize=9, legend_fontoutline=2, frameon=False,
               title="Cell Type Annotations", show=False)
    save("q8_celltype_umap.png")
    print(filtered_adata.obs["cell_type"].value_counts().to_string())
