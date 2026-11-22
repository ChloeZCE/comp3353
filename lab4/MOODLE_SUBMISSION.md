# COMP3353 Lab 4 — Moodle checkpoints (code + answers)

Run the code blocks in order in one notebook/session; each block builds on the previous one.

---

## Question 1

**Code**
```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read_10x_mtx("raw_feature_bc_matrix/")

print(adata)
print(type(adata.X), adata.X.shape)
print(adata.obs.head())
print(adata.var.head())
print("Unique barcodes:", adata.obs_names.nunique())
print("Genes:", adata.n_vars)
```

**Result**
```
AnnData object with n_obs × n_vars = 6794880 × 33538
    var: 'gene_ids', 'feature_types'
Unique barcodes: 6794880
Genes: 33538
```

**Answer**
The raw matrix contains **6,794,880 unique barcodes (droplets)** and **33,538 genes**. `.X` is a sparse (CSR) count matrix of shape 6,794,880 × 33,538. `.obs` is indexed by cell barcode. `.var` is indexed by gene symbol and also has the columns `gene_ids` (Ensembl IDs) and `feature_types` ("Gene Expression").

---

## Question 2

**Code**
```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)
adata.obs.head(10)

tc = adata.obs["total_counts"]
print("Mitochondrial genes:", adata.var["mt"].sum())
print("Non-zero total_counts:", (tc > 0).sum())
print("total_counts > 100:", (tc > 100).sum())
```

**Result**
```
Mitochondrial genes: 13
Non-zero total_counts: 234600
total_counts > 100: 1365
```

**Answer**
**234,600** barcodes have non-zero `total_counts`, and **1,365** barcodes have more than 100 UMIs.

Most of the 6.8 million barcodes are whitelist barcodes that never matched a droplet, so their count is 0. The ~234,600 barcodes with non-zero counts are mostly **empty droplets**. Each one picked up a few molecules of ambient (free-floating) RNA, or got reads from sequencing or barcode errors. The 1,365 barcodes with more than 100 UMIs are droplets that most likely contain a **real cell**. This number is close to the ~1,000 cells loaded in the experiment. Some of these may still be damaged or dying cells, or doublets, which the next QC steps remove.

---

## Question 3

**Code**
```python
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
plt.tight_layout()
plt.show()

print("Barcodes kept:", (adata.obs["total_counts"] >= UMI_THR).sum())
```

**Result**
```
Barcodes kept: 1241
```
Plot: `q3_barcode_rank.png`

**Answer**
The curve stays high (above ~4,000 UMIs) until about rank 1,000. Between ranks ~1,100 and ~1,300 it drops almost vertically, from ~2,700 UMIs to ~220. After that it flattens into a background plateau of about 10–30 UMIs, which is the empty droplets. I chose **UMI_THR = 500** because it sits in the middle of the steep drop, well above the plateau. **1,241 barcodes** pass this threshold.

---

## Question 4

**Code**
```python
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
plt.tight_layout()
plt.show()

print("Cells kept:", (mt_ratio < MT_THR).sum(), "/", len(mt_ratio))
```

**Result**
```
Cells kept: 1130 / 1241
```
Plot: `q4_mt_hist.png`

**Answer**
Most cells form a single peak at about **6–14% mitochondrial reads**, which is normal for PBMCs. Above ~16–18% the distribution turns into a long, sparse tail that reaches over 90%. These high-MT cells are likely damaged or apoptotic. I chose **MT_THR = 20%** because it cuts right after the main peak. It removes the tail and keeps **1,130 of 1,241 cells**.

---

## Question 5

**Code**
```python
UMI_THR = 500
MT_THR = 20

filtered_adata = adata[adata.obs["total_counts"] >= UMI_THR].copy()
filtered_adata = filtered_adata[filtered_adata.obs["pct_counts_mt"] < MT_THR].copy()

# Exclude genes that are rarely observed
sc.pp.filter_genes(filtered_adata, min_cells=10)
print("Filtered dimensions:", filtered_adata.shape)

# Identify top HVGs
filtered_adata.layers["counts"] = filtered_adata.X.copy()
sc.pp.highly_variable_genes(
    filtered_adata,
    layer="counts",
    flavor="seurat_v3",
    n_top_genes=2000
)

# Normalization
sc.pp.normalize_total(filtered_adata, target_sum=1e4)
sc.pp.log1p(filtered_adata)
filtered_adata.layers["log"] = filtered_adata.X.copy()
filtered_adata.raw = filtered_adata

# PCA
sc.pp.pca(filtered_adata, n_comps=40, mask_var="highly_variable")

# k-NN graph & UMAP
sc.pp.neighbors(filtered_adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(filtered_adata, min_dist=0.5, spread=1.0)

pcs = filtered_adata.obsm["X_pca"]
umap_embs = filtered_adata.obsm["X_umap"]
vr = filtered_adata.uns["pca"]["variance_ratio"]

# PCA scatter
plt.figure(figsize=(4.5, 4.5))
plt.scatter(pcs[:, 0], pcs[:, 1], s=4, c="#2b5c8f", alpha=0.7, linewidths=0)
plt.xlabel(f"PC1 ({vr[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({vr[1]*100:.1f}% variance)")
plt.title(f"PCA of PBMC cells (n={filtered_adata.n_obs})")
plt.tight_layout()
plt.show()

# UMAP scatter
plt.figure(figsize=(4.5, 4.5))
plt.scatter(umap_embs[:, 0], umap_embs[:, 1], s=4, c="#2b5c8f", alpha=0.7, linewidths=0)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title(f"UMAP of PBMC cells (n={filtered_adata.n_obs})")
plt.tight_layout()
plt.show()
```
(`flavor="seurat_v3"` needs the `scikit-misc` package: `conda install -c conda-forge scikit-misc`.)

**Result**
```
Filtered dimensions: (1130, 12519)
```
Plots: `q5_pca.png`, `q5_umap.png`

**Answer**
After filtering, the final matrix is **1,130 cells × 12,519 genes**.

- **PCA** is a *linear* dimensionality-reduction method. It projects the 2,000 highly variable genes onto a small number of orthogonal axes (principal components) that capture the most variance. This compresses and denoises the data: the top 40 PCs keep the main biological signal and drop most technical noise. Those PCs are then used to build the k-nearest-neighbour graph, which both UMAP and Leiden clustering rely on.
- **UMAP** is a *non-linear* embedding used for visualisation. It builds the k-NN graph in PC space and lays the cells out in 2-D so that each cell's neighbours stay close together. This makes the different cell populations easy to see as separate groups. However, the distances between clusters and the sizes of clusters on a UMAP should not be read literally.

---

## Question 6

**Code**
```python
for res in [0.2, 0.5, 0.8, 1.2]:
    sc.tl.leiden(filtered_adata, resolution=res, key_added=f"leiden_res_{res}")
    print(res, filtered_adata.obs[f"leiden_res_{res}"].nunique(), "clusters")

sc.pl.umap(
    filtered_adata,
    color=["leiden_res_0.2", "leiden_res_0.5", "leiden_res_0.8", "leiden_res_1.2"],
    legend_loc="on data",
    legend_fontsize=10,
    frameon=False,
    ncols=2
)

# Chosen resolution
LEIDEN_RES = 0.8
sc.tl.leiden(filtered_adata, resolution=LEIDEN_RES, key_added="leiden_clusters")
sc.pl.umap(
    filtered_adata,
    color="leiden_clusters",
    legend_loc="on data",
    legend_fontsize=10,
    frameon=False,
    title=f"Leiden clusters (resolution = {LEIDEN_RES})"
)
```

**Result**
```
0.2 7 clusters
0.5 12 clusters
0.8 13 clusters
1.2 15 clusters
```
Plot: `q6_leiden_chosen.png`

**Answer**
I chose **resolution 0.8 (13 clusters)**.
- At **0.2** (7 clusters) the data is under-clustered: CD4 and CD8 T cells and NK cells are merged into a few large groups.
- From **0.5 to 0.8** the clustering separates the subsets you would expect in PBMCs: naive vs memory T cells, cytotoxic CD8 T cells vs NK cells, naive vs memory B cells, and CD14⁺ vs FCGR3A⁺ monocytes. Each new cluster has its own marker genes (e.g. GZMK/KLRB1, BANK1, FCGR3A/MS4A7). The rare populations (platelets, plasmacytoid DCs, plasma cells) also form their own clusters.
- At **1.2** (15 clusters) the large CD14⁺ monocyte and T-cell clusters are split further into sub-clusters with no clear distinct markers, which is over-clustering.

So 0.8 gives the best balance: it resolves the biologically meaningful cell types without splitting clusters that are really one population.

---

## Question 7

**Code**
```python
LEIDEN_RES = 0.8
sc.tl.leiden(filtered_adata, resolution=LEIDEN_RES, key_added="leiden_clusters")

sc.tl.rank_genes_groups(
    filtered_adata,
    groupby="leiden_clusters",
    reference="rest",
    method="wilcoxon",
    corr_method="benjamini-hochberg",
    layer="log",
    use_raw=False,
    key_added="all_clusters",
)

markers = sc.get.rank_genes_groups_df(
    filtered_adata,
    group=None,
    key="all_clusters",
)

# Remove mitochondrial genes and housekeeping ribosomal genes
noise_prefixes = ("RPS", "RPL", "MT-", "RP")
markers_filtered = markers[
    ~markers["names"].str.startswith(noise_prefixes)
].copy()
markers_filtered.head(10)

# ---- Volcano plot for cluster 3 ----
CLUSTER_OF_INTEREST = "3"
LFC_THR, PADJ_THR = 1.0, 0.05

de = markers_filtered[markers_filtered["group"] == CLUSTER_OF_INTEREST].copy()
de["neglog10_padj"] = -np.log10(de["pvals_adj"].clip(lower=1e-300))
# |LFC| > 15 are artefacts from genes with ~zero expression in one group
de = de[de["logfoldchanges"].abs() <= 15]

up = (de["logfoldchanges"] > LFC_THR) & (de["pvals_adj"] < PADJ_THR)
down = (de["logfoldchanges"] < -LFC_THR) & (de["pvals_adj"] < PADJ_THR)
ns = ~(up | down)

plt.figure(figsize=(5.5, 5))
plt.scatter(de.loc[ns, "logfoldchanges"], de.loc[ns, "neglog10_padj"],
            s=4, c="lightgrey", label="Not significant")
plt.scatter(de.loc[up, "logfoldchanges"], de.loc[up, "neglog10_padj"],
            s=6, c="#d62728", label=f"Up ({up.sum()})")
plt.scatter(de.loc[down, "logfoldchanges"], de.loc[down, "neglog10_padj"],
            s=6, c="#1f77b4", label=f"Down ({down.sum()})")
plt.axhline(-np.log10(PADJ_THR), color="k", ls="--", lw=0.8)
plt.axvline(LFC_THR, color="k", ls="--", lw=0.8)
plt.axvline(-LFC_THR, color="k", ls="--", lw=0.8)
for _, r in de[up].sort_values("neglog10_padj", ascending=False).head(10).iterrows():
    plt.annotate(r["names"], (r["logfoldchanges"], r["neglog10_padj"]),
                 fontsize=7, xytext=(3, 2), textcoords="offset points")
plt.xlabel("log2 fold change (cluster vs rest)")
plt.ylabel("-log10(adjusted p-value)")
plt.title(f"Volcano plot: cluster {CLUSTER_OF_INTEREST} vs rest")
plt.legend(fontsize=7, loc="upper left")
plt.tight_layout()
plt.show()

# Top up-regulated genes
print(de[up].sort_values("scores", ascending=False).head(15)[
    ["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False))
```

**Result**
```
    names    scores  logfoldchanges    pvals_adj
     IGHD 18.287504        8.650855 1.302868e-70
     CD37 18.219898        2.750926 2.246017e-70
     IGHM 18.168840        8.222034 3.801705e-70
    CD79A 17.909733        6.609355 3.098766e-68
    MS4A1 17.807476        6.073802 1.548270e-67
    CD79B 17.288977        5.283802 1.188371e-63
     CD74 17.199165        4.030504 4.817785e-63
LINC00926 16.704002        5.962920 1.917935e-59
    TCL1A 16.578833        9.944152 1.378989e-58
 HLA-DQA1 15.899804        4.119612 7.957448e-54
```
Plot: `q7_volcano_cluster3.png` (259 genes up, 857 genes down at |log2FC| > 1 and adjusted p < 0.05)

**Answer**
For my cluster of interest (cluster 3), the top up-regulated genes are **IGHD, IGHM, CD79A, CD79B, MS4A1 (CD20), TCL1A, CD37, CD74, LINC00926, BANK1, CD22** and several MHC class II genes (HLA-DQA1, HLA-DQB1, HLA-DRB1, HLA-DRA).

**My guess for the cell type is naive B cells.**
- CD79A/CD79B (the signalling chains of the B-cell receptor), MS4A1 (CD20) and CD22 are general B-cell markers.
- High IGHD and IGHM together with TCL1A mean the cells express unswitched IgD⁺IgM⁺ antibodies, which is typical of *naive* B cells rather than memory B cells.
- The MHC class II genes fit the B cell's role as an antigen-presenting cell.

---

## Question 8

**Code**
```python
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

sc.pl.dotplot(
    filtered_adata,
    var_names=marker_genes,
    groupby="leiden_clusters",
    standard_scale="var",
)

sc.pl.umap(
    filtered_adata,
    color=["CD3D", "MS4A1", "CD14", "NKG7", "PPBP", "HLA-DPA1"],
    ncols=3,
    frameon=False,
    layer="log",
)

# Top DEGs per cluster
for g, df in markers_filtered.groupby("group", observed=True):
    top = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 0)].head(8)
    print(f"cluster {g}: {', '.join(top['names'])}")

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

filtered_adata.obs["cell_type"] = (
    filtered_adata.obs["leiden_clusters"].map(cluster_labels).astype("category")
)

sc.pl.umap(
    filtered_adata,
    color="cell_type",
    legend_loc="on data",
    legend_fontsize=9,
    legend_fontoutline=2,
    frameon=False,
    title="Cell Type Annotations",
)
```

**Result (top DEGs per cluster)**
```
cluster 0: S100A8, S100A9, FOS, LYZ, VCAN, CTSS, DUSP1, FCN1
cluster 1: TPT1, EEF1A1, LDHB, EEF1B2, TCF7, NPM1, SARAF, PIK3IP1
cluster 2: LDHB, TRAC, IL7R, IL32, LTB, CD3E, CD3D, CD2
cluster 3: IGHD, CD37, IGHM, CD79A, MS4A1, CD79B, CD74, LINC00926
cluster 4: GZMK, KLRB1, KLRG1, GZMA, NKG7, CTSW, IL32, CCL5
cluster 5: KLRD1, GNLY, NKG7, CTSW, KLRF1, CD7, GZMA, PRF1
cluster 6: CCL5, GZMM, NKG7, B2M, CTSW, IL32, CD8A, GZMA
cluster 7: BANK1, MS4A1, IGKC, CD79A, HLA-DQA1, CD79B, CD37, HLA-DPB1
cluster 8: LST1, AIF1, FCGR3A, COTL1, SMIM25, MS4A7, HES4, SAT1
cluster 9: CST3, HLA-DRB1, HLA-DRB5, CD74, HLA-DPA1, HLA-DPB1, HLA-DRA, S100A10
cluster 10: PPBP, TAGLN2, PF4, HIST1H2AC, PGRMC1, NRGN, CDKN2D, OAZ1
cluster 11: LILRA4, IL3RA, PPP1R14B, TCF4, SERPINF1, UGCG, IRF7, SMPD3
cluster 12: IGKC, SUB1, JCHAIN, TNFRSF17, MZB1, PPIB, ITM2C, SSR4
```
Plot: `q8_celltype_umap.png`

**Answer — cluster-to-cell-type mapping**
| Cluster | Cell type | Evidence | # cells |
|---|---|---|---|
| 0 | CD14+ Monocytes | S100A8/9, LYZ, VCAN, FCN1, CD14 | 320 |
| 1 | Naive T cells (CD4/CD8 mixed) | CD3D, CCR7, TCF7, LDHB; some CD8A/B | 179 |
| 2 | CD4+ Memory T cells | CD3D/E, IL7R, IL32, LTB; no CD8 | 177 |
| 3 | Naive B cells | IGHD, IGHM, TCL1A, MS4A1, CD79A | 133 |
| 4 | GZMK+ CD8+ / MAIT T cells | CD3D, GZMK, KLRB1, KLRG1, CCL5 | 79 |
| 5 | NK cells | GNLY, NKG7, KLRD1, KLRF1, PRF1; no CD3D | 56 |
| 6 | CD8+ Cytotoxic T cells | CD3D, CD8A/B, CCL5, GZMM, GZMA | 52 |
| 7 | Memory B cells | MS4A1, CD79A, BANK1, IGKC; low IGHD/TCL1A | 49 |
| 8 | FCGR3A+ Monocytes | FCGR3A (CD16), MS4A7, LST1, AIF1 | 28 |
| 9 | Conventional DCs | FCER1A, CST3, very high HLA class II | 28 |
| 10 | Platelets | PPBP, PF4, TUBB1, CAVIN2 | 13 |
| 11 | Plasmacytoid DCs | LILRA4, IL3RA, TCF4, IRF7 | 8 |
| 12 | Plasma cells | JCHAIN, MZB1, TNFRSF17, IGKC | 8 |

**Which clusters were challenging, and why**
- **Cluster 1:** most of its top DEGs are housekeeping or translation genes (TPT1, EEF1A1, EEF1B2, NPM1) plus naive/quiescence markers (CCR7, TCF7, LDHB). These show the cells are naive T cells, but they don't separate CD4 from CD8. CD8A/CD8B are expressed at an intermediate level, so the cluster probably mixes naive CD4 and naive CD8 T cells. I gave it the broader label "Naive T cells (CD4/CD8 mixed)".
- **Cluster 4:** it expresses CD3D, GZMK, KLRB1 (CD161), KLRG1 and IL7R. This pattern fits both GZMK⁺ CD8 effector-memory T cells and MAIT cells. Telling them apart would need extra markers (e.g. TRAV1-2, SLC4A10) or TCR data, so I used a combined label.
- **Clusters 5 and 6 (NK vs CD8 cytotoxic T cells):** both share cytotoxic genes (NKG7, GZMA, PRF1, CCL5). I separated them using CD3D/CD8A, which are present in cluster 6 and absent in cluster 5, and GNLY/KLRF1, which are much higher in cluster 5.
- **Small clusters 9–12:** these have only 8–28 cells each, so their DE statistics are weaker. Even so, each has very specific markers (FCER1A, LILRA4, JCHAIN/MZB1, PPBP). The cDC cluster also shows some monocyte genes (LYZ, S100A9), probably because cDCs and monocytes share myeloid genes, or because a few doublets are present.
