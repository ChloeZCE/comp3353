# COMP3353 Lab Session 4 — Single-cell RNA-seq (PBMC 1k v3)

All code is in `src/lab4.py` (run from `lab4/` with `python3 src/lab4.py`).
Full console output: `results/lab4_output.txt`. Plots: `results/*.png`.
Environment: scanpy 1.13, python-igraph, leidenalg, scikit-misc (needed by `flavor="seurat_v3"`).

---

## Question 1 — Raw matrix

```python
adata = sc.read_10x_mtx("data/raw_feature_bc_matrix/")
print(adata)                     # AnnData object with n_obs × n_vars = 6794880 × 33538
print(adata.obs_names.nunique(), adata.n_vars)
```

`.X` is a 6,794,880 × 33,538 sparse CSR count matrix, `.obs` is indexed by cell barcodes,
and `.var` holds gene symbols with `gene_ids` (Ensembl) and `feature_types`.

**The raw matrix has 6,794,880 unique barcodes (droplets) and 33,538 genes.**

## Question 2 — QC metrics

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
tc = adata.obs["total_counts"]
(tc > 0).sum(), (tc > 100).sum()
```

- Non-zero `total_counts`: **234,600** barcodes
- `total_counts` > 100: **1,365** barcodes
- 13 mitochondrial (`MT-`) genes were found.

Most of the 6.8 M barcodes are entries on the 10x whitelist that never matched a droplet, so they have 0 counts.
The ~234k non-zero barcodes are mostly **empty droplets**. Each picked up a few molecules of ambient (free-floating) RNA,
or got reads from sequencing or barcode errors.
The ~1,365 barcodes with more than 100 UMIs are the droplets that probably contain a **real cell**. This number is close
to the ~1,000 cells loaded in this experiment. Some of them may still be damaged or dying cells, or doublets.

## Question 3 — Barcode-rank (knee) plot

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
```

Plot: `results/q3_barcode_rank.png`

The curve stays above ~4,000 UMIs until about rank 1,000. Between ranks ~1,100 and ~1,300 it drops almost vertically
(2,659 UMIs at rank 1,100 → 220 at rank 1,300), then flattens into a background plateau at roughly 10–30 UMIs.
**UMI_THR = 500** falls in the middle of this cliff, well clear of the plateau. **1,241 barcodes** pass it.

## Question 4 — Mitochondrial fraction

```python
MT_THR = 20
filtered_adata = adata[adata.obs["total_counts"] >= UMI_THR].copy()
mt_ratio = filtered_adata.obs["pct_counts_mt"].values
plt.figure(figsize=(4, 4))
plt.hist(mt_ratio, bins=np.arange(0, 100, 2))
plt.axvline(MT_THR, color="red", ls="--", label=f"MT threshold = {MT_THR}%")
plt.xlabel("Mitochondrial counts (% of total UMIs)")
plt.ylabel("Number of cells")
plt.title("Mitochondrial read fraction")
```

Plot: `results/q4_mt_hist.png`

Most cells form one peak at **~6–14 % MT reads**, which is normal for PBMCs. After about 16–18 % the
distribution becomes a long, sparse tail that reaches over 90 %. These high-MT cells are likely damaged or apoptotic.
I chose **MT_THR = 20 %**: it cuts at the end of the main peak and keeps **1,130 of 1,241 cells**.

## Question 5 — Filtering, normalisation, PCA, UMAP

Code as given in the handout, with `UMI_THR = 500` and `MT_THR = 20`. Default parameters were used
(`n_comps=40`, `n_neighbors=15`, `n_pcs=40`, `min_dist=0.5`, `spread=1.0`).

**Final filtered dimensions: 1,130 cells × 12,519 genes.**

- **PCA** is a *linear* dimensionality reduction. It takes the 2,000 highly variable genes and finds a few orthogonal
  axes (principal components) that capture the most variance, which compresses the data and removes noise.
  The top 40 PCs keep the main biological signal and drop most technical noise. They are then used to build the
  k-NN graph, which UMAP and Leiden clustering both rely on.
- **UMAP** is a *non-linear* embedding for visualisation. It builds the k-NN graph in PC space and lays the cells out
  in 2-D so that each cell's nearest neighbours stay close together. This makes the separate cell populations easy to see.
  Distances between clusters, and cluster sizes, in UMAP space should not be read literally.

Plots: `results/q5_pca.png` (PC1 vs PC2) and `results/q5_umap.png`.

## Question 6 — Leiden resolution

| resolution | # clusters |
|---|---|
| 0.2 | 7 |
| 0.5 | 12 |
| 0.8 | 13 |
| 1.2 | 15 |

Plot of all four: `results/q6_leiden_resolutions.png`

**I chose resolution 0.8.**

- At 0.2, the major lineages are merged: CD4 and CD8 T cells and NK cells are grouped too coarsely.
- 0.5 → 0.8 separates the subsets expected in PBMCs: naive vs memory T cells, cytotoxic CD8 T vs NK cells,
  naive vs memory B cells, and CD14⁺ vs FCGR3A⁺ monocytes. Each new cluster has its own marker genes
  (e.g. GZMK/KLRB1, BANK1, FCGR3A/MS4A7).
- Rare populations (platelets, pDCs, plasma cells) also form their own clusters at 0.8.
- At 1.2, the large clusters (the CD14⁺ monocytes and T cells) are split into sub-clusters that don't have clear,
  distinct markers. This is over-clustering.

Chosen-resolution UMAP: `results/q6_leiden_chosen.png`

## Question 7 — DE and volcano plot (cluster 3)

```python
CLUSTER_OF_INTEREST = "3"
LFC_THR, PADJ_THR = 1.0, 0.05
de = markers_filtered[markers_filtered["group"] == CLUSTER_OF_INTEREST].copy()
de["neglog10_padj"] = -np.log10(de["pvals_adj"].clip(lower=1e-300))
de = de[de.logfoldchanges.abs() <= 15]   # drop zero-expression LFC artefacts (|LFC|≈30)
up   = (de.logfoldchanges >  LFC_THR) & (de.pvals_adj < PADJ_THR)
down = (de.logfoldchanges < -LFC_THR) & (de.pvals_adj < PADJ_THR)
# scatter: grey = n.s., red = up, blue = down; dashed lines at |LFC| = 1 and padj = 0.05;
# top 10 up genes labelled
plt.xlabel("log2 fold change (cluster vs rest)")
plt.ylabel("-log10(adjusted p-value)")
plt.title("Volcano plot: cluster 3 vs rest")
```

Plot: `results/q7_volcano_cluster3.png` (259 genes up, 857 down at |LFC| > 1, padj < 0.05).

I left out 1,410 genes with |LFC| > 15. Their expression is zero in one group, which makes scanpy's log fold change
blow up to about ±30. These values are artefacts, not real effects.

**Top up-regulated genes in cluster 3:** IGHD, CD37, IGHM, CD79A, MS4A1, CD79B, CD74, LINC00926, TCL1A,
HLA-DQA1/DQB1/DRB1, BANK1, CD22.

**Cell type guess: naive B cells.**
- CD79A/CD79B (B-cell receptor signalling chains), MS4A1 (CD20) and CD22 are pan-B-cell markers.
- High IGHD and IGHM together with TCL1A mean the cells express unswitched IgD⁺IgM⁺ antibodies, the signature of
  *naive* B cells.
- The MHC-II genes reflect the B cell's role as an antigen-presenting cell.

## Question 8 — Cell type annotation

I combined three lines of evidence:
- the top DEGs per cluster (printed in `results/lab4_output.txt`)
- the canonical-marker dot plot (`results/q8_dotplot.png`). I extended the marker list with FCGR3A/MS4A7 for
  FCGR3A⁺ monocytes and LILRA4/IL3RA for pDCs.
- the marker feature plots on the UMAP (`results/q8_marker_umaps.png`)

```python
cluster_labels = {
    "0": "CD14+ Monocytes",              # S100A8/9, LYZ, VCAN, FCN1, CD14
    "1": "Naive T cells (CD4/CD8 mixed)",# CD3D, CCR7, TCF7, LDHB; CD8A/B partly
    "2": "CD4+ Memory T cells",          # CD3D/E, IL7R, IL32, LTB; no CD8
    "3": "Naive B cells",                # IGHD, IGHM, TCL1A, MS4A1, CD79A
    "4": "GZMK+ CD8+ / MAIT T cells",    # CD3D, GZMK, KLRB1, KLRG1, IL7R, CCL5
    "5": "NK cells",                     # GNLY, NKG7, KLRD1, KLRF1, PRF1; no CD3D
    "6": "CD8+ Cytotoxic T cells",       # CD3D, CD8A/B, CCL5, GZMM, GZMA
    "7": "Memory B cells",               # MS4A1, CD79A, BANK1, IGKC; low IGHD/TCL1A
    "8": "FCGR3A+ Monocytes",            # FCGR3A, MS4A7, LST1, AIF1, HES4
    "9": "Conventional DCs",             # FCER1A, CST3, HLA-DR/DP high
    "10": "Platelets",                   # PPBP, PF4, TUBB1, CAVIN2
    "11": "Plasmacytoid DCs",            # LILRA4, IL3RA, TCF4, IRF7, SERPINF1
    "12": "Plasma cells",                # JCHAIN, MZB1, TNFRSF17, IGKC
}
filtered_adata.obs["cell_type"] = (
    filtered_adata.obs["leiden_clusters"].map(cluster_labels).astype("category")
)
sc.pl.umap(filtered_adata, color="cell_type", legend_loc="on data", legend_fontsize=9,
           legend_fontoutline=2, frameon=False, title="Cell Type Annotations")
```

Annotated UMAP: `results/q8_celltype_umap.png`

| Cell type | # cells |
|---|---|
| CD14+ Monocytes | 320 |
| Naive T cells (CD4/CD8 mixed) | 179 |
| CD4+ Memory T cells | 177 |
| Naive B cells | 133 |
| GZMK+ CD8+ / MAIT T cells | 79 |
| NK cells | 56 |
| CD8+ Cytotoxic T cells | 52 |
| Memory B cells | 49 |
| FCGR3A+ Monocytes | 28 |
| Conventional DCs | 28 |
| Platelets | 13 |
| Plasmacytoid DCs | 8 |
| Plasma cells | 8 |

**Challenging clusters:**

- **Cluster 1 (naive T):** its top DEGs are mostly housekeeping or translation genes (TPT1, EEF1A1, EEF1B2, NPM1)
  together with naive markers (CCR7, TCF7, LDHB). These identify the cells as naive and quiescent but don't
  separate CD4 from CD8. CD8A/CD8B are expressed at an intermediate level, so naive CD4 and naive CD8 T cells are
  probably mixed here. I gave it a broad "mixed" label.
- **Cluster 4:** it expresses CD3D, GZMK, KLRB1 (CD161), KLRG1 and IL7R. This pattern fits GZMK⁺ CD8 effector-memory
  T cells, but also MAIT cells. Telling them apart would need TRAV1-2 / SLC4A10 or TCR data, so I used a combined label.
- **Clusters 5 vs 6 (NK vs CD8 cytotoxic T):** both express NKG7, GZMA, PRF1 and CCL5. They can only be separated by
  CD3D/CD8A (present in 6, absent in 5) and by GNLY/KLRF1 (much higher in 5).
- **Clusters 9–12 (cDC, pDC, plasma cells, platelets):** these are very small (8–28 cells), so their DE statistics are
  weaker. Even so, each has highly specific markers (FCER1A, LILRA4, JCHAIN/MZB1, PPBP).
  The cDC cluster also expresses some monocyte genes (LYZ, S100A9), probably because cDCs and monocytes share
  myeloid genes, or because some doublets are present.
