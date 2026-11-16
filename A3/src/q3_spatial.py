"""COMP3353 A3 Q3: spatial transcriptomics analysis of DLPFC slide 151673.

(a) QC + preprocessing  (b) kNN graph, UMAP, Leiden, ARI
(c) spatially smoothed clustering  (d) spatial marker gene discovery

Usage:  python q3_spatial.py        (run from anywhere; paths are relative to A3/)
"""
import os
import warnings

warnings.filterwarnings("ignore")

import cv2
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "Q3")
RES = os.path.join(HERE, "..", "results")
SEED = 0

# ---- chosen parameters -----------------------------------------------------
MIN_UMI = 500          # spot filter (Layer 1 is naturally low-UMI, median ~1.7k)
MIN_CELLS = 10         # gene filter
N_HVG = 3000
N_PCS = 30
N_NEIGHBORS = 15
RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
RES_B = 0.8            # unsmoothed Leiden resolution (-> 6 clusters)
SELF_W = 0.5           # smoothing: 50% own PCs + 50% mean of hex neighbours
HEX_RADIUS = 31        # px; nearest-neighbour spacing is ~20.6 px, next ring ~35.7 px
RES_C = 0.3            # smoothed Leiden resolution (-> 6 clusters)
LAYER_ORDER = ["Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM"]

log_lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log_lines.append(s)


# ----------------------------------------------------------------- plotting
IMG = None


def tissue_image():
    global IMG
    if IMG is None:
        IMG = cv2.cvtColor(cv2.imread(os.path.join(DATA, "151673_tissue_hires_image.png")),
                           cv2.COLOR_BGR2RGB)
    return IMG


def spatial_plot(adata, ax, values=None, categories=None, title="", cmap="viridis",
                 order=None, s=11):
    """Scatter spots over the H&E image. As in the handout snippet, the
    'x' coordinate is the image row and 'y' the image column -> scatter(y, x)."""
    ax.imshow(tissue_image(), alpha=0.55)
    x, y = adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1]
    if categories is not None:
        cats = order or sorted(pd.unique(categories), key=str)
        colors = plt.get_cmap("tab10" if len(cats) <= 10 else "tab20")
        for k, cat in enumerate(cats):
            m = np.asarray(categories) == cat
            ax.scatter(y[m], x[m], s=s, color=colors(k), label=str(cat), linewidths=0)
        ax.legend(markerscale=3, fontsize=8, loc="upper right", frameon=True)
    else:
        lo, hi = np.percentile(values, [1, 99])
        sca = ax.scatter(y, x, c=values, s=s, cmap=cmap, linewidths=0, vmin=lo, vmax=hi)
        plt.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    ax.set_xlim(y.min() - 60, y.max() + 60)
    ax.set_ylim(x.max() + 60, x.min() - 60)
    ax.set_title(title)


def save_pair(adata, left, right, fname, ltitle, rtitle, rorder=None):
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    spatial_plot(adata, axes[0], categories=left, title=ltitle, order=LAYER_ORDER)
    spatial_plot(adata, axes[1], categories=right, title=rtitle, order=rorder)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, fname), dpi=110)
    plt.close(fig)


def leiden(adata, res, key):
    sc.tl.leiden(adata, resolution=res, flavor="igraph", n_iterations=2,
                 directed=False, random_state=SEED, key_added=key)


# ----------------------------------------------------------------- (a)
def part_a():
    log("=" * 70 + "\n(a) QC and preprocessing\n" + "=" * 70)
    adata = sc.read_10x_h5(os.path.join(DATA, "151673_raw_feature_bc_matrix.h5"))
    adata.var_names_make_unique()
    coord = pd.read_csv(os.path.join(DATA, "151673_spot_coordinates.csv"), index_col=0)
    ann = pd.read_csv(os.path.join(DATA, "151673_annotations.csv"), index_col=0)
    log(f"raw matrix: {adata.n_obs} barcodes x {adata.n_vars} genes")

    adata = adata[coord.index].copy()
    adata.obsm["spatial"] = coord.loc[adata.obs_names, ["x", "y"]].to_numpy()
    adata.obs["layer"] = ann.loc[adata.obs_names, "layer"]
    adata = adata[adata.obs["layer"].notna()].copy()
    adata.obs["layer"] = pd.Categorical(adata.obs["layer"], categories=LAYER_ORDER)
    log(f"spots with coordinates + annotation: {adata.n_obs}")

    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    q = np.percentile(adata.obs["total_counts"], [0, 1, 5, 50])
    log(f"UMI/spot: min {q[0]:.0f}, 1% {q[1]:.0f}, 5% {q[2]:.0f}, median {q[3]:.0f}")
    log("median UMI per layer: " + ", ".join(
        f"{k} {v:.0f}" for k, v in adata.obs.groupby("layer", observed=True)["total_counts"].median().items()))

    # QC histogram
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(np.log10(adata.obs["total_counts"]), bins=80, color="C0")
    ax.axvline(np.log10(MIN_UMI), color="r", ls="--", label=f"threshold = {MIN_UMI} UMI")
    ax.set_xlabel("log10(total UMI per spot)")
    ax.set_ylabel("# spots")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "q3a_umi_histogram.png"), dpi=130)
    plt.close(fig)

    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=MIN_UMI)
    log(f"spot filter (total UMI >= {MIN_UMI}): removed {n0 - adata.n_obs}, kept {adata.n_obs}")
    g0 = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    log(f"gene filter (expressed in >= {MIN_CELLS} spots): kept {adata.n_vars} / {g0} genes")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata                   # full log-normalised matrix X (for DE)

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
    hvg = adata.var_names[adata.var["highly_variable"]]
    top = adata.var.loc[hvg].sort_values("dispersions_norm", ascending=False).index[:20]
    log(f"HVGs kept: {len(hvg)} ; top 20 by normalised dispersion: {', '.join(top)}")
    adata = adata[:, hvg].copy()

    sc.pp.pca(adata, n_comps=N_PCS, random_state=SEED)
    vr = adata.uns["pca"]["variance_ratio"]
    log(f"PCA: {N_PCS} PCs, cumulative variance explained = {vr.sum():.3f} "
        f"(PC1 {vr[0]:.3f}, PC2 {vr[1]:.3f})")
    log(f"final matrix X: {adata.n_obs} spots x {adata.n_vars} HVGs; "
        f"X_pca: {adata.obsm['X_pca'].shape}")
    log("PC1 vs layer (mean): " + ", ".join(
        f"{k} {v:+.2f}" for k, v in pd.Series(adata.obsm['X_pca'][:, 0], index=adata.obs_names)
        .groupby(adata.obs["layer"], observed=True).mean().items()))

    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    ax[0].plot(np.arange(1, N_PCS + 1), vr, "o-")
    ax[0].set_xlabel("PC")
    ax[0].set_ylabel("variance ratio")
    ax[0].set_title("PCA elbow plot")
    for k, lab in enumerate(LAYER_ORDER):
        m = adata.obs["layer"] == lab
        ax[1].scatter(adata.obsm["X_pca"][m, 0], adata.obsm["X_pca"][m, 1], s=3,
                      color=plt.get_cmap("tab10")(k), label=lab)
    ax[1].set_xlabel("PC1")
    ax[1].set_ylabel("PC2")
    ax[1].legend(markerscale=4, fontsize=8)
    ax[1].set_title("PC1 vs PC2 coloured by annotation")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "q3a_pca.png"), dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 8))
    spatial_plot(adata, ax, values=adata.obsm["X_pca"][:, 0], cmap="coolwarm",
                 title="PC1 value per spot")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "q3a_pc1_spatial.png"), dpi=130)
    plt.close(fig)
    return adata


# ----------------------------------------------------------------- (b)
def part_b(adata):
    log("\n" + "=" * 70 + "\n(b) kNN graph, UMAP, Leiden\n" + "=" * 70)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca", metric="euclidean",
                    random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    log(f"kNN graph: k = {N_NEIGHBORS}, Euclidean distance on {N_PCS} PCs")
    log("resolution scan (resolution, #clusters, ARI):")
    for r in RESOLUTIONS:
        leiden(adata, r, "tmp")
        log(f"   {r:.1f}  {adata.obs['tmp'].nunique():2d}  "
            f"{adjusted_rand_score(adata.obs['layer'], adata.obs['tmp']):.3f}")
    del adata.obs["tmp"]
    leiden(adata, RES_B, "leiden")
    ari = adjusted_rand_score(adata.obs["layer"], adata.obs["leiden"])
    log(f"chosen resolution {RES_B}: {adata.obs['leiden'].nunique()} clusters, ARI = {ari:.4f}")
    log("contingency (rows = cluster):\n" +
        pd.crosstab(adata.obs["leiden"], adata.obs["layer"]).to_string())

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
    for a, key, title in ((ax[0], "layer", "UMAP: ground truth"),
                          (ax[1], "leiden", f"UMAP: Leiden (res {RES_B})")):
        cats = adata.obs[key].cat.categories
        for k, c in enumerate(cats):
            m = adata.obs[key] == c
            a.scatter(*adata.obsm["X_umap"][m].T, s=3, color=plt.get_cmap("tab10")(k), label=c)
        a.legend(markerscale=4, fontsize=8)
        a.set_title(title)
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "q3b_umap.png"), dpi=130)
    plt.close(fig)

    save_pair(adata, adata.obs["layer"], adata.obs["leiden"], "q3b_spatial_clusters.png",
              "Ground truth annotation", f"Leiden clusters (res {RES_B}), ARI = {ari:.3f}",
              rorder=list(adata.obs["leiden"].cat.categories))
    return ari


# ----------------------------------------------------------------- (c)
def hex_adjacency(adata):
    """Spots within HEX_RADIUS px are immediate hex-grid neighbours (<= 6)."""
    xy = adata.obsm["spatial"]
    pairs = cKDTree(xy).query_pairs(HEX_RADIUS, output_type="ndarray")
    n = adata.n_obs
    A = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return (A + A.T).tocsr()


def part_c(adata):
    log("\n" + "=" * 70 + "\n(c) Spatially smoothed clustering\n" + "=" * 70)
    A = hex_adjacency(adata)
    deg = np.asarray(A.sum(axis=1)).ravel()
    log(f"hex neighbours per spot: " +
        ", ".join(f"{k}:{v}" for k, v in sorted(pd.Series(deg).value_counts().items())))
    W = sp.diags(1.0 / np.maximum(deg, 1)) @ A                # row-normalised
    Z = adata.obsm["X_pca"]
    nb_mean = W @ Z
    has_nb = (deg > 0)[:, None]
    adata.obsm["X_smooth"] = np.where(has_nb, SELF_W * Z + (1 - SELF_W) * nb_mean, Z)

    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_smooth", random_state=SEED,
                    key_added="smooth")
    log(f"smoothed PCs: {SELF_W:.0%} own + {1 - SELF_W:.0%} mean of hex neighbours; "
        f"kNN (k = {N_NEIGHBORS}) rebuilt on smoothed PCs")
    log("resolution scan (resolution, #clusters, ARI):")
    for r in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=r, flavor="igraph", n_iterations=2, directed=False,
                     random_state=SEED, key_added="tmp", neighbors_key="smooth")
        log(f"   {r:.1f}  {adata.obs['tmp'].nunique():2d}  "
            f"{adjusted_rand_score(adata.obs['layer'], adata.obs['tmp']):.3f}")
    del adata.obs["tmp"]
    sc.tl.leiden(adata, resolution=RES_C, flavor="igraph", n_iterations=2, directed=False,
                 random_state=SEED, key_added="leiden_smooth_raw", neighbors_key="smooth")
    ari_raw = adjusted_rand_score(adata.obs["layer"], adata.obs["leiden_smooth_raw"])
    log(f"step 1 (feature smoothing), res {RES_C}: "
        f"{adata.obs['leiden_smooth_raw'].nunique()} clusters, ARI = {ari_raw:.4f}")

    # step 2: label refinement - majority vote among a spot and its hex neighbours
    labels = adata.obs["leiden_smooth_raw"].astype(int).to_numpy()
    k = labels.max() + 1
    onehot = sp.csr_matrix((np.ones(len(labels)), (np.arange(len(labels)), labels)),
                           shape=(len(labels), k))
    votes = (A + sp.identity(A.shape[0])) @ onehot
    votes = votes.toarray()
    new = labels.copy()
    own = votes[np.arange(len(labels)), labels]
    best = votes.argmax(axis=1)
    flip = votes.max(axis=1) > own          # neighbourhood majority disagrees
    new[flip] = best[flip]
    adata.obs["leiden_smooth"] = pd.Categorical(new.astype(str),
                                                categories=sorted(set(new.astype(str)), key=int))
    ari = adjusted_rand_score(adata.obs["layer"], adata.obs["leiden_smooth"])
    log(f"step 2 (neighbour majority vote on labels): {flip.sum()} spots relabelled, "
        f"{adata.obs['leiden_smooth'].nunique()} clusters, ARI = {ari:.4f}")
    log("contingency (rows = smoothed cluster):\n" +
        pd.crosstab(adata.obs["leiden_smooth"], adata.obs["layer"]).to_string())

    save_pair(adata, adata.obs["layer"], adata.obs["leiden_smooth"],
              "q3c_spatial_smoothed_clusters.png", "Ground truth annotation",
              f"Spatially smoothed Leiden (res {RES_C}), ARI = {ari:.3f}",
              rorder=list(adata.obs["leiden_smooth"].cat.categories))
    return ari


# ----------------------------------------------------------------- (d)
def match_cluster(adata, key, layer):
    """Cluster with the highest Jaccard overlap with a ground-truth layer."""
    gt = adata.obs["layer"] == layer
    best = None
    for c in adata.obs[key].cat.categories:
        m = adata.obs[key] == c
        j = (gt & m).sum() / (gt | m).sum()
        if best is None or j > best[1]:
            best = (c, j)
    return best


def top_degs(adata, groupby, group, n=10):
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=[group], reference="rest",
                            method="wilcoxon", use_raw=True, key_added="rg")
    df = sc.get.rank_genes_groups_df(adata, group=group, key="rg")
    df = df[df["logfoldchanges"] > 0]                   # up-regulated markers
    return df.head(n).reset_index(drop=True)


def part_d(adata):
    log("\n" + "=" * 70 + "\n(d) Spatial marker genes (Wilcoxon one-vs-rest, "
        "ranked by adjusted p-value / z-score, up-regulated only)\n" + "=" * 70)
    key = "leiden_smooth"
    tables = {}
    for layer in ("WM", "Layer1"):
        tables[("gt", layer)] = top_degs(adata, "layer", layer)
        c, j = match_cluster(adata, key, layer)
        tables[("cl", layer)] = top_degs(adata, key, c)
        tables[("clid", layer)] = (c, j)

    def show(df):
        return df[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(
            index=False, formatters={"scores": "{:.2f}".format,
                                     "logfoldchanges": "{:.2f}".format,
                                     "pvals_adj": "{:.2e}".format})

    for layer in ("WM", "Layer1"):
        c, j = tables[("clid", layer)]
        g = tables[("gt", layer)]
        k = tables[("cl", layer)]
        shared = [x for x in g["names"] if x in set(k["names"])]
        log(f"\n--- {layer}: ground-truth annotation '{layer}' vs rest ---\n{show(g)}")
        log(f"\n--- smoothed cluster {c} (best Jaccard match to {layer}, J = {j:.2f}) vs rest ---\n{show(k)}")
        log(f"overlap of top-10 lists: {len(shared)}/10 -> {shared}")

    wm_gene = next(x for x in tables[("gt", "WM")]["names"] if x in set(tables[("cl", "WM")]["names"]))
    l1_gene = next(x for x in tables[("gt", "Layer1")]["names"] if x in set(tables[("cl", "Layer1")]["names"]))
    log(f"\nselected representative markers: WM -> {wm_gene}, Layer1 -> {l1_gene}")

    fig, axes = plt.subplots(1, 2, figsize=(17, 8))
    for ax, g, lab in ((axes[0], wm_gene, "WM"), (axes[1], l1_gene, "Layer 1")):
        expr = adata.raw[:, g].X
        expr = np.asarray(expr.todense()).ravel() if sp.issparse(expr) else np.asarray(expr).ravel()
        spatial_plot(adata, ax, values=expr, cmap="magma_r",
                     title=f"{g} (log-normalised expression) - {lab} marker")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "q3d_marker_genes_spatial.png"), dpi=110)
    plt.close(fig)
    return tables


if __name__ == "__main__":
    os.makedirs(RES, exist_ok=True)
    sc.settings.verbosity = 0
    ad = part_a()
    part_b(ad)
    part_c(ad)
    part_d(ad)
    with open(os.path.join(RES, "q3_output.txt"), "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
