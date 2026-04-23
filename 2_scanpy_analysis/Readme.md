# 🧠 Stage 2: Scanpy Single-Cell Analysis Pipeline

> **Platform:** Python / Jupyter Notebook  
> **Folder:** `2_scanpy_analysis/`  
> **Input:** Gene-cell count matrix (`.mtx` or `.h5ad`) from Stage 1  
> **Output:** Clustered and annotated AnnData object (`.h5ad`), UMAP visualizations, marker gene tables  
> **Notebook:** `scrna (1).ipynb`

---

## 📖 Table of Contents

1. [Overview](#-overview)
2. [The Scanpy Ecosystem](#-the-scanpy-ecosystem)
3. [Full Workflow at a Glance](#-full-workflow-at-a-glance)
4. [Stage-by-Stage Breakdown](#-stage-by-stage-breakdown)
   - [Stage 1 — Data Loading](#stage-1--data-loading)
   - [Stage 2 — Quality Control](#stage-2--quality-control)
   - [Stage 3 — Normalization & Transformation](#stage-3--normalization--transformation)
   - [Stage 4 — Feature Selection](#stage-4--feature-selection)
   - [Stage 5 — Dimensionality Reduction (PCA)](#stage-5--dimensionality-reduction-pca)
   - [Stage 6 — Neighborhood Graph](#stage-6--neighborhood-graph)
   - [Stage 7 — UMAP Visualization](#stage-7--umap-visualization)
   - [Stage 8 — Leiden Clustering](#stage-8--leiden-clustering)
   - [Stage 9 — Marker Gene Detection](#stage-9--marker-gene-detection)
5. [Visual Outputs Gallery](#-visual-outputs-gallery)
6. [Biological Interpretation Guide](#-biological-interpretation-guide)
7. [Parameter Reference](#-parameter-reference)
8. [Common Pitfalls & How to Avoid Them](#-common-pitfalls--how-to-avoid-them)

---

## 📌 Overview

Stage 2 is the **computational heart** of the scRNA-seq pipeline. It takes the raw count matrix produced during preprocessing and transforms it into biologically interpretable clusters of cells — each cluster representing a distinct cell type or state.

The analysis follows a canonical workflow used across the single-cell field, as standardized by the **scverse** community and validated across thousands of published studies.

**What this stage achieves:**

- Removes low-quality cells that would distort biological signals
- Normalizes expression values to make cells directly comparable
- Identifies the most informative genes and reduces noise
- Compresses 20,000+ gene dimensions into meaningful low-dimensional representations
- Groups transcriptionally similar cells into clusters
- Identifies the marker genes that define each cluster's identity

---

##  The Scanpy Ecosystem

**Scanpy** (Single-Cell Analysis in Python) is the primary Python framework for scRNA-seq analysis and is part of the broader **scverse** ecosystem — a modular, interoperable suite of tools all built around the AnnData data format.

| Package | Role in This Pipeline |
|---------|----------------------|
| `scanpy` | Core QC, normalization, PCA, clustering, visualization |
| `anndata` | Structured container for all single-cell data and metadata |
| `leidenalg` | Graph-based clustering algorithm |
| `umap-learn` | UMAP dimensionality reduction |
| `matplotlib` / `seaborn` | Plot rendering and customization |

Scanpy is the Python equivalent of **Seurat** (R) and is widely used in high-impact publications. Its tight integration with NumPy, SciPy, and pandas makes it fast, scalable, and suitable for datasets ranging from hundreds to millions of cells.

---

## 🔁 Full Workflow at a Glance

```
Count Matrix (.mtx / .h5ad)
          │
          ▼
   ┌─────────────┐
   │  Data Load  │  ← Read 10X output into AnnData object
   └──────┬──────┘
          ▼
   ┌─────────────┐
   │  Cell QC    │  ← Filter dead cells, doublets, empty droplets
   └──────┬──────┘
          ▼
   ┌──────────────────┐
   │  Normalization   │  ← Library-size normalize → log1p transform
   └──────┬───────────┘
          ▼
   ┌──────────────────┐
   │  HVG Selection   │  ← Keep top 2,000 highly variable genes
   └──────┬───────────┘
          ▼
   ┌─────────────┐
   │  Scaling    │  ← Mean-center each gene (prep for PCA)
   └──────┬──────┘
          ▼
   ┌─────────────┐
   │    PCA      │  ← Compress to 50 principal components
   └──────┬──────┘
          ▼
   ┌────────────────────┐
   │  Neighbor Graph    │  ← Build cell-cell kNN graph in PCA space
   └──────┬─────────────┘
          ▼
   ┌─────────────┐
   │    UMAP     │  ← 2D visualization of clusters
   └──────┬──────┘
          ▼
   ┌──────────────────┐
   │  Leiden Cluster  │  ← Detect communities in the kNN graph
   └──────┬───────────┘
          ▼
   ┌──────────────────┐
   │  Marker Genes    │  ← Rank genes per cluster (Wilcoxon test)
   └──────┬───────────┘
          ▼
   Annotated AnnData (.h5ad) + Visualizations
```

---

##  Stage-by-Stage Breakdown

### Stage 1 — Data Loading

**What happens:** The raw count matrix from STARsolo (output of Galaxy preprocessing) is loaded into an **AnnData object** — the foundational data structure for all downstream analysis.

**Key concepts:**
- The matrix is cells × genes (rows = cells, columns = genes)
- Gene names are made unique to handle any duplicates in the reference genome annotation
- At this point `adata.X` contains raw integer counts

**Why it matters:** Proper loading ensures the matrix dimensions and metadata are correctly structured before any transformations are applied. Errors here cascade through the entire analysis.

---

### Stage 2 — Quality Control

**What happens:** Cells that are technically problematic — not biologically meaningful — are identified and removed before analysis begins.

**Three main categories of low-quality cells:**

| Problem | Symptom | Filter Applied |
|---------|---------|---------------|
| Empty droplets | Very few genes detected | `min_genes = 200` |
| Dead / damaged cells | High mitochondrial gene fraction | `pct_counts_mt < 20%` |
| Doublets (2 cells captured together) | Abnormally high gene count | `max_genes = 5000` |

**Mitochondrial genes** are a key QC indicator — when a cell is dying or has been lysed, cytoplasmic RNA leaks out, but RNA from mitochondria (enclosed in a separate membrane) is retained, artificially inflating the mitochondrial fraction.

####  QC Scatter Plot


<p align="center">
  <img src="https://github.com/user-attachments/assets/715b8552-c1f3-46dd-8081-32fd5c5b82c8" width="500"/>
</p>


**How to read this plot:**
- Each point is one cell, plotted by total UMI counts (x-axis) vs. number of genes detected (y-axis)
- Color represents the percentage of mitochondrial counts
- Cells in the **top-right** are high quality: many genes, many counts, low mito %
- Cells in the **bottom-left or top-left with high mito %** are flagged for removal
- After filtering, the remaining population should cluster together coherently

---

### Stage 3 — Normalization & Transformation

**What happens:** Raw counts are normalized to remove technical variation introduced by differences in sequencing depth between cells, followed by a log transformation to stabilize variance.

**Two-step process:**

**Step 1 — Library-size normalization:** Each cell's counts are scaled so that all cells have the same total count (10,000 UMIs). This removes the effect of one cell having been sequenced 2× deeper than another — a technical artifact, not biology.

**Step 2 — Log1p transformation:** `log(count + 1)` is applied to every value. This compresses the highly skewed distribution of count data (most genes: 0–5 counts; some housekeeping genes: thousands) so that no single gene dominates downstream analysis. The `+1` prevents `log(0)` which is undefined.

**Why this order matters:** Normalize first (on raw counts), then log-transform. Reversing this order produces incorrect results.

---

### Stage 4 — Feature Selection

**What happens:** Of the ~20,000–33,000 genes measured, only a subset carry meaningful biological signal. **Highly Variable Genes (HVGs)** — those that vary substantially across cells — are selected for all downstream computation.

**Why not use all genes?**
- Most genes are expressed at a relatively constant level across all cell types (housekeeping genes)
- Including them adds noise without adding information
- Selecting 2,000 HVGs dramatically speeds up PCA and neighbor computation
- It focuses the analysis precisely on the genes that differentiate cell populations

Genes are ranked by their variability relative to what would be expected by chance (dispersion vs. mean relationship), and the top 2,000 are retained.

---

### Stage 5 — Dimensionality Reduction (PCA)

**What happens:** Even after HVG selection, 2,000 genes is a very high-dimensional space. **Principal Component Analysis (PCA)** finds the mathematical directions (components) that capture the most variance in the data, compressing 2,000 dimensions into 50 principal components.

**Think of it like:** Finding the "best angle" to photograph a 3D object to capture the most information in 2D. PCA finds the gene combinations that explain the most transcriptional differences between cells.

#### 📊 PCA Sample Plot

<img width="257" height="205" alt="pca_samples" src="https://github.com/user-attachments/assets/8c655bf7-8c84-4be4-8f0f-f36f48304337" />width="500"/>
</p>


**How to read this plot:**
- Each dot is one cell, positioned by its PC1 and PC2 coordinates
- Cells that separate along PC1 have the largest source of transcriptional difference
- Clear separation of groups in PCA space suggests distinct cell populations
- This is an early confirmation that the dataset contains meaningful biological structure

####  PCA Variance Ratio (Elbow Plot)

<p align="center">
  <img src="https://github.com/user-attachments/assets/f138bf90-77ae-4ae3-8c10-2f2cf350e5f0" width="500"/>
</p>


**How to read this plot:**
- X-axis: Principal component number (1, 2, 3…)
- Y-axis: Proportion of total variance explained by that component
- The curve drops steeply, then flattens — the **"elbow"** marks where informative signal ends and noise begins
- PCs beyond the elbow are discarded; this analysis uses ~15–20 PCs

---

### Stage 6 — Neighborhood Graph

**What happens:** A **k-nearest neighbor (kNN) graph** is constructed in PCA space. In this graph, each cell is connected to its 15 most similar cells. This graph structure is the foundation for both UMAP layout and Leiden clustering.

**Why PCA space, not gene space?** PCA space is more robust to noise — the components that don't capture biological variation (the ones beyond the elbow) have already been discarded.

**Key parameters:**
- `n_neighbors = 15` — more neighbors creates a "smoother" graph, fewer preserves finer local structure
- `n_pcs = 20` — use the number of informative PCs identified from the elbow plot

The resulting graph is stored internally in the AnnData object and used directly by UMAP and Leiden in the next steps.

---

### Stage 7 — UMAP Visualization

**What happens:** **Uniform Manifold Approximation and Projection (UMAP)** takes the high-dimensional kNN graph and produces a **2D layout** where cells that are transcriptionally similar are placed close together and cells that are different are placed far apart.

> ⚠️ **Critical caveat:** UMAP is a **visualization tool only**. The 2D coordinates are not meaningful for statistical analysis — only the relative proximity of cells matters. Never perform statistical tests on UMAP coordinates.

<p align="center">
  <img src="https://github.com/user-attachments/assets/b86bc0a9-2164-42f2-81b3-205eecd26792" width="500"/>
</p>

####  UMAP with Leiden Clusters

<p align="center">
  <img src="https://github.com/user-attachments/assets/985ec65a-cd30-4a82-9b83-0b9bfd826b36" width="500"/>
</p>
**How to read UMAP plots:**
- Each point = one cell
- Distinct "islands" or "blobs" = groups of transcriptionally similar cells, often corresponding to distinct cell types
- Color in `umap_4.png` corresponds to the Leiden cluster each cell was assigned to
- Cells from the same cluster should form contiguous regions in UMAP space
- Scattered or mixed patterns may indicate technical noise or the need for parameter tuning

---

### Stage 8 — Leiden Clustering

**What happens:** **Leiden clustering** applies a graph community detection algorithm to the kNN graph, identifying groups of cells that are more densely connected to each other than to the rest of the graph. Each community = one cluster.

**Why Leiden (not Louvain)?**  
Leiden is the modern successor to the Louvain algorithm — it is mathematically guaranteed to produce well-connected communities and is more stable across runs. It is the current community standard for scRNA-seq.

**The resolution parameter controls cluster granularity:**

| Resolution | Effect | Best For |
|-----------|--------|---------|
| 0.1 – 0.3 | Few, broad clusters | Quick overview of major cell types |
| 0.4 – 0.6 | Balanced granularity | Most standard analyses |
| 0.7 – 1.5 | Many fine-grained clusters | Identifying rare subtypes |

This analysis uses `resolution = 0.5`. Cluster labels are assigned as integers (0, 1, 2…) and stored in `adata.obs['leiden']`.

---

### Stage 9 — Marker Gene Detection

**What happens:** For each Leiden cluster, a **differential expression test** (Wilcoxon rank-sum) compares that cluster's gene expression against all other cells. Genes that are significantly higher in one cluster are its **marker genes** — the molecular fingerprint that defines that cell population's identity.

####  Marker Gene Dot Plot
<p align="center">
  <img src="https://github.com/user-attachments/assets/166387cb-487c-484e-bdf9-4ee316075a5f" width="700"/>
</p>

**How to read the dot plot:**
- **Rows** = Leiden clusters (0, 1, 2…)
- **Columns** = Top marker genes per cluster
- **Dot size** = Fraction of cells in the cluster expressing this gene (larger = more cells express it)
- **Dot color** = Mean expression level (darker = higher expression)
- The ideal marker gene has a **large, dark dot** in one cluster and **small, pale dots** in all others
-


**Using markers to annotate clusters:**

After identifying marker genes, cross-reference them against published cell type databases (e.g., PanglaoDB, CellMarker, the Human Cell Atlas) to assign biological identities to each cluster. For PBMC data, expected markers include:

| Cell Type | Key Marker Genes |
|-----------|----------------|
| CD4⁺ T cells | CD3D, IL7R, CD4 |
| CD8⁺ T cells | CD3D, CD8A, CD8B |
| B cells | MS4A1, CD79A |
| NK cells | GNLY, NKG7 |
| CD14⁺ Monocytes | CD14, LYZ |
| CD16⁺ Monocytes | FCGR3A, MS4A7 |
| Dendritic cells | FCER1A, CST3 |

---

##  Visual Outputs Gallery

| File | Plot Type | Biological Meaning |
|------|-----------|-------------------|
| `scatter_plot.png` | QC scatter | Shows cell quality distribution before/after filtering |
| `pca_samples.png` | PCA projection | Reveals major axes of transcriptional variation |
| `pca_variance_ratio.png` | Scree / elbow plot | Guides selection of informative PCs |
| `umap.png` | UMAP (unlabeled) | 2D cell landscape, shape of population structure |
| `umap_4.png` | UMAP (Leiden colored) | Cluster assignments overlaid on UMAP |
| `dot_plot.png` | Marker gene dot plot | Cluster-defining genes and their expression specificity |

---

##  Biological Interpretation Guide

### What makes a good cluster?

A biologically meaningful cluster should have:
1. A consistent spatial region in UMAP (not scattered)
2. Multiple high-confidence marker genes with known biological relevance
3. A plausible cell number relative to the expected biology of the tissue

### When to merge clusters?

If two clusters share nearly identical marker genes and overlap in UMAP space, they likely represent the same cell type over-split by a high resolution value. Lower the `resolution` parameter or manually merge.

### When to split clusters?

If a cluster appears large and heterogeneous in UMAP (elongated, or with visible sub-structure), and biological motivation exists, try increasing resolution or sub-clustering that cluster specifically.

---

##  Parameter Reference

| Parameter | Value Used | Exploration Range | What It Controls |
|-----------|-----------|-----------------|-----------------|
| `min_genes` | 200 | 100 – 500 | Empty droplet removal |
| `max_genes` | 5000 | 3000 – 8000 | Doublet removal |
| `pct_counts_mt` | < 20% | 10% – 25% | Dead cell removal |
| `n_top_genes` (HVG) | 2000 | 1000 – 5000 | Signal vs. noise tradeoff |
| `n_comps` (PCA) | 50 | 30 – 100 | PCA depth |
| `n_pcs` (neighbors) | 20 | 10 – 50 | Based on elbow |
| `n_neighbors` | 15 | 5 – 50 | Local structure resolution |
| `resolution` (Leiden) | 0.5 | 0.1 – 2.0 | Cluster granularity |


```python
adata.write('results/adata_post_qc.h5ad')
adata.write('results/adata_clustered.h5ad')
```

---

##  Outputs

- `scrna (1).ipynb` — Complete annotated analysis notebook
- `visuals/scatter_plot.png` — QC visualization
- `visuals/pca_samples.png` — PCA projection
- `visuals/pca_variance_ratio.png` — PCA elbow plot
- `visuals/umap.png` — UMAP embedding
- `visuals/umap_4.png` — UMAP with cluster labels
- `visuals/dot_plot.png` — Marker gene dot plot
- `adata_clustered.h5ad` *(recommended save)* — Final annotated object

---

⬅️ [Back: Preprocessing](../1_preprocessing_10X_galaxy/README.md) | [🏠 Main README](../README.md) | ➡️ [Next: AnnData](../3_AnnData/README.md)
