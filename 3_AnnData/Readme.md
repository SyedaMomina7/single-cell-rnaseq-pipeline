# 🧩 Stage 3: AnnData — Structured Single-Cell Data Management

> **Platform:** Python / Jupyter Notebook  
> **Folder:** `3_AnnData/`  
> **Notebooks:** `Getting_started_with_anndata.ipynb` · `Getting_started_with_the_anndata_package.ipynb`  
> **Input:** Any `.h5ad` file or count matrix (typically from Stage 1 or Stage 2)  
> **Output:** Annotated, structured `.h5ad` datasets ready for analysis or sharing

---

## 📖 Table of Contents

1. [Overview](#-overview)
2. [What Is AnnData?](#-what-is-anndata)
3. [The AnnData Data Structure](#-the-anndata-data-structure)
4. [Core Components In Depth](#-core-components-in-depth)
   - [`.X` — The Expression Matrix](#x--the-expression-matrix)
   - [`.obs` — Cell Metadata](#obs--cell-metadata)
   - [`.var` — Gene Metadata](#var--gene-metadata)
   - [`.obsm` — Cell Embeddings](#obsm--cell-embeddings)
   - [`.varm` — Gene Loadings](#varm--gene-loadings)
   - [`.obsp` — Cell-Cell Graphs](#obsp--cell-cell-graphs)
   - [`.uns` — Unstructured Metadata](#uns--unstructured-metadata)
   - [`.layers` — Data Layer Storage](#layers--data-layer-storage)
5. [Visual Explorations](#-visual-explorations)
6. [Key Operations](#-key-operations)
   - [Subsetting](#subsetting)
   - [Concatenation](#concatenation)
   - [Working with Layers](#working-with-layers)
   - [Views vs Copies](#views-vs-copies)
7. [Reading and Writing Files](#-reading-and-writing-files)
8. [AnnData in the scverse Ecosystem](#-anndata-in-the-scverse-ecosystem)
9. [Best Practices](#-best-practices)
10. [Outputs](#-outputs)

---

##  Overview

Stage 3 focuses on **AnnData** — the foundational data structure that underpins the entire scverse ecosystem and is used throughout this pipeline. While Scanpy (Stage 2) uses AnnData internally, this stage explores it explicitly and in depth: what it stores, how to manipulate it, how to compare raw vs. processed data, and how to handle embeddings and metadata.

Understanding AnnData is essential for:
- Reproducing analyses reliably
- Sharing datasets with collaborators in a standardized format
- Integrating data from multiple samples or studies
- Building custom analysis steps on top of Scanpy's output

---

##  What Is AnnData?

**AnnData** (Annotated Data) is a Python data container designed specifically for matrix-shaped data with row and column annotations — making it a perfect fit for single-cell data, where:

- **Rows** = cells (observations)
- **Columns** = genes (variables/features)
- **Annotations** = any metadata, embeddings, or graphs associated with cells or genes

AnnData was developed as part of the **scverse** project and is the standard data format used by Scanpy, scVI-tools, squidpy, cellrank, and dozens of other tools. Datasets are saved as **`.h5ad` files** — a compact, efficient binary format based on HDF5.

**Why AnnData over a plain matrix?**

A plain matrix (NumPy or pandas DataFrame) stores only the numbers. AnnData stores the matrix *and* all the context that makes it interpretable: which cells came from which patient, which genes are highly variable, what the PCA coordinates are, what the kNN graph looks like, what cluster each cell belongs to — all in one coherent, self-describing object.

---

##  The AnnData Data Structure

```
AnnData object
│
├── .X          ←  Main expression matrix  (n_cells × n_genes)
│
├── .obs         ←  Cell annotations  (n_cells × n_obs_features)
│                    e.g., batch, cluster, cell type, QC metrics
│
├── .var         ←  Gene annotations  (n_genes × n_var_features)
│                    e.g., highly_variable, mean, dispersion
│
├── .obsm        ←  Cell embeddings  (dict of 2D arrays)
│                    e.g., 'X_pca', 'X_umap'
│
├── .varm        ←  Gene-level matrices  (dict of 2D arrays)
│                    e.g., 'PCs' (gene loadings)
│
├── .obsp        ←  Cell-cell sparse graphs  (dict of sparse matrices)
│                    e.g., 'distances', 'connectivities'
│
├── .uns         ←  Unstructured metadata  (dict)
│                    e.g., color palettes, analysis parameters, DE results
│
└── .layers      ←  Additional data matrices  (dict of matrices, same shape as .X)
                     e.g., 'raw_counts', 'log1p_norm', 'scaled'
```

At any point in an analysis, you can inspect the full object:

```python
print(adata)
# AnnData object with n_obs × n_vars = 5432 × 2000
#     obs: 'n_genes', 'n_counts', 'pct_counts_mt', 'leiden', 'cell_type'
#     var: 'highly_variable', 'means', 'dispersions_norm'
#     uns: 'leiden_colors', 'rank_genes_groups'
#     obsm: 'X_pca', 'X_umap'
#     obsp: 'distances', 'connectivities'
#     layers: 'raw_counts', 'log1p_norm'
```

---

##  Core Components In Depth

### `.X` — The Expression Matrix

`.X` is the **primary data matrix** — the heart of the AnnData object. Its content changes as the analysis progresses:

| Analysis Stage | What `.X` Contains |
|---------------|-------------------|
| After loading | Raw integer counts |
| After normalization | Library-size normalized values |
| After log1p | Log-normalized values |
| After scaling | Mean-centered, variance-stabilized values |

**Best practice:** Because `.X` gets overwritten at each step, save raw counts and normalized values in `.layers` before proceeding:

```python
adata.layers['raw_counts'] = adata.X.copy()
# ... after normalization ...
adata.layers['log1p_norm'] = adata.X.copy()
```

---

### `.obs` — Cell Metadata

`.obs` is a **pandas DataFrame** with one row per cell. It accumulates metadata throughout the analysis:

| Column | Added By | Meaning |
|--------|---------|---------|
| `n_genes_by_counts` | QC step | Number of genes detected |
| `total_counts` | QC step | Total UMI count (library size) |
| `pct_counts_mt` | QC step | % mitochondrial reads |
| `leiden` | Clustering | Cluster assignment (0, 1, 2…) |
| `cell_type` | Annotation | Biological cell type label |

Access and filter:
```python
# Inspect metadata
adata.obs.head()

# Filter to T cells only
t_cells = adata[adata.obs['cell_type'] == 'CD4+ T cells']
```

---

### `.var` — Gene Metadata

`.var` is a **pandas DataFrame** with one row per gene:

| Column | Added By | Meaning |
|--------|---------|---------|
| `highly_variable` | HVG selection | True if gene is in top HVGs |
| `means` | QC metrics | Mean expression across all cells |
| `dispersions_norm` | HVG selection | Normalized dispersion score |
| `mt` | QC step | True if gene is mitochondrial |

---

### `.obsm` — Cell Embeddings

`.obsm` stores **2D arrays of cell coordinates** — one row per cell, multiple columns per embedding:

| Key | Shape | Stored By |
|-----|-------|-----------|
| `X_pca` | (n_cells, n_pcs) | `sc.tl.pca()` |
| `X_umap` | (n_cells, 2) | `sc.tl.umap()` |
| `X_tsne` | (n_cells, 2) | `sc.tl.tsne()` (if used) |

Access UMAP coordinates directly:
```python
umap_coords = adata.obsm['X_umap']  # numpy array, shape (n_cells, 2)
```

---

### `.varm` — Gene Loadings

`.varm` stores gene-level matrices, most commonly the **PCA loadings** (how much each gene contributes to each PC):

```python
pca_loadings = adata.varm['PCs']  # shape (n_genes, n_pcs)
```

This is useful for identifying which genes drive separation along specific principal components.

---

### `.obsp` — Cell-Cell Graphs

`.obsp` stores **sparse cell × cell matrices** representing graph connectivity. After running `sc.pp.neighbors()`:

| Key | Meaning |
|-----|---------|
| `distances` | Actual distances to kNN neighbors |
| `connectivities` | UMAP/Leiden-ready affinity weights |

These are typically not accessed directly but are used internally by UMAP and Leiden.

---

### `.uns` — Unstructured Metadata

`.uns` is a flexible dictionary for anything that doesn't fit neatly into a matrix:

```python
adata.uns['leiden_colors']        # Color palette for clusters
adata.uns['rank_genes_groups']    # Marker gene results (dict of arrays)
adata.uns['pca']                  # PCA variance ratios
adata.uns['neighbors']            # Neighbor computation parameters
```

When performing custom analyses, use `.uns` to store results or configuration parameters alongside the data.

---

### `.layers` — Data Layer Storage

`.layers` allows multiple versions of the expression matrix to coexist, all with the same shape as `.X`:

```python
adata.layers['raw_counts']   # Original integer counts
adata.layers['log1p_norm']   # Normalized and log-transformed
adata.layers['scaled']       # Mean-centered (if saved)
```

Layers are essential for **visualization** — when making UMAP feature plots, you want to show log-normalized expression (interpretable) rather than scaled values (which can be negative and are not meaningful for visualization).

---

## 📊 Visual Explorations

### 🔹 Raw vs. Normalized Counts Comparison

<img width="611" height="151" alt="cpm_vs_raw_counts" src="https://github.com/user-attachments/assets/bbfb3ddf-30c7-4ebf-b080-a319c16c2ae6" />

**What this plot shows:**

This side-by-side comparison visualizes the effect of library-size normalization. Raw counts show high variability between cells purely due to differences in sequencing depth — a technical artifact. After normalization (scaling to counts per million, or a fixed total of 10,000 counts), cells become directly comparable. Key observations:

- **Before normalization:** The same gene can appear to have 5 counts in one cell and 50 in another, simply because one cell was sequenced 10× deeper
- **After normalization:** Relative expression differences reflect true biology, not technical variation
- **After log1p:** The distribution compresses, making moderate and highly expressed genes visible on the same scale

This visualization provides justification for the normalization step and is a critical quality check in any published analysis.

---

### 🔹 Distance Matrix

<img width="379" height="283" alt="distance_matrix_raw" src="https://github.com/user-attachments/assets/427a422a-1e58-429c-9127-f3bda98353bf" />


**What this plot shows:**

A **pairwise distance matrix** (heatmap) showing how similar every cell is to every other cell based on their gene expression profiles. Each row and column represents one cell; color intensity represents expression distance (darker = more similar, lighter = more different).

Key patterns to look for:
- **Diagonal blocks** of high similarity → distinct cell populations that will become clusters
- **Off-diagonal blocks** of moderate similarity → related cell types with shared gene programs
- A completely uniform matrix → no structure, possibly indicative of failed sequencing or extreme batch effects
- This matrix is computed on **raw** data (before dimensionality reduction), providing an unbiased view of transcriptional similarity

This is the precursor to the kNN graph used in Leiden clustering — in practice, computing all pairwise distances for large datasets is computationally prohibitive, which is why PCA-space approximations and kNN graphs are used instead.

---

### 🔹 Embeddings Visualization

<img width="730" height="386" alt="embeddings_plot" src="https://github.com/user-attachments/assets/01b92036-2eb3-48e6-a6fd-f9698f2aca07" />


**What this plot shows:**

A comparison of different **dimensionality reduction embeddings** applied to the same dataset. The `.obsm` slot stores multiple embeddings (PCA, UMAP, t-SNE), and this plot displays them side by side to illustrate how each captures different aspects of the data structure.

Key differences between embeddings:

| Embedding | Strength | Weakness |
|-----------|---------|---------|
| **PCA** | Globally accurate, linear | Cannot capture non-linear structure |
| **UMAP** | Captures local clusters well, fast | Distances not globally meaningful |
| **t-SNE** | Good local cluster separation | Slow, clusters may be artificially far apart |

For scRNA-seq, **UMAP is the standard visualization** — it is faster than t-SNE, better preserves global structure, and is more reproducible. This plot demonstrates why: compare how well each embedding separates the known cell populations.

---

### 🔹 Sorted Heatmap of Gene Expression

<img width="324" height="248" alt="sorted heatmap" src="https://github.com/user-attachments/assets/7edbd7bd-3158-4731-a23f-054809a31ba9" />

**What this plot shows:**

A **hierarchically sorted heatmap** of gene expression across all cells, with genes as columns and cells as rows (or the transpose). Cells are sorted by cluster assignment; genes are sorted by which cluster they are most highly expressed in.

How to interpret:
- **Horizontal bands** of high expression (yellow/red) = a cluster with a characteristic set of marker genes
- **Vertical stripes** = genes expressed broadly across all cells (housekeeping genes — these would typically be excluded from this plot)
- **Sparse regions** (dark/purple) = low or absent expression — the typical state for most genes in most cells (single-cell data is ~90% zeros, a property called **sparsity**)
- The sorted order reveals the **block-diagonal structure** of the expression matrix that corresponds to distinct cell types

This visualization is one of the most powerful ways to visually confirm that clustering has identified biologically meaningful populations.

---

## ⚙️ Key Operations

### Subsetting

AnnData supports NumPy-style slicing. The result is a **view** (see below):

```python
# Select specific cells (by index or boolean mask)
t_cells = adata[adata.obs['cell_type'] == 'T cells']

# Select specific genes
hvg_data = adata[:, adata.var['highly_variable']]

# Select both simultaneously
subset = adata[adata.obs['leiden'] == '0', ['CD3D', 'CD4', 'IL7R']]
```

---

### Concatenation

Merging multiple samples into one AnnData object:

```python
import anndata as ad

# Concatenate along the cell axis (stack cells)
combined = ad.concat(
    [adata_sample1, adata_sample2],
    label='sample',          # Adds a 'sample' column to .obs
    keys=['patient_A', 'patient_B']
)
```

After concatenation, batch correction is typically needed before joint analysis.

---

### Working with Layers

```python
# Store raw counts before normalization
adata.layers['raw'] = adata.X.copy()

# Store normalized data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log_norm'] = adata.X.copy()

# For visualization, use log_norm (not scaled .X)
sc.pl.umap(adata, color='CD3D', layer='log_norm')
```

---

### Views vs Copies

This is one of the most important and commonly misunderstood aspects of AnnData:

| Operation | Result | Modifies Original? | Memory |
|-----------|--------|--------------------|--------|
| `adata[mask]` | **View** | Yes (changes propagate) | Minimal |
| `adata[mask].copy()` | **Copy** | No (independent) | Full copy |

**When to use a copy:** Any time you want to modify a subset without affecting the original object:

```python
# This modifies the original (dangerous)
subset = adata[adata.obs['leiden'] == '0']

# This creates an independent object (safe)
subset = adata[adata.obs['leiden'] == '0'].copy()
```

---

## 💾 Reading and Writing Files

AnnData uses the **`.h5ad` format** — an HDF5-based binary format that efficiently stores sparse matrices, metadata DataFrames, and nested dictionaries in a single file.

```python
# Save (write)
adata.write('data/adata_final.h5ad', compression='gzip')

# Load (read)
import anndata as ad
adata = ad.read_h5ad('data/adata_final.h5ad')

# Read from 10X Genomics output
adata = sc.read_10x_mtx('path/to/filtered_feature_bc_matrix/')

# Read from CSV (slower, for compatibility)
adata = sc.read_csv('counts.csv').T  # Transpose if genes × cells
```

**File size considerations:** A dataset with 10,000 cells × 33,000 genes might be 500 MB in memory as a dense matrix, but only 20–50 MB on disk in sparse `.h5ad` format. Always use `compression='gzip'` when saving.

---

## 🌐 AnnData in the scverse Ecosystem

AnnData serves as the **shared data format** that allows all scverse tools to interoperate:

```
         ┌───────────────┐
         │    AnnData    │  ← The universal format
         └───────┬───────┘
                 │
    ┌────────────┼────────────┐
    ▼            ▼            ▼
 Scanpy      scVI-tools    squidpy
(Analysis)  (Deep learning) (Spatial)
    │            │            │
    └────────────┴────────────┘
                 │
           Single .h5ad file
         (portable, shareable)
```

This means a dataset analyzed in Scanpy can be directly loaded into scVI-tools for deep learning-based integration, or into squidpy for spatial analysis — no conversion needed.

---


##  Outputs

| File | Description |
|------|-------------|
| `Getting_started_with_anndata.ipynb` | Hands-on tutorial covering core AnnData operations |
| `Getting_started_with_the_anndata_package.ipynb` | Extended exploration including subsetting, concatenation, and layers |
| `visuals/cpm_vs_raw_counts.png` | Normalization effect visualization |
| `visuals/distance_matrix_raw.png` | Pairwise cell distance heatmap |
| `visuals/embeddings_plot.png` | Comparison of dimensionality reduction methods |
| `visuals/sorted_heatmap.png` | Gene expression heatmap sorted by cluster |
| `*.h5ad` | Annotated AnnData output files |

---

## 📚 Further Reading

- [AnnData Documentation](https://anndata.readthedocs.io/en/latest/)
- [scverse Ecosystem Overview](https://scverse.org/)



