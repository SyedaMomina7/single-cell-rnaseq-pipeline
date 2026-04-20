# End-to-End Single-Cell RNA-Seq Analysis Workflow

**Galaxy (10x preprocessing) · Scanpy · AnnData · scverse**

---

##  Overview

This repository presents a **comprehensive and reproducible single-cell RNA sequencing (scRNA-seq) workflow**, covering the full pipeline from **raw sequencing reads (FASTQ)** to **biologically meaningful cell clusters and annotations**.

The workflow integrates:

* **Galaxy (GTN)** for upstream preprocessing of raw 10x Genomics data
* **Scanpy (Python)** for downstream statistical analysis and clustering
* **AnnData (scverse ecosystem)** for structured data representation and manipulation

This hybrid design reflects **real-world scRNA-seq analysis pipelines**, combining **user-friendly preprocessing** with **flexible computational analysis**.

---

##  Workflow Summary

```
Raw FASTQ (10x Genomics)
        ↓
[Galaxy Preprocessing]
(Barcode extraction, alignment, UMI counting)
        ↓
Gene-Cell Count Matrix → AnnData (.h5ad)
        ↓
[Scanpy Analysis]
(QC → Normalization → HVG → PCA → UMAP → Clustering)
        ↓
[AnnData Exploration]
(Data structure, metadata handling, annotation)
```

---

##  Repository Structure

```
scRNA-seq-pipeline/
│
├── 1_preprocessing_10X_galaxy/
│   └── scRNA-galaxy-files.zip
│
├── 2_scanpy_analysis/
│   └── basic-scrna-tutorial.ipynb
│
├── 3_anndata_learning/
│   ├── anndata1.ipynb
│   └── anndata2.ipynb
│
├── README.md
└── requirements.txt
```

---

## Stage 1 — Upstream Preprocessing (Galaxy)

**Platform:** Galaxy Training Network (GTN)
**Objective:** Transform raw sequencing reads into a structured gene-cell count matrix

### Methodology

* Input: Paired-end FASTQ files

  * **R1** → cell barcodes + UMIs
  * **R2** → transcript sequences

* Alignment performed using **STARsolo**, a widely used alternative to Cell Ranger

* Key operations:

  * Barcode identification and filtering
  * UMI deduplication
  * Gene-level quantification

* Quality control performed using **MultiQC**, evaluating:

  * Mapping efficiency
  * Number of detected cell barcodes
  * Filtering thresholds

* Output converted into **AnnData (`.h5ad`) format**, enabling compatibility with Scanpy

### Output Artifacts

* Gene-cell count matrix
* Barcode and feature annotations
* `.h5ad` file for downstream analysis

---

##  Stage 2 — Downstream Analysis (Scanpy)

**Notebook:** `basic-scrna-tutorial.ipynb`

This stage performs the **core statistical and computational analysis** of scRNA-seq data.

### Analytical Pipeline

#### 🔹 Quality Control (QC)

* Computation of per-cell metrics:

  * Number of detected genes (`n_genes_by_counts`)
  * Total counts (`total_counts`)
  * Mitochondrial gene percentage (`pct_counts_mt`)
* Removal of:

  * Low-quality cells
  * Rarely expressed genes

#### 🔹 Normalization

* Library size normalization to equalize sequencing depth
* Logarithmic transformation for variance stabilization

#### 🔹 Feature Selection

* Identification of **Highly Variable Genes (HVGs)**
* Focus on biologically informative genes

#### 🔹 Dimensionality Reduction

* Principal Component Analysis (PCA)
* Selection of informative components

#### 🔹 Graph Construction

* k-nearest neighbor (kNN) graph
* Captures cell-to-cell similarity structure

#### 🔹 Clustering & Visualization

* UMAP for nonlinear embedding
* Leiden algorithm for community detection

#### 🔹 Marker Gene Identification

* Statistical testing (e.g., Wilcoxon rank-sum)
* Identification of cluster-specific gene signatures

### Outputs

* Processed AnnData object
* Low-dimensional embeddings (UMAP)
* Cluster assignments and marker genes

---

##  Stage 3 — AnnData Exploration & Data Handling

**Notebooks:**

* `anndata1.ipynb`
* `anndata2.ipynb`

This stage focuses on **deep understanding and manipulation of the AnnData structure**, which underpins modern single-cell analysis workflows.

### Core Concepts

#### 🔹 Data Structure

* `.X` → expression matrix
* `.obs` → cell-level metadata
* `.var` → gene-level metadata
* `.obsm` → embeddings (e.g., PCA, UMAP)
* `.layers` → alternative data representations

#### 🔹 Data Manipulation

* Subsetting cells and genes
* Adding annotations (e.g., cell types)
* Working with categorical variables for efficiency

#### 🔹 Advanced Features

* Sparse matrix optimization
* Backed mode for memory-efficient file handling
* View vs copy behavior in slicing operations

#### 🔹 File Operations

* Reading and writing `.h5ad` files
* Managing large-scale datasets

### Outputs

* Structured and annotated AnnData objects
* Demonstrations of efficient data handling practices

---

##  Skills & Concepts Demonstrated

* End-to-end scRNA-seq workflow design
* Preprocessing of 10x Genomics data
* Statistical analysis and clustering of single-cell data
* Dimensionality reduction and visualization techniques
* Efficient handling of high-dimensional biological datasets
* Practical understanding of AnnData data model

---

##  Environment Setup

```bash
pip install scanpy anndata numpy pandas matplotlib scipy
```

---

##  Reproducibility & Execution

1. Perform preprocessing using Galaxy (external platform)
2. Execute analysis notebook:

   ```
   basic-scrna-tutorial.ipynb
   ```
3. Run AnnData exploration notebooks:

   ```
   anndata1.ipynb
   anndata2.ipynb
   ```

---

##  References

* Galaxy Training Network — scRNA-seq preprocessing
* Scanpy documentation (scverse)
* AnnData documentation
* scverse tutorials and best practices

---

## 👩‍💻 Author

**Syeda Momina**

---

##  Remarks

* This project demonstrates a **hybrid analytical workflow**, integrating graphical and programmatic tools.
* The separation into stages ensures **modularity, clarity, and reproducibility**.
* The implementation reflects **standard practices in modern single-cell transcriptomics analysis**.

---
