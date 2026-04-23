# 🧬 Stage 1: Pre-processing of 10X Single-Cell RNA Datasets

> **Platform:** Galaxy Training Network (GTN)  
> **Folder:** `1_preprocessing_10X_galaxy/`  
> **Input:** Raw paired-end FASTQ files (10X Genomics Chromium)  
> **Output:** Gene-cell count matrix (`.mtx`), QC reports, AnnData (`.h5ad`)

---

## 📖 Table of Contents

1. [Background & Motivation](#-background--motivation)
2. [10X Genomics Data Structure](#-10x-genomics-data-structure)
3. [Key Concepts](#-key-concepts)
4. [Pipeline Overview](#-pipeline-overview)
5. [Step-by-Step Workflow](#-step-by-step-workflow)
   - [Step 1 — Data Input (FASTQ)](#step-1--data-input-fastq)
   - [Step 2 — Alignment with STARsolo](#step-2--alignment-with-starsolo)
   - [Step 3 — Cell Calling with DropletUtils](#step-3--cell-calling-with-dropletutils)
   - [Step 4 — Quality Control with MultiQC](#step-4--quality-control-with-multiqc)
   - [Step 5 — Count Matrix Export](#step-5--count-matrix-export)
6. [Output Files Explained](#-output-files-explained)
7. [Quality Metrics Interpretation](#-quality-metrics-interpretation)
8. [Common Issues & Troubleshooting](#-common-issues--troubleshooting)
9. [Galaxy Workflow Tips](#-galaxy-workflow-tips)

---

##  Background & Motivation

Single-cell RNA sequencing (scRNA-seq) using the **10X Genomics Chromium** platform has become the gold standard for profiling transcriptomes at single-cell resolution. Each cell is captured in a gel bead-in-emulsion (GEM) droplet, tagged with a unique **cell barcode**, and all RNA molecules are tagged with **Unique Molecular Identifiers (UMIs)**.

Before any biological analysis can begin, raw sequencing data must be:

1. **Aligned** to a reference genome/transcriptome to identify which gene each read came from
2. **Demultiplexed** by cell barcode to assign reads to individual cells
3. **Deduplicated** using UMIs to remove PCR amplification artifacts
4. **Quality-controlled** to ensure data integrity before downstream analysis

This stage handles all of the above using **Galaxy** — a web-based, graphical bioinformatics platform — making the preprocessing accessible without requiring command-line expertise.

---

## 🧬 10X Genomics Data Structure

### FASTQ File Pairs

10X Genomics data always comes as **paired-end reads** split across two files:

| File | Read | Content | Length |
|------|------|---------|--------|
| `*_R1.fastq.gz` | Read 1 | Cell barcode (16 bp) + UMI (12 bp) | 28 bp |
| `*_R2.fastq.gz` | Read 2 | cDNA sequence (gene-derived) | 90–150 bp |

- **Read 1** is not aligned — it carries the cell identity information
- **Read 2** is the biologically informative read that maps to the transcriptome

### The 10X Chromium Chemistry

```
Cell Capture → GEM Formation → Barcoded Reverse Transcription
       ↓
 [Cell Barcode 16bp] + [UMI 12bp] + [polyT] + [cDNA]
       ↓
    PCR Amplification → Library Prep → Sequencing
```

Each unique **cell barcode** identifies one cell. Each unique **UMI** within a barcode-gene combination represents one original mRNA molecule.

---

## 🔑 Key Concepts

### STARsolo vs Cell Ranger

| Feature | STARsolo | Cell Ranger |
|---------|----------|-------------|
| **Developer** | Alexander Dobin (STAR) | 10X Genomics |
| **Speed** |  Very fast |  Slower |
| **License** | Open source (MIT) | Proprietary (requires 10X account) |
| **Galaxy support** |  Yes |  Limited |
| **Output format** | `.mtx`, `.h5ad` | `.h5`, filtered/raw folders |
| **Accuracy** | Comparable | Comparable |

STARsolo is the preferred open-source alternative to Cell Ranger and is fully integrated into Galaxy.

### Count Matrix Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| **Market Exchange** | `.mtx` | Sparse matrix format; 3 files: matrix, barcodes, features |
| **HDF5** | `.h5` / `.h5ad` | Binary format; single file; faster I/O; AnnData native |
| **CSV/TSV** | `.csv` | Dense; not recommended for large datasets |

For this pipeline, we use `.mtx` output from Galaxy and import it into AnnData (`.h5ad`) for downstream analysis.

---

## 🗺️ Pipeline Overview

```
Paired FASTQ (R1 + R2)
        │
        ▼
  ┌─────────────┐
  │  STARsolo   │  ← Barcode extraction, UMI deduplication, genome alignment
  └──────┬──────┘
         │ Raw count matrix (.mtx) + logs
         ▼
  ┌─────────────────┐
  │  DropletUtils   │  ← Knee-point detection, empty droplet filtering
  └──────┬──────────┘
         │ Filtered cell barcodes
         ▼
  ┌─────────────┐
  │   MultiQC   │  ← Aggregated QC report across all samples
  └──────┬──────┘
         │ QC-verified count matrix
         ▼
  Gene-Cell Count Matrix (.mtx / .h5ad)
```

---

## 🔬 Step-by-Step Workflow

### Step 1 — Data Input (FASTQ)

**What to do in Galaxy:**
1. Navigate to **Upload Data** in Galaxy
2. Upload your `*_R1.fastq.gz` and `*_R2.fastq.gz` files
3. Set the datatype to `fastqsanger.gz`

**Important notes:**
- Ensure R1 and R2 files are correctly paired (same sample name prefix)
- For multiple samples, upload all pairs; Galaxy supports batch processing
- Verify file integrity with MD5 checksums if provided by your sequencing facility

**File size expectations:**
- A typical 10X dataset (~5,000 cells, ~50M reads) will be ~5–15 GB per FASTQ file

---

### Step 2 — Alignment with STARsolo

**Tool:** RNA STARsolo  
**Galaxy tool ID:** `rna_starsolo`

#### Key Parameters

| Parameter | Value | Explanation |
|-----------|-------|-------------|
| **Chemistry** | 10x Chromium v3 | Specifies barcode + UMI lengths |
| **Genome** | GRCh38 / mm10 | Must match your organism |
| **Cell barcode whitelist** | 10X 3M-february-2018.txt | Known valid barcodes |
| **Solo type** | CB_UMI_Simple | Standard 10X mode |
| **Output** | Gene counts (raw) | Includes all detected barcodes |

#### What STARsolo Does

1. **Barcode extraction:** Reads R1 to identify cell barcodes; corrects 1-bp mismatches against whitelist
2. **UMI extraction:** Reads remaining bases of R1 as UMI
3. **Alignment:** Maps R2 to genome using STAR's splice-aware alignment
4. **Gene assignment:** Assigns each aligned read to a gene using genome annotation (GTF)
5. **UMI deduplication:** Collapses PCR duplicates using UMI sequences
6. **Matrix generation:** Produces a barcode × gene count matrix

#### Output Files from STARsolo

```
STARsolo_output/
├── matrix.mtx          ← Sparse count matrix (rows=genes, cols=barcodes)
├── barcodes.tsv        ← List of all detected cell barcodes
├── features.tsv        ← Gene IDs and names
└── Log.final.out       ← Alignment statistics summary
```

**Alignment log file** (`RNA STARSolo on dataset 1-6_ log].txt`) contains:
- Total reads processed
- Uniquely mapped reads (%)
- Multi-mapped reads (%)
- Reads mapped to too many loci
- Unmapped reads (%)

>  **Good alignment:** Uniquely mapped reads > 60% is acceptable; > 75% is ideal for human/mouse data.

---

### Step 3 — Cell Calling with DropletUtils

**Tool:** DropletUtils  
**Purpose:** Separate real cells from empty droplets

#### The Empty Droplet Problem

In a 10X experiment, the vast majority of droplets do not capture a cell — they are empty but may still contain ambient RNA. The raw count matrix includes **all barcodes** detected, including:
-  Real cells (high UMI count)
-  Empty droplets (very low UMI count from ambient RNA)
-  Doublets (two cells in one droplet — handled later by Scrublet/DoubletFinder)

#### The Knee Plot
<img width="480" height="480" alt="Galaxy28- DropletUtils Plot on dataset 14-16" src="https://github.com/user-attachments/assets/1596cf43-1b4f-401d-91f8-0d7bc229cdd8" />


**How to read this plot:**
- **X-axis:** Barcode rank (sorted by decreasing UMI count)
- **Y-axis:** Total UMI count per barcode (log scale)
- **Knee point:** Sharp inflection — separates real cells (plateau region) from empty droplets (steep drop)
- **Inflection point:** Second change in slope; alternative threshold

**Interpretation of this dataset:**
- Barcodes to the left of the knee (high UMI counts) → **real cells**
- Barcodes to the right (low UMI counts) → **empty droplets / debris**
- A clean, sharp knee indicates high-quality data with good cell capture

#### DropletUtils Algorithm

DropletUtils uses two complementary approaches:

1. **`barcodeRanks()`** — Generates the knee plot; identifies knee/inflection thresholds
2. **`emptyDrops()`** — Statistical test comparing each barcode's expression profile to the ambient RNA profile; retains barcodes significantly different from ambient (FDR < 0.001)

The `emptyDrops()` approach is more sensitive than simple UMI thresholding and is recommended for datasets with continuous UMI distributions.

---

### Step 4 — Quality Control with MultiQC

**Tool:** MultiQC  
**Input:** STARsolo log files + DropletUtils outputs  
**Output:** `MultiQC on dataset 13_ Stats].tabular`

MultiQC aggregates quality metrics from multiple tools and samples into a single comprehensive report.

#### Key Metrics to Evaluate

| Metric | Threshold | Meaning |
|--------|-----------|---------|
| **% Uniquely Mapped** | > 60% | Reads aligning to exactly one genomic location |
| **% Mapped to Multiple Loci** | < 20% | Reads mapping to repetitive regions |
| **% Unmapped (too short)** | < 10% | Read quality/trimming indicator |
| **Duplication Rate** | Variable | Expected to be high in scRNA-seq due to UMI deduplication |
| **Average Input Read Length** | ~90–150 bp | Should match sequencing parameters |

#### Interpreting the MultiQC Stats Table

The `.tabular` file contains per-sample statistics. Key columns:

```
Sample | Total Reads | Uniquely Mapped % | Mismatch Rate | ...
```

>  **Warning signs:** Uniquely mapped reads < 50%, high mismatch rates (> 5%), or very low read counts per sample indicate a quality problem and should be investigated before proceeding.

---

### Step 5 — Count Matrix Export

After QC validation, export the following files from Galaxy:

1. **Raw count matrix:** `RNA STARSolo on dataset 1-6_ Matrix Gene Counts raw].mtx`
2. **Barcodes list:** `barcodes.tsv`
3. **Features/genes list:** `features.tsv`

These three files together constitute the complete count matrix and are the inputs for **Stage 2 (Scanpy Analysis)**.

The `.mtx` format is a sparse matrix representation:
```
%% MatrixMarket matrix coordinate integer general
%%
33538 6794 8517005      ← (genes, barcodes, non-zero entries)
5       1   1            ← (gene_index, barcode_index, count)
17      1   1
...
```

---

##  Output Files Explained

| File | Format | Description | Used In |
|------|--------|-------------|---------|
| `*.mtx` | Sparse matrix | Raw gene × cell count matrix | Stage 2 (Scanpy) |
| `barcodes.tsv` | Tab-separated | Cell barcode sequences | Stage 2 (Scanpy) |
| `features.tsv` | Tab-separated | Ensembl gene IDs + gene names | Stage 2 (Scanpy) |
| `Log.final.out` | Plain text | STARsolo alignment statistics | QC review |
| `MultiQC_Stats.tabular` | Tab-separated | Aggregated QC metrics table | QC review |
| `DropletUtils_plot.png` | PNG image | Knee/barcode rank plot | QC review |

---

## 📈 Quality Metrics Interpretation

### What Good Data Looks Like

| Metric | Good | Acceptable | Concerning |
|--------|------|------------|------------|
| Uniquely mapped reads | > 75% | 60–75% | < 60% |
| Cells detected | 1,000–10,000 | 500–15,000 | < 200 or > 20,000 |
| Median genes/cell | 1,500–4,000 | 800–6,000 | < 500 |
| Median UMIs/cell | 3,000–15,000 | 1,000–30,000 | < 500 |
| % Mitochondrial genes | < 10–15% | 15–25% | > 25% (cell stress/death) |

### Knee Plot Quality Indicators

- **Sharp knee** → Clean experiment, clear cell/empty droplet separation 
- **Gradual slope** → Possible multiplet contamination or low-quality cells 
- **No knee visible** → Potential sequencing or library prep failure 

---

##  Common Issues & Troubleshooting

### Issue: Very low mapping rate (< 50%)

**Possible causes:**
- Wrong genome/species selected
- Poor FASTQ quality (check with FastQC)
- Contamination with other organisms
- Using wrong STARsolo chemistry (v2 vs v3)

**Solution:** Re-run with correct genome annotation; verify chemistry version matches your 10X kit

---

### Issue: No clear knee in DropletUtils plot

**Possible causes:**
- Very low cell capture efficiency
- High ambient RNA contamination
- Wrong UMI length specified

**Solution:** Adjust `emptyDrops()` FDR threshold; inspect the UMI count distribution; consult library QC metrics from Bioanalyzer/TapeStation

---

### Issue: Fewer cells than expected

**Possible causes:**
- Overly stringent DropletUtils filtering
- Low-quality sample with many dying cells
- Cell viability < 80% at time of capture

**Solution:** Relax FDR threshold in `emptyDrops()`; review wet-lab QC before re-sequencing

---

## 🌐 Galaxy Workflow Tips

1. **Use the GTN tutorial** as your guide: [Single Cell 1M Tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/)
2. **Save your history** after each major step — Galaxy histories serve as an electronic lab notebook
3. **Use workflow extraction** to convert your history into a reusable `.ga` workflow file
4. **Set meaningful names** for each dataset in your history to avoid confusion
5. **Use Galaxy's dataset collections** for processing multiple samples in parallel

---

⬅️ [Back to Main README](../README.md) | ➡️ [Next: Scanpy Analysis](../2_scanpy_analysis/README.md)
