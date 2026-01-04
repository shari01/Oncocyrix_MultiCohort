# Oncocyrix Multicohort

Oncocyrix Multicohort is a production-grade, Scanpy-based single-cell RNA-seq analysis pipeline designed for 10x Genomics data. 

It supports single-sample and multi-cohort studies, with optional CAR-T–aware analysis, cell-type annotation, pseudobulk DESeq2, and multi-database pathway enrichment.

The pipeline is built for research reproducibility, structured outputs, and large-scale cohort integration, including pre/post, multi-patient, and multi-condition study designs.

---

## Overview

This repository provides an end-to-end single-cell RNA-seq workflow that combines:

- Standard scRNA-seq best practices (QC, clustering, DE)
- Cohort-aware integration and comparison
- Bulk-grade statistics via pseudobulk aggregation
- Immuno-oncology–focused CAR-T state modeling

It is suitable for both exploratory single-cell analysis and translational / clinical research pipelines.

---

## Key Features

### Core scRNA-seq Analysis
- 10x Genomics loader (raw `mtx/tsv/gz`)
- Robust QC & filtering (mitochondrial %, gene counts)
- Normalization, log1p, HVG selection
- Dimensionality reduction (PCA, UMAP, optional t-SNE)
- Clustering & trajectory inference (Leiden, DPT)

### Multi-Cohort & Integration
- Single-sample and multi-cohort modes
- Batch correction & integration (BBKNN)
- Group-wise comparisons (e.g., Pre vs Post, Tumor vs Normal)

### Cell-Type Annotation
- Metadata-driven annotation (if provided)
- ML-based annotation via CellTypist (optional)
- Cell-type–specific and cluster-specific marker discovery

### Differential Expression
- Single-cell DE (Scanpy, Wilcoxon)
- Group-wise DE within clusters or cell types
- Pseudobulk aggregation
- DESeq2 via rpy2 for bulk-grade inference

### Pathway Enrichment
- Multi-database enrichment:
  - GO BP / MF / CC
  - KEGG
  - Reactome
  - WikiPathways
- Publication-ready plots
- Semantic pathway deduplication (MiniLM + FAISS, optional)

---

## CAR-T Analysis (Optional)

When enabled, the pipeline becomes CAR-T aware, providing biologically interpretable functional state modeling.

### Enabling CAR-T
CAR-T analysis is controlled in code, not via CLI:

```python
DO_CART_SCORING = True
```

### Supported CAR-T States
- TStemCM_like
- TPEX_like
- TEX_terminal
- Effector_TEFF
- Proliferating_T
- Terminal_diff

### CAR-T Outputs
- CAR-T gene signature scoring per cell
- Refined CAR-T state classification
- CAR-T UMAPs and score visualizations
- Patient × phase × state × gene summaries
- Pre/Post delta tables
- CAR-T state marker genes

---
## use pip install oncocyrix_multicohort --mode multi --multi-base-dir "folder location where 10x multi samples are avb" 



## Installation

### Core installation
```bash
pip install .
```

### Full installation (recommended)
```bash
pip install ".[all]"
```

### R dependencies (for DESeq2)
Requires R ≥ 4.0 with:
```r
install.packages(c("DESeq2", "ggplot2", "pheatmap"))
```

---

## Package Structure

```text
oncocyrix-multicohort/
├── pyproject.toml
├── README.md
└── oncocyrix_multicohort/
    ├── __init__.py
    └── pipeline.py
```

### CLI entry point
```text
oncocyrix-multicohort → oncocyrix_multicohort.pipeline:main
```

---

## Usage

### Single-sample mode
```bash
oncocyrix-multicohort \
  --mode single \
  --single-10x-dir /path/to/10x/sample \
  --out-name SC_ANALYSIS_RESULTS
```

### Multi-cohort mode
```bash
oncocyrix-multicohort \
  --mode multi \
  --multi-base-dir /path/to/GSE208653_RAW \
  --out-name SC_ANALYSIS_RESULTS
```

> See `cli.md` for full CLI usage and metadata requirements.

---

## Output Structure

```text
SC_ANALYSIS_RESULTS/
├── 00_analysis_summary
├── 01_qc_and_filtering
├── 02_highly_variable_genes
├── 03_dimensionality_reduction_and_embeddings
├── 04_clustering_and_cell_states
├── 05_celltype_analysis
├── 05_CART_analysis
├── 06_groupwise_deg
├── 07_pathway_enrichment
├── 08_pseudobulk
├── 09_reference_summary
└── *.h5ad
```

---

## Requirements

- Python ≥ 3.9 (tested up to 3.12)
- R ≥ 4.0 (for DESeq2)
- Recommended RAM ≥ 32 GB

---

## Citation

Malik S.  
*Oncocyrix Multicohort: A CAR-T–aware single-cell RNA-seq analysis framework.*

---

## Author

Sheryar Malik  
Bioinformatics Scientist

---

## License

MIT License
