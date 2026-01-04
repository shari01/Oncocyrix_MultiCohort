<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
</head>
<body>

<h1>Oncocyrix Multicohort</h1>

<p>
<strong>Oncocyrix Multicohort</strong> is a production-grade single-cell RNA-seq analysis pipeline
designed for large-scale <strong>10x Genomics</strong> datasets, supporting both
single-sample and multi-cohort study designs with optional <strong>CAR-T–aware analysis</strong>.
</p>

<p>
The pipeline is built for reproducibility, structured outputs, and translational research,
covering cohort integration, cell-type annotation, differential expression,
pseudobulk statistics, and multi-database pathway enrichment.
</p>

<hr/>

<h2>🔗 PyPI Package</h2>
<p>
<strong>Project page:</strong>
<a href="https://pypi.org/project/oncocyrix-multicohort/" target="_blank">
https://pypi.org/project/oncocyrix-multicohort/
</a>
</p>

<p>
<strong>Install via pip:</strong>
</p>

<pre><code>pip install oncocyrix-multicohort</code></pre>

<hr/>

<h2>Overview</h2>

<p>
This repository provides an end-to-end single-cell RNA-seq workflow that combines:
</p>

<ul>
  <li>Standard scRNA-seq best practices (QC, clustering, DE)</li>
  <li>Multi-sample and multi-cohort integration</li>
  <li>Bulk-grade statistics via pseudobulk aggregation</li>
  <li>Immuno-oncology–focused CAR-T functional state modeling</li>
</ul>

<p>
It is suitable for both exploratory single-cell analysis and
translational / clinical research pipelines.
</p>

<hr/>

<h2>Key Features</h2>

<h3>Core scRNA-seq Analysis</h3>
<ul>
  <li>10x Genomics loader (matrix.mtx / barcodes / features)</li>
  <li>Robust QC & filtering (mitochondrial %, gene counts)</li>
  <li>Normalization, log1p, HVG selection</li>
  <li>PCA, UMAP, optional t-SNE</li>
  <li>Clustering & trajectory inference (Leiden, DPT)</li>
</ul>

<h3>Multi-Cohort & Integration</h3>
<ul>
  <li>Single-sample and multi-cohort modes</li>
  <li>Batch correction and integration (BBKNN)</li>
  <li>Group-wise comparisons (Pre/Post, Tumor/Normal)</li>
</ul>

<h3>Cell-Type Annotation</h3>
<ul>
  <li>Metadata-driven annotation (if provided)</li>
  <li>Machine-learning annotation via CellTypist (optional)</li>
  <li>Cell-type–specific and cluster-specific marker discovery</li>
</ul>

<h3>Differential Expression</h3>
<ul>
  <li>Single-cell DE (Wilcoxon, Scanpy)</li>
  <li>Group-wise DE within clusters or cell types</li>
  <li>Pseudobulk aggregation</li>
  <li>DESeq2 integration via rpy2</li>
</ul>

<h3>Pathway Enrichment</h3>
<ul>
  <li>GO (BP / MF / CC)</li>
  <li>KEGG</li>
  <li>Reactome</li>
  <li>WikiPathways</li>
  <li>Semantic pathway deduplication (MiniLM + FAISS, optional)</li>
</ul>

<hr/>

<h2>CAR-T Analysis (Optional)</h2>

<p>
When enabled, the pipeline performs CAR-T–aware functional state modeling,
designed for immuno-oncology and cell therapy studies.
</p>

<h3>Supported CAR-T States</h3>
<ul>
  <li>TStemCM_like</li>
  <li>TPEX_like</li>
  <li>TEX_terminal</li>
  <li>Effector_TEFF</li>
  <li>Proliferating_T</li>
  <li>Terminal_diff</li>
</ul>

<h3>CAR-T Outputs</h3>
<ul>
  <li>Signature gene scoring per cell</li>
  <li>Refined CAR-T state classification</li>
  <li>CAR-T UMAP visualizations</li>
  <li>Patient × phase × state × gene summaries</li>
  <li>Pre/Post delta tables</li>
  <li>CAR-T state marker genes</li>
</ul>

<hr/>

<h2>Usage</h2>

<h3>Multi-cohort (with CAR-T)</h3>
<pre><code>
oncocyrix-multicohort --mode multi \
--multi-base-dir "C:\Users\shery\Downloads\oncocyrix_multicohort\GSE208653_RAW-v2" \
--do-cart
</code></pre>

<h3>Multi-cohort (without CAR-T)</h3>
<pre><code>
oncocyrix-multicohort --mode multi \
--multi-base-dir "C:\Users\shery\Downloads\oncocyrix_multicohort\GSE208653_RAW-v2"
</code></pre>

<hr/>

<h2>Output Structure</h2>
<h2>Output-results</h2>
https://drive.google.com/drive/folders/1hURfwj2z8Vqqn0QIdR1omaxPqWwfSDhs?usp=sharing

<pre><code>
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
</code></pre>

<hr/>

<h2>Requirements</h2>
<ul>
  <li>Python ≥ 3.9 (tested up to 3.12)</li>
  <li>R ≥ 4.0 (for DESeq2)</li>
  <li>Recommended RAM ≥ 32 GB</li>
</ul>

<hr/>

<h2>Author</h2>
<p>
<strong>Sheryar Malik</strong><br/>
Bioinformatics Scientist
</p>

<hr/>

<h2>Citation</h2>
<p>
Malik S.<br/>
<i>Oncocyrix Multicohort: A CAR-T–aware single-cell RNA-seq analysis framework.</i>
</p>

<hr/>

<h2>License</h2>
<p>MIT License</p>

</body>
</html>
