# CLI Usage Guide (cli.md)

This document explains how to run the Scanpy-based single-cell RNA-seq
pipeline via the command line (CLI), with or without CAR-T analysis.

The pipeline supports:
- Single 10x dataset analysis
- Multi-sample / multi-cohort analysis
- Optional CAR-T functional state scoring
- Optional pathway enrichment and pseudobulk DESeq2

---

## 1. Supported CLI Arguments

The script supports the following CLI arguments:

- `--mode {single,multi}`
- `--single-10x-dir PATH`
- `--multi-base-dir PATH`
- `--out-name NAME`
- `--no-pathway-clustering`

### Argument descriptions

- `--mode`
  - `single`: analyze one 10x dataset
  - `multi`: analyze multiple samples and perform a combined analysis

- `--single-10x-dir`
  - Path to a single 10x Genomics folder
  - Required when `--mode single`

- `--multi-base-dir`
  - Base directory containing group/sample folders
  - Required when `--mode multi`

- `--out-name`
  - Custom name for the output directory
  - Default: `SC_ANALYSIS_RESULTS`

- `--no-pathway-clustering`
  - Disables pathway enrichment (Enrichr / gseapy)
  - Useful for quick exploratory runs

---

## 2. CAR-T Scoring (IMPORTANT)

There is **NO CLI flag** for CAR-T analysis.

CAR-T scoring is controlled by a hardcoded variable in the script:

```python
DO_CART_SCORING = True
