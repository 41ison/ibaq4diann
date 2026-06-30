# ibaq4diann

<!-- badges: start -->

[![R-CMD-check](https://github.com/41ison/ibaq4diann/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/41ison/ibaq4diann/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**ibaq4diann** computes **intensity-Based Absolute Quantification
(iBAQ)** values from [DIA-NN](https://github.com/vdemichev/DiaNN) search
results (`report.parquet`) and a reference proteome FASTA. It implements
the iBAQ method of [Schwanhäusser *et al.*
(2011)](https://doi.org/10.1038/nature10098) and provides plotting
utilities for exploratory data analysis.

## Overview

    iBAQ = Σ(precursor intensities per protein) / (# theoretical tryptic peptides)

The package provides:

| Function                    | Description                                                |
|-----------------------------|------------------------------------------------------------|
| `compute_ibaq()`            | **Main function** — compute iBAQ values from DIA-NN output |
| `run_ibaq_pipeline()`       | Convenience wrapper that also writes results to TSV        |
| `filter_diann_report()`     | Quality-control filtering of DIA-NN report tibbles         |
| `plot_ibaq_rank()`          | Ranked protein abundance plot per sample                   |
| `plot_ibaq_distribution()`  | Violin/box plot of iBAQ distributions                      |
| `plot_ibaq_correlation()`   | Sample-to-sample correlation heatmap                       |
| `plot_ibaq_dynamic_range()` | Dynamic range visualisation per sample                     |

## Installation

``` r
# Install from GitHub (requires remotes)
# install.packages("remotes")
remotes::install_github("41ison/ibaq4diann")
```

> **Note**: `ibaq4diann` depends on `Biostrings` from Bioconductor. If
> it is not already installed, run:
>
> ``` r
> if (!requireNamespace("BiocManager", quietly = TRUE))
>   install.packages("BiocManager")
> BiocManager::install("Biostrings")
> ```

## Quick start

``` r
library(ibaq4diann)

# Compute iBAQ values (wide format, log2-transformed by default)
ibaq <- compute_ibaq(
  parquet_path = "report.parquet",
  fasta_path   = "proteome.fasta"
)

head(ibaq)

# Save to TSV
readr::write_tsv(ibaq, "iBAQ_results.tsv")
```

## Visualisation

``` r
# Ranked protein abundance per sample
plot_ibaq_rank(ibaq, top_n_labels = 10)

# iBAQ value distributions (violin + box)
plot_ibaq_distribution(ibaq, plot_type = "both")

# Sample correlation heatmap
plot_ibaq_correlation(ibaq)

# Dynamic range per sample
plot_ibaq_dynamic_range(ibaq)
```

## Parameters at a glance

``` r
compute_ibaq(
  parquet_path          = "report.parquet",       # DIA-NN output
  fasta_path            = "proteome.fasta",        # Reference proteome
  intensity_col         = "Precursor.Normalised",  # or "Precursor.Quantity"
  q_value_cutoff        = 0.01,   # Precursor FDR
  pg_q_value_cutoff     = 0.01,   # Protein group FDR
  lib_pg_q_value_cutoff = 0.01,   # Library PG FDR
  proteotypic_only      = TRUE,   # Unique peptides only
  lfq_quality_cutoff    = 0.50,   # LFQ quality threshold
  min_peptide_len       = 6,      # In-silico digest: min peptide length
  max_peptide_len       = 30,     # In-silico digest: max peptide length
  max_missed_cleavages  = 0,      # Missed cleavages allowed
  log2_transform        = TRUE,   # Apply log2(iBAQ + 1)
  output_long           = FALSE   # Wide (default) or long format
)
```