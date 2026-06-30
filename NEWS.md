# NEWS.md

# ibaq4diann 0.1.0

## New features

* `compute_ibaq()`: Main function to compute iBAQ values from DIA-NN
  `report.parquet` and a reference proteome FASTA.
* `run_ibaq_pipeline()`: Convenience wrapper that computes iBAQ values and
  writes results to a TSV file.
* `filter_diann_report()`: Quality-control filtering of DIA-NN report tibbles
  with configurable FDR, proteotypicity, and LFQ quality thresholds.
* `plot_ibaq_rank()`: Ranked protein abundance plot per sample (Whittaker plot).
* `plot_ibaq_distribution()`: Violin and/or box plot of iBAQ value
  distributions across samples.
* `plot_ibaq_correlation()`: Sample-to-sample Pearson/Spearman correlation
  heatmap.
* `plot_ibaq_dynamic_range()`: Dynamic range visualisation with per-sample
  min–median–max segments and protein count annotations.
