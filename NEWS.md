# NEWS.md

# ibaq4diann 0.2.0

## New features

* `compute_ibaq()` gains a validated `id_col` argument to choose which DIA-NN
  identifier column defines a protein: `"Protein.Group"` (default),
  `"Protein.Ids"`, `"Protein.Names"` or `"Genes"`. The FASTA is matched by
  accession, entry name or `GN=` gene symbol accordingly. Multi-member
  identifiers (`;`-separated, including `cRAP-` contaminants) are matched by
  their first member found in the FASTA. Messages report how the number of
  proteins changes with the chosen identifier (e.g. protein groups merged
  under one gene, rows with empty identifiers dropped).
* `parse_fasta()` gains an `id_type` argument (`"accession"`, `"entry_name"`,
  `"gene"`), dropping entries without the requested identifier and keeping
  the first occurrence of duplicates.

## Improvements

* `plot_ibaq_pie()` redesigned for legibility: slices are labelled beside the
  pie with leader lines (`"NAME (x.x%)"`), the "Other" slice is grey and
  labelled inside when large, and the redundant legend is hidden by default
  (`show_legend`). New `min_label_pct` argument.
* `plot_ibaq_rank()` places top-N labels in a column beside the curve instead
  of stacking them at the panel edge.
* Plotting functions now recognise all identifier columns as metadata, so
  results computed with any `id_col` convert correctly from wide to long.

## Bug fixes

* Output annotation columns are now collapsed per `protein_id` instead of
  taking the first row, so no annotation is silently lost when several
  report rows map to one identifier.

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
