#' Run the full iBAQ pipeline and write results to disk
#'
#' @description
#' Convenience wrapper that calls [compute_ibaq()] with the specified
#' arguments, writes the resulting iBAQ table to a tab-separated file, and
#' invisibly returns the tibble for further use in R.
#'
#' This function is intended for straightforward, end-to-end runs where the
#' user wants to compute iBAQ values and immediately save them to disk without
#' further customisation. For finer control over individual steps, call
#' [compute_ibaq()] directly.
#'
#' @param parquet_path `character(1)`. Path to the DIA-NN `report.parquet`
#'   file. Default: `"report.parquet"` (current working directory).
#' @param fasta_path `character(1)`. Path to the reference proteome FASTA.
#'   Default: `"proteome.fasta"` (current working directory).
#' @param output_path `character(1)`. Destination path for the output TSV file.
#'   Parent directories must already exist. Default: `"iBAQ_results.tsv"`.
#' @param ... Additional arguments forwarded verbatim to [compute_ibaq()].
#'   Use these to override any default filtering, digestion, or transformation
#'   parameters (e.g. `log2_transform = FALSE`, `output_long = TRUE`).
#'
#' @return Invisibly returns the iBAQ [tibble][tibble::tibble] (same as
#'   [compute_ibaq()]). The primary side-effect is writing a TSV file to
#'   `output_path`.
#'
#' @seealso [compute_ibaq()] for full parameter documentation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # ── Minimal call (all defaults) ─────────────────────────────────────────────
#' run_ibaq_pipeline(
#'   parquet_path = "report.parquet",
#'   fasta_path   = "UP000031575.fasta"
#' )
#'
#' # ── Save to a custom path with extra parameters ─────────────────────────────
#' result <- run_ibaq_pipeline(
#'   parquet_path   = "data/report.parquet",
#'   fasta_path     = "data/proteome.fasta",
#'   output_path    = "results/iBAQ_results.tsv",
#'   log2_transform = FALSE,
#'   output_long    = TRUE
#' )
#'
#' head(result)
#' }
run_ibaq_pipeline <- function(
  parquet_path = "report.parquet",
  fasta_path   = "proteome.fasta",
  output_path  = "iBAQ_results.tsv",
  ...
) {
  result <- compute_ibaq(
    parquet_path = parquet_path,
    fasta_path   = fasta_path,
    ...
  )

  message("\n[run_ibaq_pipeline] Writing results to: ", output_path)
  readr::write_tsv(result, output_path)
  message("[run_ibaq_pipeline] Done. ", nrow(result), " proteins saved.")

  invisible(result)
}
