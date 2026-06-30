#' Apply quality-control filters to a DIA-NN report
#'
#' @description
#' Filters a DIA-NN report data frame (loaded from `report.parquet`) according
#' to standard quality thresholds used in proteomics data analysis. The default
#' thresholds follow DIA-NN best-practice recommendations and the
#' iBAQ-for-DIANN reference implementation.
#'
#' Filters are applied only when the corresponding column exists in `report`,
#' so the function is robust to different DIA-NN output versions.
#'
#' @param report A `data.frame` or `tibble` produced by
#'   [arrow::read_parquet()] on a DIA-NN `report.parquet` file.
#' @param q_value_cutoff `numeric(1)` or `NULL`. Maximum allowed precursor-level
#'   false discovery rate (`Q.Value` column). Set to `NULL` to skip this filter.
#'   Default: `0.01` (1 % FDR).
#' @param pg_q_value_cutoff `numeric(1)` or `NULL`. Maximum allowed
#'   protein-group-level FDR (`PG.Q.Value` column). Set to `NULL` to skip.
#'   Default: `0.01`.
#' @param lib_pg_q_value_cutoff `numeric(1)` or `NULL`. Maximum allowed
#'   library protein-group-level FDR (`Lib.PG.Q.Value` column). Set to `NULL`
#'   to skip. Default: `0.01`.
#' @param proteotypic_only `logical(1)`. If `TRUE` (default), retains only rows
#'   where `Proteotypic == 1`, i.e. unique peptides unambiguously mapped to a
#'   single protein. Recommended for accurate iBAQ quantification.
#' @param lfq_quality_cutoff `numeric(1)` or `NULL`. Minimum LFQ quality score
#'   (`PG.MaxLFQ.Quality` column). Set to `NULL` or `0` to skip. Default: `0.5`.
#'
#' @return A filtered `tibble` with the same columns as `report` but fewer rows.
#'   A summary of rows removed is printed to the console via [message()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- arrow::read_parquet("report.parquet")
#'
#' # Default filters (1 % FDR, proteotypic, LFQ quality >= 0.5)
#' filtered <- filter_diann_report(report)
#'
#' # Relaxed: no LFQ quality filter, allow non-proteotypic peptides
#' filtered_relaxed <- filter_diann_report(
#'   report,
#'   lfq_quality_cutoff = NULL,
#'   proteotypic_only = FALSE
#' )
#' }
filter_diann_report <- function(
  report,
  q_value_cutoff = 0.01,
  pg_q_value_cutoff = 0.01,
  lib_pg_q_value_cutoff = 0.01,
  proteotypic_only = TRUE,
  lfq_quality_cutoff = 0.50
) {
  n_before <- nrow(report)
  message("[filter_diann_report] Rows before filtering: ", n_before)

  result <- report

  if (!is.null(q_value_cutoff) && "Q.Value" %in% names(result)) {
    result <- dplyr::filter(result, .data$Q.Value <= q_value_cutoff)
  }
  if (!is.null(pg_q_value_cutoff) && "PG.Q.Value" %in% names(result)) {
    result <- dplyr::filter(result, .data$PG.Q.Value <= pg_q_value_cutoff)
  }
  if (!is.null(lib_pg_q_value_cutoff) && "Lib.PG.Q.Value" %in% names(result)) {
    result <- dplyr::filter(result, .data$Lib.PG.Q.Value <= lib_pg_q_value_cutoff)
  }
  if (isTRUE(proteotypic_only) && "Proteotypic" %in% names(result)) {
    result <- dplyr::filter(result, .data$Proteotypic == 1)
  }
  if (
    !is.null(lfq_quality_cutoff) &&
    lfq_quality_cutoff > 0 &&
    "PG.MaxLFQ.Quality" %in% names(result)
  ) {
    result <- dplyr::filter(result, .data$PG.MaxLFQ.Quality >= lfq_quality_cutoff)
  }

  message(
    "[filter_diann_report] Rows after filtering:  ",
    nrow(result),
    " (removed ",
    n_before - nrow(result),
    " rows)"
  )
  result
}
