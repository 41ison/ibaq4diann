#' DIA-NN protein identifier columns supported by `compute_ibaq()`
#'
#' Maps each supported `id_col` to the FASTA identifier type (`id_type` in
#' [parse_fasta()]) that it must be matched against.
#'
#' @keywords internal
#' @noRd
diann_id_cols <- c(
  Protein.Group = "accession",
  Protein.Ids = "accession",
  Protein.Names = "entry_name",
  Genes = "gene"
)

#' Non-sample (metadata) columns that can appear in a `compute_ibaq()` result
#'
#' Used by the plotting functions to tell annotation columns apart from
#' sample columns when converting a wide-format iBAQ table to long format.
#'
#' @return `character` vector of column names.
#' @keywords internal
#' @noRd
ibaq_meta_cols <- function() {
  c(
    "protein_id",
    names(diann_id_cols),
    "n_theoretical_peptides"
  )
}

#' Convert a wide- or long-format iBAQ table to long format
#'
#' @param ibaq A `data.frame` returned by [compute_ibaq()] (wide or long).
#' @param samples `character` or `NULL`. Samples to keep.
#' @return A long-format tibble with columns `sample` and `iBAQ` plus any
#'   metadata columns present in `ibaq`.
#' @keywords internal
#' @noRd
ibaq_to_long <- function(ibaq, samples = NULL) {
  meta_present <- dplyr::intersect(ibaq_meta_cols(), names(ibaq))

  if (!"sample" %in% names(ibaq)) {
    sample_cols <- dplyr::setdiff(names(ibaq), meta_present)
    if (!is.null(samples)) {
      sample_cols <- dplyr::intersect(sample_cols, samples)
    }
    tidyr::pivot_longer(
      ibaq,
      cols = dplyr::all_of(sample_cols),
      names_to = "sample",
      values_to = "iBAQ"
    )
  } else {
    if (!is.null(samples)) {
      ibaq <- dplyr::filter(ibaq, .data$sample %in% samples)
    }
    ibaq
  }
}

#' Match (possibly multi-member) report identifiers to theoretical peptide counts
#'
#' @description
#' DIA-NN identifier columns may contain several members separated by `;`
#' (e.g. `"P0DP23;P0DP24;P0DP25"`, or `"P04040;cRAP-P04040"` when a
#' contaminant database was appended). Each identifier is split on `;` and
#' its members are looked up in `theo_peptides` in order; the theoretical
#' peptide count of the **first member found** (normally the leading protein
#' of the group) is used as the iBAQ denominator.
#'
#' @param ids `character`. Unique identifiers from the report `id_col`.
#' @param theo_peptides Tibble with columns `protein_id` and
#'   `n_theoretical_peptides`, as returned by
#'   [build_theoretical_peptide_counts()].
#'
#' @return A tibble with one row per element of `ids` and columns:
#'   \describe{
#'     \item{`protein_id`}{The original identifier.}
#'     \item{`fasta_id`}{The member that matched the FASTA (`NA` if none).}
#'     \item{`n_theoretical_peptides`}{Count for `fasta_id` (`NA` if none).}
#'   }
#' @keywords internal
#' @noRd
match_theoretical_peptides <- function(ids, theo_peptides) {
  counts <- stats::setNames(
    theo_peptides$n_theoretical_peptides,
    theo_peptides$protein_id
  )

  members <- strsplit(ids, ";", fixed = TRUE)
  fasta_id <- vapply(
    members,
    function(m) {
      hit <- m[m %in% names(counts)]
      if (length(hit) == 0) NA_character_ else hit[[1]]
    },
    FUN.VALUE = character(1)
  )

  tibble::tibble(
    protein_id = ids,
    fasta_id = fasta_id,
    n_theoretical_peptides = unname(counts[fasta_id])
  )
}
