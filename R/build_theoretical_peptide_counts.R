#' Build a table of theoretical tryptic peptide counts for a proteome
#'
#' @description
#' Applies [count_theoretical_peptides()] to every sequence in a named
#' character vector (as returned by [parse_fasta()]) and returns a tidy tibble
#' with one row per protein. This table provides the denominators for the iBAQ
#' calculation.
#'
#' @param sequences Named `character` vector. Names are protein accessions,
#'   values are amino-acid sequences. Typically the output of [parse_fasta()].
#' @param protease_regex `character(1)`. Perl-compatible regex defining the
#'   cleavage rule. Default: `"(?<=[KR])(?!P)"` (trypsin, no cut before Pro).
#' @param min_len `integer(1)`. Minimum peptide length to count. Default: `6`.
#' @param max_len `integer(1)`. Maximum peptide length to count. Default: `30`.
#' @param max_missed `integer(1)`. Maximum missed cleavages. Default: `0`.
#'
#' @return A [tibble][tibble::tibble] with two columns:
#'   \describe{
#'     \item{`protein_id`}{`character` — protein accession.}
#'     \item{`n_theoretical_peptides`}{`integer` — number of theoretically
#'       observable tryptic peptides within the specified length window.}
#'   }
#'
#' @seealso [parse_fasta()], [compute_ibaq()]
#' @keywords internal
build_theoretical_peptide_counts <- function(
  sequences,
  protease_regex = "(?<=[KR])(?!P)",
  min_len = 6,
  max_len = 30,
  max_missed = 0
) {
  message(
    "[build_theoretical_peptide_counts] Digesting ",
    length(sequences),
    " sequences in silico ..."
  )

  counts <- vapply(
    sequences,
    count_theoretical_peptides,
    FUN.VALUE = integer(1),
    protease_regex = protease_regex,
    min_len = min_len,
    max_len = max_len,
    max_missed = max_missed
  )

  tibble::tibble(
    protein_id = names(counts),
    n_theoretical_peptides = as.integer(counts)
  )
}
