#' Count theoretically observable tryptic peptides for a protein sequence
#'
#' @description
#' Performs an *in-silico* digest of a single amino-acid sequence using a
#' configurable cleavage rule (default: trypsin), then counts the number of
#' **unique** fully-cleaved peptides whose length falls within
#' `[min_len, max_len]`. This count is the denominator in the iBAQ formula.
#'
#' The default trypsin rule cuts after lysine (K) or arginine (R) **unless**
#' the next residue is proline (P), i.e., the regex look-behind/ahead
#' `(?<=[KR])(?!P)`.
#'
#' @param sequence `character(1)`. A single amino-acid sequence as a plain
#'   character string (e.g. `"MKTIIALSYIFCLVFA"`). `NA` or an empty string
#'   returns `NA_integer_`.
#' @param protease_regex `character(1)`. A Perl-compatible regular expression
#'   that matches **cut sites** (positions between residues where cleavage
#'   occurs). Default: `"(?<=[KR])(?!P)"` (trypsin, no-cleavage before Pro).
#' @param min_len `integer(1)`. Minimum peptide length (number of residues) to
#'   include in the count. Default: `6`.
#' @param max_len `integer(1)`. Maximum peptide length (number of residues) to
#'   include in the count. Default: `30`.
#' @param max_missed `integer(1)`. Maximum number of missed cleavages allowed.
#'   `0` = fully cleaved peptides only. Default: `0`.
#'
#' @return `integer(1)`. The number of unique theoretically observable peptides
#'   within `[min_len, max_len]` after *in-silico* digestion. Returns
#'   `NA_integer_` when `sequence` is `NA` or empty.
#'
#' @seealso [build_theoretical_peptide_counts()], [compute_ibaq()]
#' @keywords internal
count_theoretical_peptides <- function(
  sequence,
  protease_regex = "(?<=[KR])(?!P)",
  min_len = 6,
  max_len = 30,
  max_missed = 0
) {
  if (is.na(sequence) || nchar(sequence) == 0) {
    return(NA_integer_)
  }

  fragments <- unlist(strsplit(sequence, protease_regex, perl = TRUE))
  fragments <- fragments[nchar(fragments) > 0]

  n <- length(fragments)
  if (n == 0) {
    return(0L)
  }

  peptides <- character(0)
  for (mc in 0:max_missed) {
    if (mc == 0) {
      peptides <- c(peptides, fragments)
    } else {
      # Concatenate (mc+1) consecutive fragments to simulate missed cleavages
      for (i in seq_len(n - mc)) {
        peptides <- c(peptides, paste(fragments[i:(i + mc)], collapse = ""))
      }
    }
  }

  len <- nchar(peptides)
  n_observable <- sum(len >= min_len & len <= max_len, na.rm = TRUE)
  n_observable
}
