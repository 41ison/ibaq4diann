#' Parse a FASTA file into a named character vector
#'
#' @description
#' Reads a (possibly gzipped) FASTA file and returns a named character vector
#' of amino-acid sequences. Protein identifiers are extracted from the FASTA
#' header line using a flexible regular expression that handles both UniProt
#' `sp`/`tr` style headers (`>db|ACCESSION|name ...`) and plain headers
#' (`>ACCESSION description ...`).
#'
#' @param fasta_path `character(1)`. Path to the FASTA file. May be plain text
#'   or gzip-compressed (`.gz`). The file must exist and be readable.
#' @param id_pattern `character(1)` or `NULL`. A regular expression with **two**
#'   capture groups used to extract the protein accession from each header line.
#'   - Group 1 should capture the accession from UniProt-style headers
#'     (`>db|ACCESSION|name`).
#'   - Group 2 should capture the accession from plain headers
#'     (`>ACCESSION description`).
#'   When `NULL` (default), a built-in pattern that handles both formats is
#'   used automatically.
#'
#' @return A named `character` vector where:
#'   - **names** are protein accession identifiers parsed from the FASTA headers.
#'   - **values** are the corresponding amino-acid sequences as plain strings.
#'
#' @seealso [build_theoretical_peptide_counts()], [compute_ibaq()]
#' @keywords internal
parse_fasta <- function(
  fasta_path,
  id_pattern = NULL
) {
  # Default pattern handles:
  #   >db|ACCESSION|name ...  (UniProt sp/tr)  -> captures group 1: ACCESSION
  #   >ACCESSION description ...               -> captures group 2: ACCESSION
  if (is.null(id_pattern)) {
    id_pattern <- "^[a-z]{2}[|]([^|]+)[|]|^([^ ]+)"
  }

  message("[parse_fasta] Reading: ", fasta_path)
  seqs <- Biostrings::readAAStringSet(fasta_path)
  headers <- names(seqs) # Biostrings strips the leading '>'

  # Match two-group regex: group 2 = UniProt accession, group 3 = plain header
  m <- stringr::str_match(headers, id_pattern)
  ids <- ifelse(!is.na(m[, 2]), m[, 2], m[, 3])

  # Safety fallback for headers that don't match either group
  still_missing <- is.na(ids)
  if (any(still_missing)) {
    ids[still_missing] <- stringr::str_extract(
      headers[still_missing], "^\\S+"
    )
    warning(
      sum(still_missing),
      " header(s) could not be parsed; used first whitespace-delimited token."
    )
  }

  sequences <- as.character(seqs)
  names(sequences) <- ids
  message("[parse_fasta] Parsed ", length(sequences), " protein sequences.")
  sequences
}
