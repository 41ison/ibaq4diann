#' Parse a FASTA file into a named character vector
#'
#' @description
#' Reads a (possibly gzipped) FASTA file and returns a named character vector
#' of amino-acid sequences. The names are protein identifiers extracted from
#' each header line; which identifier is used is controlled by `id_type`, so
#' the FASTA can be keyed by the same type of identifier as the DIA-NN report
#' column chosen in [compute_ibaq()].
#'
#' The parser understands UniProt `sp`/`tr` style headers
#' (`>db|ACCESSION|ENTRY_NAME description OS=... GN=GENE ...`) and plain
#' headers (`>ACCESSION description ...`).
#'
#' @param fasta_path `character(1)`. Path to the FASTA file. May be plain text
#'   or gzip-compressed (`.gz`). The file must exist and be readable.
#' @param id_type `character(1)`. Which identifier to use as the sequence
#'   name. One of:
#'   - `"accession"` (default): UniProt accession (`ACCESSION` in
#'     `>db|ACCESSION|ENTRY_NAME`), or the first whitespace-delimited token of
#'     a plain header. Matches the `Protein.Group` / `Protein.Ids` columns of a
#'     DIA-NN report.
#'   - `"entry_name"`: UniProt entry name (`ENTRY_NAME` in
#'     `>db|ACCESSION|ENTRY_NAME`, e.g. `ALBU_HUMAN`). Matches the
#'     `Protein.Names` column.
#'   - `"gene"`: gene symbol taken from the `GN=` field of the header.
#'     Matches the `Genes` column. Entries without a `GN=` field are dropped.
#' @param id_pattern `character(1)` or `NULL`. Only used when
#'   `id_type = "accession"`. A regular expression with **two** capture groups
#'   used to extract the accession from each header line:
#'   - Group 1 should capture the accession from UniProt-style headers
#'     (`>db|ACCESSION|name`).
#'   - Group 2 should capture the accession from plain headers
#'     (`>ACCESSION description`).
#'   When `NULL` (default), a built-in pattern that handles both formats is
#'   used automatically.
#'
#' @details
#' Identifiers are not guaranteed to be unique in a FASTA file. In particular,
#' several entries (isoforms, unreviewed duplicates) may share the same gene
#' symbol. When duplicates occur, only the **first** entry in file order is
#' kept and a message reports how many were discarded. UniProt proteome
#' downloads list reviewed (`sp`) entries before unreviewed (`tr`) ones, so
#' the reviewed canonical entry is usually the one retained.
#'
#' @return A named `character` vector where:
#'   - **names** are protein identifiers of the requested `id_type`.
#'   - **values** are the corresponding amino-acid sequences as plain strings.
#'
#' @seealso [build_theoretical_peptide_counts()], [compute_ibaq()]
#' @keywords internal
parse_fasta <- function(
  fasta_path,
  id_type = c("accession", "entry_name", "gene"),
  id_pattern = NULL
) {
  id_type <- match.arg(id_type)

  message("[parse_fasta] Reading: ", fasta_path)
  seqs <- Biostrings::readAAStringSet(fasta_path)
  headers <- names(seqs) # Biostrings strips the leading '>'

  ids <- switch(
    id_type,
    accession = extract_fasta_accession(headers, id_pattern),
    entry_name = extract_fasta_entry_name(headers),
    gene = extract_fasta_gene(headers)
  )

  sequences <- as.character(seqs)
  names(sequences) <- ids

  # Drop entries for which no identifier could be extracted (e.g. no GN= field)
  no_id <- is.na(ids) | ids == ""
  if (any(no_id)) {
    message(
      "[parse_fasta] ",
      sum(no_id),
      " entr",
      if (sum(no_id) == 1) "y has" else "ies have",
      " no '",
      id_type,
      "' identifier and ",
      if (sum(no_id) == 1) "was" else "were",
      " dropped."
    )
    sequences <- sequences[!no_id]
  }

  # Keep the first occurrence of duplicated identifiers
  dup <- duplicated(names(sequences))
  if (any(dup)) {
    message(
      "[parse_fasta] ",
      sum(dup),
      " entr",
      if (sum(dup) == 1) "y shares" else "ies share",
      " a '",
      id_type,
      "' identifier with an earlier entry; keeping the first occurrence."
    )
    sequences <- sequences[!dup]
  }

  message(
    "[parse_fasta] Parsed ",
    length(sequences),
    " protein sequences keyed by '",
    id_type,
    "'."
  )
  sequences
}

# ---- Header parsers ----------------------------------------------------------

#' Extract UniProt accessions (or first token) from FASTA headers
#' @param headers `character`. FASTA header lines without the leading `>`.
#' @param id_pattern `character(1)` or `NULL`. See [parse_fasta()].
#' @return `character` vector of accessions.
#' @keywords internal
#' @noRd
extract_fasta_accession <- function(headers, id_pattern = NULL) {
  # Default pattern handles:
  #   >db|ACCESSION|name ...  (UniProt sp/tr)  -> group 1: ACCESSION
  #   >ACCESSION description ...               -> group 2: ACCESSION
  if (is.null(id_pattern)) {
    id_pattern <- "^[a-z]{2}[|]([^|]+)[|]|^([^ ]+)"
  }

  m <- stringr::str_match(headers, id_pattern)
  ids <- ifelse(!is.na(m[, 2]), m[, 2], m[, 3])

  # Safety fallback for headers that don't match either group
  still_missing <- is.na(ids)
  if (any(still_missing)) {
    ids[still_missing] <- stringr::str_extract(
      headers[still_missing],
      "^\\S+"
    )
    warning(
      sum(still_missing),
      " header(s) could not be parsed; used first whitespace-delimited token."
    )
  }
  ids
}

#' Extract UniProt entry names (e.g. `ALBU_HUMAN`) from FASTA headers
#' @inheritParams extract_fasta_accession
#' @return `character` vector of entry names; the first whitespace-delimited
#'   token for non-UniProt headers.
#' @keywords internal
#' @noRd
extract_fasta_entry_name <- function(headers) {
  ids <- stringr::str_match(headers, "^[a-z]{2}[|][^|]+[|](\\S+)")[, 2]
  plain <- is.na(ids)
  if (any(plain)) {
    ids[plain] <- stringr::str_extract(headers[plain], "^\\S+")
  }
  ids
}

#' Extract gene symbols from the `GN=` field of FASTA headers
#' @inheritParams extract_fasta_accession
#' @return `character` vector of gene symbols; `NA` where no `GN=` field exists.
#' @keywords internal
#' @noRd
extract_fasta_gene <- function(headers) {
  stringr::str_match(headers, "(?:^|\\s)GN=(\\S+)")[, 2]
}
