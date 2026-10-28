#' Compute iBAQ values from a DIA-NN report and a reference proteome FASTA
#'
#' @title Compute iBAQ Values from DIA-NN Output
#'
#' @description
#' Calculates **intensity-Based Absolute Quantification (iBAQ)** values from
#' DIA-NN search results (`report.parquet`) and a reference proteome FASTA
#' file.
#'
#' The iBAQ method normalises the summed precursor intensities of a protein by
#' the number of theoretically observable tryptic peptides, making protein
#' abundances comparable regardless of protein size:
#'
#' \deqn{iBAQ = \frac{\sum \text{precursor intensities}}{\text{# theoretical peptides}}}
#'
#' @details
#' ## Workflow
#' 1. **Load** the DIA-NN `report.parquet` using [arrow::read_parquet()].
#' 2. **Filter** low-quality precursors via [filter_diann_report()].
#' 3. **Aggregate** precursor intensities to the protein × sample level.
#'    Each unique stripped peptide sequence is counted once per sample (even
#'    if observed with multiple charge states or modifications), then summed
#'    per protein per sample.
#' 4. **Parse** the FASTA and count theoretically observable tryptic peptides
#'    per protein using [build_theoretical_peptide_counts()].
#' 5. **Divide** the summed intensity by the theoretical peptide count.
#' 6. **Optionally** log2-transform: `log2(iBAQ + 1)` (adds a pseudo-count of
#'    1 to avoid `log2(0)`).
#' 7. **Return** a wide- or long-format tibble.
#'
#' ## References
#' - Schwanhäusser B et al. (2011). Global quantification of mammalian gene
#'   expression control. *Nature* **473**:337–342.
#'   \doi{10.1038/nature10098}
#' - 41ison/iBAQ-for-DIANN (GitHub reference implementation)
#' - bigbio/ibaqpy (Python reference implementation)
#'
#' @param parquet_path `character(1)`. Path to the DIA-NN `report.parquet`
#'   file. The file must exist and be readable.
#' @param fasta_path `character(1)`. Path to the reference proteome FASTA file
#'   (plain text or gzip-compressed). Must contain the same protein sequences
#'   used during the DIA-NN database search.
#' @param intensity_col `character(1)`. Name of the intensity column in the
#'   report to use for quantification. Common choices:
#'   - `"Precursor.Normalised"` (default) — cross-run normalised intensities.
#'   - `"Precursor.Quantity"` — raw precursor intensities.
#' @param q_value_cutoff `numeric(1)`. Maximum precursor-level FDR (`Q.Value`).
#'   Default: `0.01` (1 %). Set to `NULL` to skip.
#' @param pg_q_value_cutoff `numeric(1)`. Maximum protein-group-level FDR
#'   (`PG.Q.Value`). Default: `0.01`. Set to `NULL` to skip.
#' @param lib_pg_q_value_cutoff `numeric(1)`. Maximum library protein-group
#'   FDR (`Lib.PG.Q.Value`). Default: `0.01`. Set to `NULL` to skip.
#' @param proteotypic_only `logical(1)`. Keep only proteotypic (unique-to-one-
#'   protein) peptides? Default: `TRUE`. Recommended for unambiguous
#'   quantification.
#' @param lfq_quality_cutoff `numeric(1)`. Minimum LFQ quality score
#'   (`PG.MaxLFQ.Quality`). Default: `0.5`. Set to `NULL` or `0` to skip.
#' @param protease_regex `character(1)`. Perl-compatible regex defining the
#'   in-silico cleavage rule. Default: `"(?<=[KR])(?!P)"` (trypsin, no cut
#'   immediately before Pro).
#' @param min_peptide_len `integer(1)`. Minimum peptide length (residues) for
#'   the theoretical peptide count. Default: `6`.
#' @param max_peptide_len `integer(1)`. Maximum peptide length (residues) for
#'   the theoretical peptide count. Default: `30`.
#' @param max_missed_cleavages `integer(1)`. Maximum number of missed cleavages
#'   allowed during in-silico digestion. Default: `0` (fully tryptic peptides
#'   only).
#' @param log2_transform `logical(1)`. Apply `log2(iBAQ + 1)` transformation
#'   to the final iBAQ values? Default: `TRUE`. Set to `FALSE` to return raw
#'   (linear-scale) iBAQ values.
#' @param output_long `logical(1)`. If `FALSE` (default), return a wide-format
#'   tibble with one column per sample. If `TRUE`, return a long-format tibble
#'   with columns `sample` and `iBAQ`.
#' @param id_col `character(1)`. Which DIA-NN identifier column defines a
#'   "protein" for quantification. One of:
#'   - `"Protein.Group"` (default): protein group, named by its leading
#'     accession (e.g. `"P04040"`). Matched to FASTA **accessions**.
#'   - `"Protein.Ids"`: all accessions in the group (e.g.
#'     `"P0DP23;P0DP24;P0DP25"`). Matched to FASTA accessions; the leading
#'     member found in the FASTA supplies the denominator.
#'   - `"Protein.Names"`: UniProt entry names (e.g. `"ALBU_HUMAN"`). Matched
#'     to FASTA **entry names**.
#'   - `"Genes"`: gene symbols (e.g. `"ALB"`). Matched to the `GN=` field of
#'     FASTA headers.
#'
#'   See **Choice of protein identifier** in Details — the number of output
#'   rows depends on this choice.
#' @param fasta_id_pattern `character(1)` or `NULL`. Custom regex (two capture
#'   groups) for extracting protein accessions from FASTA headers. Only used
#'   when `id_col` is `"Protein.Group"` or `"Protein.Ids"`. When `NULL`
#'   (default), a built-in pattern handles both UniProt and plain headers
#'   automatically. See [parse_fasta()] for details.
#'
#' @section Choice of protein identifier:
#' DIA-NN reports four identifier columns that describe the same protein
#' group at different levels of granularity, so the **number of proteins in
#' the output differs** depending on `id_col`:
#'
#' - `Protein.Ids` is the most granular: groups that share a leading protein
#'   but differ in their other members (e.g. isoform sets) are kept separate.
#' - `Protein.Group` and `Protein.Names` are one-to-one with each other.
#' - `Genes` is the coarsest: distinct protein groups from the same gene
#'   (e.g. isoforms with distinct accessions) are **merged** and their
#'   intensities summed. Rows with an empty `Genes` value are dropped with a
#'   message, since they cannot be matched to the FASTA.
#'
#' Identifiers may contain several members separated by `;` (including
#' contaminant-prefixed members such as `cRAP-P04040`). Members are looked up
#' in the FASTA in order and the first match — normally the leading protein —
#' supplies the theoretical peptide count. Identifiers with no member in the
#' FASTA are excluded with a warning.
#'
#' When `id_col = "Genes"`, several FASTA entries may share a gene symbol;
#' the first entry in file order is used (see [parse_fasta()]).
#'
#' @return
#' A [tibble][tibble::tibble]. The structure depends on `output_long`:
#'
#' **Wide format** (`output_long = FALSE`, default):
#' \describe{
#'   \item{`protein_id`}{Protein identifier (values of `id_col`).}
#'   \item{`Protein.Group`, `Protein.Names`, `Genes`}{The other identifier
#'     columns (whichever are present in the report and not used as
#'     `id_col`). When several report values map to one `protein_id` (e.g.
#'     `id_col = "Genes"`), the distinct values are joined with `;`.}
#'   \item{`n_theoretical_peptides`}{Number of theoretical tryptic peptides
#'     used as iBAQ denominator.}
#'   \item{`<sample_1>`, `<sample_2>`, ...}{One column per DIA-NN run
#'     (`Run` column), containing the iBAQ value for that protein in that
#'     sample. `NA` when the protein was not detected in a sample.}
#' }
#'
#' **Long format** (`output_long = TRUE`):
#' \describe{
#'   \item{`protein_id`}{Protein identifier.}
#'   \item{`Protein.Group`, `Protein.Names`, `Genes`}{Annotation columns, as
#'     above.}
#'   \item{`n_theoretical_peptides`}{Theoretical peptide count.}
#'   \item{`sample`}{Run/sample name.}
#'   \item{`iBAQ`}{iBAQ value for this protein–sample pair.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # ── Basic usage with default parameters ────────────────────────────────────
#' ibaq <- compute_ibaq(
#'   parquet_path = "report.parquet",
#'   fasta_path   = "UP000031575.fasta"
#' )
#'
#' head(ibaq)
#' dim(ibaq)
#'
#' # ── Custom intensity column and relaxed filters ─────────────────────────────
#' ibaq_raw <- compute_ibaq(
#'   parquet_path = "report.parquet",
#'   fasta_path = "UP000031575.fasta",
#'   intensity_col = "Precursor.Quantity",   # raw (non-normalised)
#'   q_value_cutoff = 0.01,
#'   pg_q_value_cutoff = 0.01,
#'   lib_pg_q_value_cutoff = 0.01,
#'   proteotypic_only = TRUE,
#'   lfq_quality_cutoff = NULL,              # skip LFQ quality filter
#'   min_peptide_len = 7,
#'   max_peptide_len = 25,
#'   max_missed_cleavages = 1,
#'   log2_transform = FALSE,                  # return linear iBAQ
#'   output_long = FALSE
#' )
#'
#' # ── Quantify at the gene level instead of the protein group ────────────────
#' ibaq_genes <- compute_ibaq(
#'   parquet_path = "report.parquet",
#'   fasta_path = "UP000031575.fasta",
#'   id_col = "Genes"
#' )
#'
#' # ── Long format output (tidy) ───────────────────────────────────────────────
#' ibaq_long <- compute_ibaq(
#'   parquet_path = "report.parquet",
#'   fasta_path = "UP000031575.fasta",
#'   output_long = TRUE
#' )
#'
#' # ── Save results to TSV ─────────────────────────────────────────────────────
#' readr::write_tsv(ibaq, "iBAQ_results.tsv")
#' }
compute_ibaq <- function(
  parquet_path,
  fasta_path,
  intensity_col = "Precursor.Normalised",
  q_value_cutoff = 0.01,
  pg_q_value_cutoff = 0.01,
  lib_pg_q_value_cutoff = 0.01,
  proteotypic_only = TRUE,
  lfq_quality_cutoff = 0.50,
  protease_regex = "(?<=[KR])(?!P)",
  min_peptide_len = 6,
  max_peptide_len = 30,
  max_missed_cleavages = 0,
  log2_transform = TRUE,
  output_long = FALSE,
  id_col = c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes"),
  fasta_id_pattern = NULL
) {
  stopifnot(
    "parquet_path must be a character string" = is.character(parquet_path),
    "fasta_path must be a character string" = is.character(fasta_path),
    file.exists(parquet_path),
    file.exists(fasta_path)
  )
  id_col <- match.arg(id_col)
  fasta_id_type <- diann_id_cols[[id_col]]

  ## ---- Step 1: Load parquet --------------------------------------------------
  message("\n=== Step 1: Loading DIA-NN report ===")
  report <- arrow::read_parquet(parquet_path)
  message(
    "[compute_ibaq] Loaded ",
    nrow(report),
    " rows, ",
    ncol(report),
    " columns."
  )

  # Validate required columns
  required_cols <- c(id_col, "Run", "Stripped.Sequence", intensity_col)
  missing_cols <- dplyr::setdiff(required_cols, names(report))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in report.parquet: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ## ---- Step 2: Quality filtering ---------------------------------------------
  message("\n=== Step 2: Quality filtering ===")
  report_filtered <- filter_diann_report(
    report,
    q_value_cutoff = q_value_cutoff,
    pg_q_value_cutoff = pg_q_value_cutoff,
    lib_pg_q_value_cutoff = lib_pg_q_value_cutoff,
    proteotypic_only = proteotypic_only,
    lfq_quality_cutoff = lfq_quality_cutoff
  )

  if (nrow(report_filtered) == 0) {
    stop("No rows passed quality filters. Please relax filter thresholds.")
  }

  ## ---- Step 3: Aggregate intensities to protein × sample --------------------
  message(
    "\n=== Step 3: Aggregating intensities per protein \u00d7 sample ",
    "(id_col = '",
    id_col,
    "') ==="
  )

  report_filtered <- dplyr::rename(
    report_filtered,
    protein_id = dplyr::all_of(id_col),
    intensity = dplyr::all_of(intensity_col)
  )

  # Rows without an identifier (e.g. empty Genes) cannot be matched to the
  # FASTA and are dropped here
  no_id <- is.na(report_filtered$protein_id) | report_filtered$protein_id == ""
  if (any(no_id)) {
    message(
      "[compute_ibaq] Dropping ",
      sum(no_id),
      " precursor rows with an empty '",
      id_col,
      "' value."
    )
    report_filtered <- report_filtered[!no_id, ]
  }

  # Annotation columns: the other identifier columns present in the report.
  # Several report values may map to a single protein_id (e.g. two protein
  # groups from one gene), so distinct values are collapsed with ";".
  annot_cols <- dplyr::intersect(
    dplyr::setdiff(c("Protein.Group", "Protein.Names", "Genes"), id_col),
    names(report_filtered)
  )
  annot_lookup <- report_filtered |>
    dplyr::distinct(.data$protein_id, dplyr::pick(dplyr::all_of(annot_cols))) |>
    dplyr::group_by(.data$protein_id) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(annot_cols),
        \(x) paste(unique(x[!is.na(x) & x != ""]), collapse = ";")
      ),
      .groups = "drop"
    )

  n_ids <- dplyr::n_distinct(report_filtered$protein_id)
  n_report_groups <- if ("Protein.Group" %in% names(report_filtered)) {
    dplyr::n_distinct(report_filtered$Protein.Group)
  } else {
    n_ids
  }
  if (id_col != "Protein.Group" && n_ids != n_report_groups) {
    message(
      "[compute_ibaq] ",
      n_report_groups,
      " protein groups in the filtered report correspond to ",
      n_ids,
      " distinct '",
      id_col,
      "' identifiers",
      if (n_ids < n_report_groups) {
        "; intensities of groups sharing an identifier are summed."
      } else {
        "."
      }
    )
  }

  # Sum intensity per unique stripped peptide per protein × sample, then
  # sum across unique peptides. Using Stripped.Sequence deduplication ensures
  # each peptide is counted once per sample even if observed across multiple
  # charge states or modifications.
  protein_intensity <- report_filtered |>
    dplyr::group_by(.data$Run, .data$protein_id, .data$Stripped.Sequence) |>
    dplyr::summarise(
      peptide_intensity = sum(.data$intensity, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$Run, .data$protein_id) |>
    dplyr::summarise(
      sum_intensity = sum(.data$peptide_intensity, na.rm = TRUE),
      .groups = "drop"
    )

  message(
    "[compute_ibaq] Aggregated ",
    dplyr::n_distinct(protein_intensity$protein_id),
    " proteins across ",
    dplyr::n_distinct(protein_intensity$Run),
    " samples."
  )

  ## ---- Step 4: Parse FASTA and count theoretical peptides -------------------
  message(
    "\n=== Step 4: Parsing FASTA & counting theoretical tryptic peptides ==="
  )
  sequences <- parse_fasta(
    fasta_path,
    id_type = fasta_id_type,
    id_pattern = fasta_id_pattern
  )
  theo_peptides <- build_theoretical_peptide_counts(
    sequences,
    protease_regex = protease_regex,
    min_len = min_peptide_len,
    max_len = max_peptide_len,
    max_missed = max_missed_cleavages
  )

  # Match each (possibly ';'-separated) identifier to the FASTA
  detected_proteins <- unique(protein_intensity$protein_id)
  matched <- match_theoretical_peptides(detected_proteins, theo_peptides)

  n_found <- sum(!is.na(matched$fasta_id))
  n_missing <- length(detected_proteins) - n_found
  n_via_member <- sum(
    !is.na(matched$fasta_id) & matched$fasta_id != matched$protein_id
  )

  if (n_found == 0) {
    stop(
      "None of the ",
      length(detected_proteins),
      " detected '",
      id_col,
      "' identifiers were found in the FASTA (matched by '",
      fasta_id_type,
      "'). Verify that the FASTA matches the search database."
    )
  }
  if (n_missing > 0) {
    warning(
      n_missing,
      " detected protein(s) not found in FASTA. ",
      "They will be excluded from the output. ",
      "Verify that the FASTA matches the search database."
    )
    message(
      "[compute_ibaq] Proteins not in FASTA (first 10): ",
      paste(
        utils::head(matched$protein_id[is.na(matched$fasta_id)], 10),
        collapse = ", "
      )
    )
  }
  message(
    "[compute_ibaq] ",
    n_found,
    " / ",
    length(detected_proteins),
    " detected proteins matched in FASTA",
    if (n_via_member > 0) {
      paste0(" (", n_via_member, " via a member of a multi-protein identifier)")
    },
    "."
  )

  ## ---- Step 5: Calculate iBAQ -----------------------------------------------
  message("\n=== Step 5: Computing iBAQ values ===")

  ibaq_long <- protein_intensity |>
    dplyr::left_join(
      dplyr::select(matched, "protein_id", "n_theoretical_peptides"),
      by = "protein_id"
    ) |>
    dplyr::left_join(annot_lookup, by = "protein_id") |>
    dplyr::filter(
      !is.na(.data$n_theoretical_peptides),
      .data$n_theoretical_peptides > 0 # guard against division by zero
    ) |>
    dplyr::mutate(
      iBAQ = .data$sum_intensity / .data$n_theoretical_peptides
    )

  if (isTRUE(log2_transform)) {
    ibaq_long <- dplyr::mutate(ibaq_long, iBAQ = log2(.data$iBAQ + 1))
    message("[compute_ibaq] iBAQ values have been log2(iBAQ + 1) transformed.")
  } else {
    message("[compute_ibaq] Returning raw (non-transformed) iBAQ values.")
  }

  ## ---- Step 6: Shape output --------------------------------------------------
  message("\n=== Step 6: Shaping output ===")

  annot_final_cols <- c("protein_id", annot_cols, "n_theoretical_peptides")

  if (isTRUE(output_long)) {
    result <- ibaq_long |>
      dplyr::select(dplyr::all_of(c(annot_final_cols, "Run", "iBAQ"))) |>
      dplyr::rename(sample = "Run") |>
      dplyr::arrange(.data$protein_id, .data$sample)
    message(
      "[compute_ibaq] Output: long format \u2014 ",
      nrow(result),
      " rows."
    )
  } else {
    result <- ibaq_long |>
      dplyr::select(dplyr::all_of(c(annot_final_cols, "Run", "iBAQ"))) |>
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(annot_final_cols),
        names_from = "Run",
        values_from = "iBAQ"
      ) |>
      dplyr::arrange(.data$protein_id)

    n_prot <- nrow(result)
    n_samples <- ncol(result) - length(annot_final_cols)
    message(
      "[compute_ibaq] Output: wide format \u2014 ",
      n_prot,
      " proteins \u00d7 ",
      n_samples,
      " samples."
    )
  }

  message("\n\u2705 iBAQ calculation complete!")
  result
}
