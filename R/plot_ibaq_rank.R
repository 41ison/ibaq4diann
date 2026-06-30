#' Plot ranked protein abundances per sample
#'
#' @description
#' Visualises protein abundances as a **rank-abundance (Whittaker) plot**:
#' proteins are ranked from most to least abundant along the x-axis, and their
#' iBAQ values are shown on the y-axis. One panel is drawn per sample when
#' multiple samples are selected (`facet_wrap`).
#'
#' This plot is useful for assessing the dynamic range of quantification and
#' for identifying the most abundant proteins in each sample.
#'
#' @param ibaq A `data.frame` or `tibble` returned by [compute_ibaq()].
#'   Accepted in both **wide** format (default output, one column per sample)
#'   and **long** format (`output_long = TRUE`, columns `sample` and `iBAQ`).
#' @param samples `character` vector or `NULL`. Names of samples (columns in
#'   wide format, or values in the `sample` column of long format) to include.
#'   When `NULL` (default), all samples are plotted.
#' @param log2_values `logical(1)`. Are the iBAQ values already
#'   log2-transformed (i.e., `log2_transform = TRUE` was used in
#'   [compute_ibaq()])? This only affects the y-axis label. Default: `TRUE`.
#' @param top_n_labels `integer(1)`. Number of top-ranked proteins (by iBAQ
#'   value, per sample) to annotate with text labels. Set to `0` to suppress
#'   all labels. Default: `10`.
#' @param label_col `character(1)`. Name of the column in `ibaq` to use as
#'   the protein label text. Falls back to `"protein_id"` if the specified
#'   column is absent. Default: `"Protein.Names"`.
#' @param point_size `numeric(1)`. Size of the scatter points. Default: `0.8`.
#' @param point_alpha `numeric(1)`. Opacity (0–1) of the scatter points.
#'   Default: `0.6`.
#' @param color_palette `character` or `NULL`. A named character vector mapping
#'   sample names to colours, or `NULL` to use the default ggplot2 colour
#'   scale. Default: `NULL`.
#'
#' @return A [`ggplot2`][ggplot2::ggplot] object. Use [ggplot2::ggsave()] to
#'   save to disk, or display it directly in an interactive R session.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ibaq <- compute_ibaq("report.parquet", "proteome.fasta")
#'
#' # ── All samples, default labels ─────────────────────────────────────────────
#' p <- plot_ibaq_rank(ibaq)
#' print(p)
#'
#' # ── Subset to two samples, label top 5 proteins ────────────────────────────
#' p2 <- plot_ibaq_rank(
#'   ibaq,
#'   samples      = c("sample_A", "sample_B"),
#'   top_n_labels = 5,
#'   label_col    = "Genes"
#' )
#' ggplot2::ggsave("rank_abundance.pdf", p2, width = 10, height = 5)
#' }
plot_ibaq_rank <- function(
  ibaq,
  samples      = NULL,
  log2_values  = TRUE,
  top_n_labels = 10,
  label_col    = "Protein.Names",
  point_size   = 0.8,
  point_alpha  = 0.6,
  color_palette = NULL
) {
  # ---- Determine metadata columns -------------------------------------------
  meta_cols <- c("protein_id", "Protein.Names", "Genes", "n_theoretical_peptides")
  meta_present <- intersect(meta_cols, names(ibaq))

  # ---- Convert wide → long if needed ----------------------------------------
  if (!"sample" %in% names(ibaq)) {
    sample_cols <- setdiff(names(ibaq), meta_present)
    if (!is.null(samples)) {
      sample_cols <- intersect(sample_cols, samples)
    }
    ibaq_long <- tidyr::pivot_longer(
      ibaq,
      cols      = dplyr::all_of(sample_cols),
      names_to  = "sample",
      values_to = "iBAQ"
    )
  } else {
    ibaq_long <- ibaq
    if (!is.null(samples)) {
      ibaq_long <- dplyr::filter(ibaq_long, .data$sample %in% samples)
    }
  }

  # ---- Resolve label column --------------------------------------------------
  if (!label_col %in% names(ibaq_long)) {
    label_col <- "protein_id"
  }

  # ---- Rank proteins within each sample (highest iBAQ = rank 1) -------------
  ibaq_long <- ibaq_long |>
    dplyr::filter(!is.na(.data$iBAQ)) |>
    dplyr::group_by(.data$sample) |>
    dplyr::arrange(dplyr::desc(.data$iBAQ), .by_group = TRUE) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()

  # ---- Identify top-N proteins for labelling ---------------------------------
  label_data <- if (top_n_labels > 0) {
    dplyr::filter(ibaq_long, .data$rank <= top_n_labels)
  } else {
    ibaq_long[0, ]
  }

  # ---- y-axis label ----------------------------------------------------------
  y_label <- if (isTRUE(log2_values)) "log\u2082(iBAQ + 1)" else "iBAQ"

  # ---- Build plot ------------------------------------------------------------
  p <- ggplot2::ggplot(
    ibaq_long,
    ggplot2::aes(
      x     = .data$rank,
      y     = .data$iBAQ,
      colour = .data$sample
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::facet_wrap(~ sample, scales = "free_x") +
    ggplot2::labs(
      title    = "Ranked Protein Abundances",
      subtitle = "Proteins ranked from most to least abundant (iBAQ)",
      x        = "Protein rank",
      y        = y_label,
      colour   = "Sample"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "#2d3561", colour = NA),
      strip.text        = ggplot2::element_text(colour = "white", face = "bold"),
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "none"
    )

  if (!is.null(color_palette)) {
    p <- p + ggplot2::scale_colour_manual(values = color_palette)
  }

  if (top_n_labels > 0 && nrow(label_data) > 0) {
    p <- p + ggplot2::geom_text(
      data    = label_data,
      mapping = ggplot2::aes(label = .data[[label_col]]),
      size    = 2.5,
      hjust   = -0.15,
      colour  = "grey20",
      show.legend = FALSE
    )
  }

  p
}
