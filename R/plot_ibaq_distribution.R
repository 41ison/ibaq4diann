#' Plot iBAQ value distributions across samples
#'
#' @description
#' Draws a violin plot, box plot, or both overlaid for the distribution of
#' iBAQ values across samples. Useful for quality-control checks such as
#' assessing between-sample normalisation and detecting outlier runs.
#'
#' @param ibaq A `data.frame` or `tibble` returned by [compute_ibaq()], in
#'   either wide or long format.
#' @param samples `character` vector or `NULL`. Samples to include. When
#'   `NULL` (default), all samples are shown.
#' @param plot_type `character(1)`. One of `"violin"` (default), `"boxplot"`,
#'   or `"both"` (violin with overlaid box plot).
#' @param log2_values `logical(1)`. Are values already log2-transformed?
#'   Affects only the y-axis label. Default: `TRUE`.
#' @param fill_alpha `numeric(1)`. Fill opacity for violin/box shapes (0–1).
#'   Default: `0.6`.
#' @param color_palette `character` or `NULL`. Named vector mapping sample
#'   names to fill colours, or `NULL` to use the default ggplot2 palette.
#'   Default: `NULL`.
#'
#' @return A [`ggplot2`][ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ibaq <- compute_ibaq("report.parquet", "proteome.fasta")
#'
#' # Violin plot for all samples
#' plot_ibaq_distribution(ibaq)
#'
#' # Box plot for selected samples
#' plot_ibaq_distribution(
#'   ibaq,
#'   samples    = c("sample_A", "sample_B", "sample_C"),
#'   plot_type  = "boxplot"
#' )
#'
#' # Both violin + box
#' p <- plot_ibaq_distribution(ibaq, plot_type = "both")
#' ggplot2::ggsave("ibaq_distributions.pdf", p, width = 10, height = 6)
#' }
plot_ibaq_distribution <- function(
    ibaq,
    samples = NULL,
    plot_type = c("violin", "boxplot", "both"),
    log2_values = TRUE,
    fill_alpha = 0.6,
    color_palette = NULL) {
  plot_type <- match.arg(plot_type)

  # ---- Wide → long conversion ------------------------------------------------
  meta_cols <- c("protein_id", "Protein.Names", "Genes", "n_theoretical_peptides")
  meta_present <- intersect(meta_cols, names(ibaq))

  if (!"sample" %in% names(ibaq)) {
    sample_cols <- setdiff(names(ibaq), meta_present)
    if (!is.null(samples)) sample_cols <- intersect(sample_cols, samples)
    ibaq_long <- tidyr::pivot_longer(
      ibaq,
      cols = dplyr::all_of(sample_cols),
      names_to = "sample",
      values_to = "iBAQ"
    )
  } else {
    ibaq_long <- ibaq
    if (!is.null(samples)) {
      ibaq_long <- dplyr::filter(ibaq_long, .data$sample %in% samples)
    }
  }

  ibaq_long <- dplyr::filter(ibaq_long, !is.na(.data$iBAQ))

  y_label <- if (isTRUE(log2_values)) "log\u2082(iBAQ + 1)" else "iBAQ"

  p <- ggplot2::ggplot(
    ibaq_long,
    ggplot2::aes(x = .data$sample, y = .data$iBAQ, fill = .data$sample)
  ) +
    ggplot2::labs(
      title    = "iBAQ Value Distributions per Sample",
      subtitle = "Violin/box plots of protein abundance distributions",
      x        = "Sample",
      y        = y_label,
      fill     = "Sample"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none"
    )

  if (plot_type %in% c("violin", "both")) {
    p <- p + ggplot2::geom_violin(alpha = fill_alpha, trim = FALSE)
  }
  if (plot_type %in% c("boxplot", "both")) {
    box_width <- if (plot_type == "both") 0.15 else 0.5
    p <- p + ggplot2::geom_boxplot(
      alpha         = if (plot_type == "both") 0.9 else fill_alpha,
      width         = box_width,
      outlier.size  = 0.5,
      outlier.alpha = 0.4
    )
  }

  if (!is.null(color_palette)) {
    p <- p + ggplot2::scale_fill_manual(values = color_palette)
  }

  p
}
