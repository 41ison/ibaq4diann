#' Plot the dynamic range and protein detection coverage across samples
#'
#' @description
#' Displays the **dynamic range** of iBAQ quantification per sample: the
#' minimum, median, and maximum iBAQ value are annotated, and the fraction of
#' quantified proteins is indicated. This is a useful QC visualisation for
#' assessing proteome coverage and signal dynamic range across runs.
#'
#' Each sample is represented as a vertical segment spanning its iBAQ value
#' range, with a point at the median. A rug on the y-axis shows the full
#' distribution.
#'
#' @param ibaq A `data.frame` or `tibble` returned by [compute_ibaq()], in
#'   either wide or long format.
#' @param samples `character` vector or `NULL`. Samples to include. When
#'   `NULL` (default), all samples are shown.
#' @param log2_values `logical(1)`. Are values already log2-transformed?
#'   Affects only the y-axis label. Default: `TRUE`.
#' @param show_n_proteins `logical(1)`. Annotate each sample bar with the
#'   number of quantified proteins. Default: `TRUE`.
#' @param point_size `numeric(1)`. Size of the median point. Default: `3`.
#' @param segment_linewidth `numeric(1)`. Width of the range segment line.
#'   Default: `1.2`.
#' @param color_palette `character` or `NULL`. Named vector mapping sample
#'   names to colours, or `NULL` to use the default ggplot2 palette.
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
#' # Default dynamic range plot
#' plot_ibaq_dynamic_range(ibaq)
#'
#' # Subset to specific samples, hide protein count labels
#' p <- plot_ibaq_dynamic_range(
#'   ibaq,
#'   samples = c("sample_A", "sample_B"),
#'   show_n_proteins = FALSE
#' )
#' ggplot2::ggsave("dynamic_range.pdf", p, width = 8, height = 6)
#' }
plot_ibaq_dynamic_range <- function(
  ibaq,
  samples = NULL,
  log2_values = TRUE,
  show_n_proteins = TRUE,
  point_size = 3,
  segment_linewidth = 1.2,
  color_palette = NULL
) {
  # ---- Wide → long conversion ------------------------------------------------
  meta_cols <- c(
    "protein_id",
    "Protein.Names",
    "Genes",
    "n_theoretical_peptides"
  )
  meta_present <- dplyr::intersect(meta_cols, names(ibaq))

  if (!"sample" %in% names(ibaq)) {
    sample_cols <- dplyr::setdiff(names(ibaq), meta_present)
    if (!is.null(samples)) {
      sample_cols <- dplyr::intersect(sample_cols, samples)
    }
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

  # ---- Compute per-sample summary stats --------------------------------------
  summary_df <- ibaq_long |>
    dplyr::group_by(.data$sample) |>
    dplyr::summarise(
      ymin = min(.data$iBAQ, na.rm = TRUE),
      ymax = max(.data$iBAQ, na.rm = TRUE),
      ymedian = stats::median(.data$iBAQ, na.rm = TRUE),
      n_quant = dplyr::n(),
      .groups = "drop"
    )

  y_label <- if (isTRUE(log2_values)) "log\u2082(iBAQ + 1)" else "iBAQ"
  y_label_range <- paste0(y_label, " range")

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = .data$sample, colour = .data$sample)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data$sample,
        y = .data$ymin,
        yend = .data$ymax
      ),
      linewidth = segment_linewidth
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$ymedian),
      size = point_size,
      shape = 18
    ) +
    ggplot2::labs(
      title = "Dynamic Range of iBAQ Quantification",
      subtitle = paste0(
        "Vertical segments show min\u2013max range; diamond = median. ",
        y_label
      ),
      x = "Sample",
      y = y_label_range,
      colour = "Sample"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )

  if (isTRUE(show_n_proteins)) {
    p <- p +
      ggplot2::geom_text(
        data = summary_df,
        mapping = ggplot2::aes(
          y = .data$ymax,
          label = paste0("n=", .data$n_quant)
        ),
        vjust = -0.5,
        size = 3,
        colour = "grey30",
        show.legend = FALSE
      )
  }

  if (!is.null(color_palette)) {
    p <- p + ggplot2::scale_colour_manual(values = color_palette)
  }

  p
}
