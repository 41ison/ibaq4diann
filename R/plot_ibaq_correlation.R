#' Plot a sample-to-sample correlation heatmap of iBAQ values
#'
#' @description
#' Computes pairwise Pearson (or Spearman) correlations of log-transformed iBAQ
#' values across all selected samples and displays the result as a colour-coded
#' heatmap. This is a standard quality-control visualisation for assessing
#' reproducibility between biological or technical replicates.
#'
#' @param ibaq A `data.frame` or `tibble` returned by [compute_ibaq()], in
#'   either wide or long format.
#' @param samples `character` vector or `NULL`. Samples to include. When `NULL`
#'   (default), all samples are used.
#' @param method `character(1)`. Correlation method: `"pearson"` (default) or
#'   `"spearman"`. Passed to [stats::cor()].
#' @param log2_values `logical(1)`. Are values already log2-transformed?
#'   When `FALSE`, a `log2(x + 1)` transformation is applied before computing
#'   correlations. Default: `TRUE`.
#' @param low_color `character(1)`. Colour for low correlation values.
#'   Default: `"#f7fbff"` (near-white blue).
#' @param high_color `character(1)`. Colour for high correlation values.
#'   Default: `"#08306b"` (dark blue).
#' @param show_values `logical(1)`. If `TRUE`, print the correlation
#'   coefficient inside each tile. Default: `TRUE`.
#' @param digits `integer(1)`. Number of decimal places for displayed
#'   correlation values. Default: `2`.
#'
#' @return A [`ggplot2`][ggplot2::ggplot] object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ibaq <- compute_ibaq("report.parquet", "proteome.fasta")
#'
#' # Default Pearson correlation heatmap
#' plot_ibaq_correlation(ibaq)
#'
#' # Spearman correlation, no value labels
#' p <- plot_ibaq_correlation(
#'   ibaq,
#'   method      = "spearman",
#'   show_values = FALSE
#' )
#' ggplot2::ggsave("correlation_heatmap.pdf", p, width = 8, height = 7)
#' }
plot_ibaq_correlation <- function(
  ibaq,
  samples      = NULL,
  method       = c("pearson", "spearman"),
  log2_values  = TRUE,
  low_color    = "#f7fbff",
  high_color   = "#08306b",
  show_values  = TRUE,
  digits       = 2
) {
  method <- match.arg(method)

  # ---- Wide → long conversion ------------------------------------------------
  meta_cols    <- c("protein_id", "Protein.Names", "Genes", "n_theoretical_peptides")
  meta_present <- intersect(meta_cols, names(ibaq))

  if (!"sample" %in% names(ibaq)) {
    sample_cols <- setdiff(names(ibaq), meta_present)
    if (!is.null(samples)) sample_cols <- intersect(sample_cols, samples)
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

  # ---- Build wide matrix for correlation (proteins × samples) ---------------
  ibaq_wide <- ibaq_long |>
    dplyr::select("protein_id", "sample", "iBAQ") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "iBAQ")

  mat <- as.matrix(ibaq_wide[, -1])
  rownames(mat) <- ibaq_wide$protein_id

  if (!isTRUE(log2_values)) {
    mat <- log2(mat + 1)
  }

  # ---- Compute pairwise correlations ----------------------------------------
  cor_mat <- stats::cor(mat, use = "pairwise.complete.obs", method = method)

  # ---- Tidy the correlation matrix for ggplot2 ------------------------------
  cor_df <- as.data.frame(cor_mat)
  cor_df$sample_x <- rownames(cor_df)
  cor_long <- tidyr::pivot_longer(
    cor_df,
    cols      = -"sample_x",
    names_to  = "sample_y",
    values_to = "correlation"
  )

  # Fix factor order so the diagonal runs top-left → bottom-right
  samp_order <- rownames(cor_mat)
  cor_long$sample_x <- factor(cor_long$sample_x, levels = samp_order)
  cor_long$sample_y <- factor(cor_long$sample_y, levels = rev(samp_order))

  p <- ggplot2::ggplot(
    cor_long,
    ggplot2::aes(x = .data$sample_x, y = .data$sample_y, fill = .data$correlation)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(
      low    = low_color,
      high   = high_color,
      limits = c(NA, 1),
      name   = paste0(tools::toTitleCase(method), "\ncorrelation")
    ) +
    ggplot2::labs(
      title    = "Sample-to-Sample Correlation Heatmap",
      subtitle = paste0(
        tools::toTitleCase(method),
        " correlation of iBAQ values across samples"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid        = ggplot2::element_blank()
    )

  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = round(.data$correlation, digits)),
      size = 3,
      colour = "grey20"
    )
  }

  p
}
