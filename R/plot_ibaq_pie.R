#' Plot proportional iBAQ abundance as a pie chart
#'
#' @description
#' Draws a **pie chart** of proportional iBAQ (riBAQ) values per sample.
#' Each protein's contribution is expressed as a percentage of the total
#' sample abundance:
#'
#' \deqn{\text{riBAQ}_i = \frac{\text{iBAQ}_i}{\sum_j \text{iBAQ}_j} \times 100}
#'
#' Only the `top_n` most-abundant proteins are shown as individual slices;
#' the remaining proteins are grouped into an "Other" category. Optionally,
#' protein names (or gene symbols) can be displayed on the slices using
#' [ggrepel::geom_text_repel()].
#'
#' @param ibaq A `data.frame` or `tibble` returned by [compute_ibaq()], in
#'   either wide or long format.
#' @param samples `character` vector or `NULL`. Samples to include. When
#'   `NULL` (default), all samples are plotted.
#' @param log2_values `logical(1)`. Are the iBAQ values already
#'   log2-transformed (i.e., `log2_transform = TRUE` was used in
#'   [compute_ibaq()])? When `TRUE` (default), the function back-transforms
#'   via `2^x - 1` before computing proportions.
#' @param top_n `integer(1)`. Number of top proteins (by proportional
#'   abundance) to display as individual slices. All remaining proteins are
#'   collapsed into an "Other" slice. Default: `10`.
#' @param label_col `character(1)`. Column in `ibaq` to use for protein
#'   labels. Falls back to `"protein_id"` if the column is absent.
#'   Default: `"Protein.Names"`.
#' @param show_labels `logical(1)`. Show percentage labels on slices?
#'   Default: `TRUE`.
#' @param show_names `logical(1)`. Show protein name labels for the top-N
#'   slices using [ggrepel::geom_text_repel()]? Default: `TRUE`.
#' @param color_palette `character` or `NULL`. A character vector of colours
#'   (length >= `top_n + 1`), or `NULL` to use a default pastel palette.
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
#' # Default pie chart for all samples
#' plot_ibaq_pie(ibaq)
#'
#' # Show top 5 proteins, label by gene name
#' plot_ibaq_pie(
#'   ibaq,
#'   top_n = 5,
#'   label_col = "Genes",
#'   samples = c("sample_A", "sample_B")
#' )
#' }
plot_ibaq_pie <- function(
  ibaq,
  samples = NULL,
  log2_values = TRUE,
  top_n = 10,
  label_col = "Protein.Names",
  show_labels = TRUE,
  show_names = TRUE,
  color_palette = NULL
) {
  # ---- Determine metadata columns -------------------------------------------
  meta_cols <- c(
    "protein_id",
    "Protein.Names",
    "Genes",
    "n_theoretical_peptides"
  )
  meta_present <- dplyr::intersect(meta_cols, names(ibaq))

  # ---- Convert wide -> long if needed ----------------------------------------
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

  # ---- Resolve label column --------------------------------------------------
  if (!label_col %in% names(ibaq_long)) {
    label_col <- "protein_id"
  }

  # ---- Remove NAs ------------------------------------------------------------
  ibaq_long <- dplyr::filter(ibaq_long, !is.na(.data$iBAQ))

  # ---- Back-transform from log2 if needed ------------------------------------
  if (isTRUE(log2_values)) {
    ibaq_long <- dplyr::mutate(
      ibaq_long,
      iBAQ_linear = 2^.data$iBAQ - 1
    )
  } else {
    ibaq_long <- dplyr::mutate(ibaq_long, iBAQ_linear = .data$iBAQ)
  }

  # Ensure no negative values after back-transformation
  ibaq_long <- dplyr::mutate(
    ibaq_long,
    iBAQ_linear = pmax(.data$iBAQ_linear, 0)
  )

  # ---- Compute proportional iBAQ (riBAQ) per sample -------------------------
  ibaq_long <- ibaq_long |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(
      total_iBAQ = sum(.data$iBAQ_linear, na.rm = TRUE),
      pct = .data$iBAQ_linear / .data$total_iBAQ * 100
    ) |>
    dplyr::ungroup()

  # ---- Identify top-N proteins per sample ------------------------------------
  top_proteins <- ibaq_long |>
    dplyr::group_by(.data$sample) |>
    dplyr::arrange(dplyr::desc(.data$pct), .by_group = TRUE) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()

  # Assign label: top-N keep their name, others become "Other"
  top_proteins <- top_proteins |>
    dplyr::mutate(
      slice_label = dplyr::if_else(
        .data$rank <= top_n,
        as.character(.data[[label_col]]),
        "Other"
      )
    )

  # ---- Aggregate "Other" proteins per sample ---------------------------------
  pie_data <- top_proteins |>
    dplyr::group_by(.data$sample, .data$slice_label) |>
    dplyr::summarise(
      pct = sum(.data$pct, na.rm = TRUE),
      .groups = "drop"
    )

  # ---- Order slices: top-N by descending pct, "Other" last -------------------
  pie_data <- pie_data |>
    dplyr::group_by(.data$sample) |>
    dplyr::arrange(
      .data$slice_label == "Other",
      dplyr::desc(.data$pct),
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      slice_label = factor(
        .data$slice_label,
        levels = unique(.data$slice_label)
      )
    ) |>
    dplyr::ungroup()

  # ---- Compute label positions (midpoint of cumulative pct) ------------------
  pie_data <- pie_data |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(
      cum_pct = cumsum(.data$pct),
      y_pos = .data$cum_pct - .data$pct / 2
    ) |>
    dplyr::ungroup()

  # ---- Default colour palette ------------------------------------------------
  n_slices <- length(unique(pie_data$slice_label))
  if (is.null(color_palette)) {
    color_palette <- grDevices::hcl.colors(n_slices, palette = "Set 2")
  }

  # ---- Build plot ------------------------------------------------------------
  p <- ggplot2::ggplot(
    pie_data,
    ggplot2::aes(
      x = "",
      y = .data$pct,
      fill = .data$slice_label
    )
  ) +
    ggplot2::geom_col(width = 1, colour = "white", linewidth = 0.3) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~sample) +
    ggplot2::scale_fill_manual(
      values = color_palette,
      name = "Protein"
    ) +
    ggplot2::labs(
      title = "Proportional iBAQ Abundance",
      subtitle = paste0(
        "Top ", top_n,
        " proteins shown individually; remaining grouped as \u201cOther\u201d"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40")
    )

  # ---- Percentage labels on slices -------------------------------------------
  if (isTRUE(show_labels)) {
    # Only label slices >= 1 % so tiny slivers stay clean
    label_df <- dplyr::filter(pie_data, .data$pct >= 1)
    p <- p +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(
          y = .data$y_pos,
          label = paste0(round(.data$pct, 1), "%")
        ),
        size = 2.5,
        colour = "grey20",
        show.legend = FALSE
      )
  }

  # ---- Protein name labels via ggrepel ---------------------------------------
  if (isTRUE(show_names)) {
    rlang::check_installed(
      "ggrepel",
      reason = "to label protein slices on the pie chart."
    )
    name_df <- dplyr::filter(
      pie_data,
      .data$slice_label != "Other",
      .data$pct >= 1
    )
    p <- p +
      ggrepel::geom_text_repel(
        data = name_df,
        ggplot2::aes(y = .data$y_pos, label = .data$slice_label),
        size = 2.2,
        colour = "grey30",
        nudge_x = 0.6,
        show.legend = FALSE,
        max.overlaps = Inf,
        segment.colour = "grey60",
        segment.size = 0.2,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0
      )
  }

  p
}
