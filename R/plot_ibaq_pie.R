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
#' the remaining proteins are grouped into a neutral grey "Other" slice.
#' Slices are annotated with the protein name and/or percentage. Labels are
#' placed beside the pie with leader lines (via
#' [ggrepel::geom_text_repel()]) so that they remain legible even when many
#' slices are thin.
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
#' @param show_labels `logical(1)`. Include the percentage in slice labels?
#'   Default: `TRUE`.
#' @param show_names `logical(1)`. Include the protein name in slice labels?
#'   Default: `TRUE`.
#' @param min_label_pct `numeric(1)`. Slices smaller than this percentage are
#'   not labelled, to avoid clutter. Default: `1`.
#' @param show_legend `logical(1)`. Show a colour legend? Default: `FALSE`,
#'   as the on-plot labels already identify each slice. Useful when
#'   `show_names = FALSE`.
#' @param color_palette `character` or `NULL`. Colours for the protein
#'   slices: either a named vector (names are protein labels, optionally
#'   including `"Other"`) or an unnamed vector that is recycled over the
#'   proteins in order of first appearance. When `NULL` (default), a
#'   qualitative palette is used and "Other" is drawn in grey.
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
  min_label_pct = 1,
  show_legend = FALSE,
  color_palette = NULL
) {
  # ---- Convert wide -> long if needed ----------------------------------------
  ibaq_long <- ibaq_to_long(ibaq, samples)

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

  # ---- Compute proportional iBAQ (riBAQ) and top-N per sample ---------------
  pie_data <- ibaq_long |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(
      pct = .data$iBAQ_linear / sum(.data$iBAQ_linear, na.rm = TRUE) * 100
    ) |>
    dplyr::arrange(dplyr::desc(.data$pct), .by_group = TRUE) |>
    dplyr::mutate(
      slice_label = dplyr::if_else(
        dplyr::row_number() <= top_n,
        as.character(.data[[label_col]]),
        "Other"
      )
    ) |>
    dplyr::group_by(.data$sample, .data$slice_label) |>
    dplyr::summarise(pct = sum(.data$pct, na.rm = TRUE), .groups = "drop") |>
    # Order slices: top-N by descending pct, "Other" last
    dplyr::arrange(
      .data$sample,
      .data$slice_label == "Other",
      dplyr::desc(.data$pct)
    )

  # ---- Slice geometry (angles clockwise from 12 o'clock) --------------------
  pie_data <- pie_data |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(
      slice = dplyr::row_number(),
      theta_end = cumsum(.data$pct) / 100 * 2 * pi,
      theta_start = .data$theta_end - .data$pct / 100 * 2 * pi,
      theta_mid = (.data$theta_start + .data$theta_end) / 2
    ) |>
    dplyr::ungroup()

  # ---- Colour scale: proteins get a qualitative palette, "Other" is grey ----
  protein_levels <- unique(pie_data$slice_label[
    pie_data$slice_label != "Other"
  ])
  fill_levels <- c(protein_levels, "Other")

  if (is.null(color_palette)) {
    fill_values <- stats::setNames(
      grDevices::hcl.colors(length(protein_levels), palette = "Dark 3"),
      protein_levels
    )
  } else if (!is.null(names(color_palette))) {
    fill_values <- color_palette
  } else {
    fill_values <- stats::setNames(
      rep_len(color_palette, length(protein_levels)),
      protein_levels
    )
  }
  if (!"Other" %in% names(fill_values)) {
    fill_values[["Other"]] <- "grey85"
  }

  pie_data <- dplyr::mutate(
    pie_data,
    slice_label = factor(.data$slice_label, levels = fill_levels)
  )

  # ---- Polygon vertices for each slice ---------------------------------------
  polygons <- pie_data |>
    dplyr::rowwise() |>
    dplyr::reframe(
      sample = .data$sample,
      slice = .data$slice,
      slice_label = .data$slice_label,
      theta = c(
        .data$theta_start,
        seq(
          .data$theta_start,
          .data$theta_end,
          length.out = max(3L, ceiling(.data$pct * 2))
        ),
        .data$theta_end
      ),
      r = 1
    ) |>
    dplyr::mutate(
      x = .data$r * sin(.data$theta),
      y = .data$r * cos(.data$theta)
    )

  # Prepend the centre point to close each wedge
  centres <- dplyr::distinct(
    polygons,
    .data$sample,
    .data$slice,
    .data$slice_label
  ) |>
    dplyr::mutate(theta = NA_real_, r = 0, x = 0, y = 0)
  polygons <- dplyr::bind_rows(centres, polygons) |>
    dplyr::arrange(.data$sample, .data$slice, .data$r) |>
    dplyr::mutate(group = paste(.data$sample, .data$slice, sep = "_"))

  # ---- Build plot ------------------------------------------------------------
  p <- ggplot2::ggplot(polygons) +
    ggplot2::geom_polygon(
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$group,
        fill = .data$slice_label
      ),
      colour = "white",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = fill_levels,
      name = "Protein",
      drop = FALSE
    ) +
    ggplot2::coord_equal(
      xlim = c(-2.4, 2.4),
      ylim = c(-1.25, 1.25),
      expand = FALSE
    ) +
    ggplot2::facet_wrap(~sample) +
    ggplot2::labs(
      title = "Proportional iBAQ abundance (riBAQ)",
      subtitle = paste0(
        "Top ",
        top_n,
        " proteins per sample; remaining grouped as \u201cOther\u201d"
      )
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 11,
        margin = ggplot2::margin(b = 4)
      ),
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        colour = "grey40",
        size = 10
      ),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  # ---- Slice labels beside the pie with leader lines ------------------------
  if (isTRUE(show_labels) || isTRUE(show_names)) {
    rlang::check_installed("ggrepel", reason = "to label pie chart slices.")

    label_df <- pie_data |>
      dplyr::filter(.data$pct >= min_label_pct) |>
      dplyr::mutate(
        label = dplyr::case_when(
          isTRUE(show_names) & isTRUE(show_labels) ~
            paste0(
              .data$slice_label,
              " (",
              formatC(.data$pct, format = "f", digits = 1),
              "%)"
            ),
          isTRUE(show_names) ~ as.character(.data$slice_label),
          TRUE ~ paste0(formatC(.data$pct, format = "f", digits = 1), "%")
        ),
        # Leader line starts just inside the rim of the slice
        x = 0.95 * sin(.data$theta_mid),
        y = 0.95 * cos(.data$theta_mid),
        right_side = .data$x >= 0,
        # Labels sit in a column on the side of the pie their slice faces
        x_lab = dplyr::if_else(.data$right_side, 1.4, -1.4),
        y_lab = 1.1 * .data$y,
        hjust = dplyr::if_else(.data$right_side, 0, 1)
      )

    # A large "Other" slice is labelled inside the pie rather than beside it
    inside_df <- dplyr::filter(
      label_df,
      .data$slice_label == "Other",
      .data$pct >= 20
    )
    label_df <- dplyr::anti_join(
      label_df,
      inside_df,
      by = c("sample", "slice_label")
    )

    if (nrow(inside_df) > 0) {
      p <- p +
        ggplot2::geom_text(
          data = inside_df,
          ggplot2::aes(
            x = 0.55 * sin(.data$theta_mid),
            y = 0.55 * cos(.data$theta_mid),
            label = .data$label
          ),
          size = 2.8,
          colour = "grey20",
          show.legend = FALSE
        )
    }

    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        nudge_x = label_df$x_lab - label_df$x,
        nudge_y = label_df$y_lab - label_df$y,
        hjust = label_df$hjust,
        direction = "y",
        size = 2.6,
        colour = "grey20",
        segment.colour = "grey60",
        segment.size = 0.25,
        min.segment.length = 0,
        box.padding = 0.15,
        point.padding = 0,
        max.overlaps = Inf,
        seed = 1,
        show.legend = FALSE
      )
  }

  p
}
