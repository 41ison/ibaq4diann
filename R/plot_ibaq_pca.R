#' PCA plot of iBAQ values with KNN-based Louvain clustering
#'
#' @description
#' Performs **principal component analysis (PCA)** on the iBAQ protein
#' abundance matrix (samples as observations, proteins as features) and
#' displays the results as a 2D scatter plot. Samples are coloured by
#' cluster assignments derived from **K-nearest-neighbour (KNN) graph**
#' construction followed by **Louvain community detection**.
#'
#' This visualisation is useful for assessing sample-level similarity,
#' detecting batch effects, and identifying outlier runs.
#'
#' @param ibaq A `data.frame` or `tibble` returned by [compute_ibaq()], in
#'   either wide or long format.
#' @param samples `character` vector or `NULL`. Samples to include. When
#'   `NULL` (default), all samples are used.
#' @param log2_values `logical(1)`. Are values already log2-transformed?
#'   Affects only the subtitle annotation. Default: `TRUE`.
#' @param k `integer(1)`. Number of nearest neighbours for the KNN graph.
#'   Default: `5`. If `k` is greater than or equal to the number of samples,
#'   it is automatically reduced to `n_samples - 1`.
#' @param components `integer(2)`. Which principal components to plot.
#'   Default: `c(1, 2)` (PC1 vs PC2).
#' @param point_size `numeric(1)`. Size of the scatter points.
#'   Default: `3`.
#' @param label_samples `logical(1)`. If `TRUE`, annotate each point with
#'   the sample name using [ggrepel::geom_text_repel()]. Default: `TRUE`.
#' @param scale `logical(1)`. Whether to scale variables to unit variance
#'   before PCA. Default: `TRUE`.
#' @param color_palette `character` or `NULL`. A character vector of colours
#'   for cluster assignments, or `NULL` to use a default palette.
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
#' # Default PCA with KNN clustering
#' plot_ibaq_pca(ibaq)
#'
#' # Custom k and no sample labels
#' plot_ibaq_pca(ibaq, k = 3, label_samples = FALSE)
#'
#' # Plot PC2 vs PC3
#' plot_ibaq_pca(ibaq, components = c(2, 3))
#' }
plot_ibaq_pca <- function(
  ibaq,
  samples = NULL,
  log2_values = TRUE,
  k = 5,
  components = c(1, 2),
  point_size = 3,
  label_samples = TRUE,
  scale = TRUE,
  color_palette = NULL
) {
  rlang::check_installed("FNN", reason = "to build the KNN graph for clustering.")
  rlang::check_installed("igraph", reason = "to perform Louvain community detection.")

  # ---- Validate components ---------------------------------------------------
  stopifnot(
    "components must be an integer vector of length 2" =
      length(components) == 2 && all(components == as.integer(components)),
    "components must be positive" = all(components > 0)
  )

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

  # ---- Build samples x proteins matrix ---------------------------------------
  ibaq_wide <- ibaq_long |>
    dplyr::select("protein_id", "sample", "iBAQ") |>
    tidyr::pivot_wider(names_from = "protein_id", values_from = "iBAQ")

  sample_names <- ibaq_wide$sample
  mat <- as.matrix(ibaq_wide[, -1])
  rownames(mat) <- sample_names

  # Replace NA with 0 for PCA

  mat[is.na(mat)] <- 0

  n_samples <- nrow(mat)

  # ---- Check minimum sample requirement --------------------------------------
  if (n_samples < 3) {
    stop(
      "PCA requires at least 3 samples but only ",
      n_samples, " were found."
    )
  }

  # ---- Remove zero-variance columns (proteins) to avoid scaling issues -------
  col_vars <- apply(mat, 2, stats::var, na.rm = TRUE)
  mat <- mat[, col_vars > 0, drop = FALSE]

  if (ncol(mat) < 2) {
    stop("Fewer than 2 proteins with non-zero variance. Cannot perform PCA.")
  }

  # ---- PCA -------------------------------------------------------------------
  pca_res <- stats::prcomp(mat, center = TRUE, scale. = scale)

  max_pc <- max(components)
  if (max_pc > ncol(pca_res$x)) {
    stop(
      "Requested PC", max_pc,
      " but only ", ncol(pca_res$x),
      " components are available."
    )
  }

  var_explained <- summary(pca_res)$importance["Proportion of Variance", ]
  pc_x <- components[1]
  pc_y <- components[2]

  scores <- as.data.frame(pca_res$x[, components, drop = FALSE])
  scores$sample <- sample_names
  colnames(scores)[1:2] <- c("PC_x", "PC_y")

  # ---- KNN graph + Louvain clustering ----------------------------------------
  effective_k <- min(k, n_samples - 1)

  if (effective_k < 1) {
    # Only 1 or 2 samples — assign all to one cluster
    scores$cluster <- factor(1)
  } else {
    knn_res <- FNN::get.knn(pca_res$x, k = effective_k)

    # Build adjacency matrix from KNN results
    adj_mat <- matrix(0, nrow = n_samples, ncol = n_samples)
    for (i in seq_len(n_samples)) {
      neighbours <- knn_res$nn.index[i, ]
      adj_mat[i, neighbours] <- 1
    }
    # Make symmetric (mutual KNN)
    adj_mat <- pmax(adj_mat, t(adj_mat))

    g <- igraph::graph_from_adjacency_matrix(
      adj_mat,
      mode = "undirected",
      diag = FALSE
    )
    louvain_res <- igraph::cluster_louvain(g)
    scores$cluster <- factor(igraph::membership(louvain_res))
  }

  # ---- Axis labels -----------------------------------------------------------
  x_lab <- paste0(
    "PC", pc_x, " (",
    round(var_explained[pc_x] * 100, 1), "% variance)"
  )
  y_lab <- paste0(
    "PC", pc_y, " (",
    round(var_explained[pc_y] * 100, 1), "% variance)"
  )

  transform_note <- if (isTRUE(log2_values)) {
    "log\u2082-transformed iBAQ values"
  } else {
    "raw iBAQ values"
  }

  # ---- Build plot ------------------------------------------------------------
  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(
      x = .data$PC_x,
      y = .data$PC_y,
      colour = .data$cluster
    )
  ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(
      title = "PCA of iBAQ Protein Abundances",
      subtitle = paste0(
        "Louvain clustering on KNN graph (k=", effective_k,
        ") | ", transform_note
      ),
      x = x_lab,
      y = y_lab,
      colour = "Cluster"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40")
    )

  if (!is.null(color_palette)) {
    p <- p + ggplot2::scale_colour_manual(values = color_palette)
  } else {
    n_clusters <- length(unique(scores$cluster))
    default_pal <- grDevices::hcl.colors(
      max(n_clusters, 3),
      palette = "Dark 3"
    )
    p <- p + ggplot2::scale_colour_manual(values = default_pal)
  }

  # ---- Sample labels via ggrepel ---------------------------------------------
  if (isTRUE(label_samples)) {
    rlang::check_installed(
      "ggrepel",
      reason = "to label samples on the PCA plot."
    )
    p <- p +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = .data$sample),
        size = 3,
        colour = "grey20",
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        min.segment.length = 0,
        segment.colour = "grey60",
        segment.size = 0.2,
        show.legend = FALSE
      )
  }

  p
}
