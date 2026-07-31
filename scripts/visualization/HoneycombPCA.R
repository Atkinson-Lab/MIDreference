# ============================================================
# HoneycombPCA
#
# Given PLINK2 --pca output (.eigenvec / .eigenval), automatically
# choose the number of k-means clusters, assign clusters, and plot:
#   (1) PC1 vs PC2 colored by cluster, shaped by region/population
#   (2) a "honeycomb" circle-packing plot of cluster composition
#       faceted by region/population
#
# Region/population grouping per sample is the FID column, used as-is
# (no parsing). Supply `region_map` instead if you'd rather provide
# region/population from a separate lookup keyed by sample ID.
#
# Requires: dplyr, ggplot2, cluster, packcircles, cowplot
# ============================================================

.check_packages_honeycomb <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      ". Install with install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), "))"
    )
  }
}

#' Default colorblind-friendly palette, recycled/extended as needed
.default_honeycomb_colors <- function(n) {
  base_palette <- c(
    "#E69F00", "#009E73", "#56B4E9", "#882255", "#CC79A7",
    "#F0E442", "#003366", "#D55E00", "#0072B2", "#999999"
  )
  if (n <= length(base_palette)) {
    base_palette[1:n]
  } else {
    grDevices::colorRampPalette(base_palette)(n)
  }
}

#' Choose the number of k-means clusters via within-group SS + silhouette width
#'
#' Runs k-means for k = 1..kmax, computes total within-cluster sum of
#' squares for each, and uses PAM silhouette width (k = 2..kmax) to pick
#' the best k. Adapted from `00-DetermineK.R` (Constantinescu et al. 2022).
#'
#' @param mydata Numeric matrix/data frame (e.g. the first `pcs` PCs).
#' @param kmax Maximum number of clusters to try. Default 10.
#' @param num_iter `iter.max` passed to `kmeans()`. Default 100.
#' @param diagnostic_plot_file Optional file path to save a WSS +
#'   silhouette diagnostic plot to. If NULL, no diagnostic plot is made.
#' @return A list with `K` (all kmeans fits), `BestK` (chosen k), `WSS`,
#'   `Silhouette`, and `ClusteredData` (`mydata` with a `Cluster` column).
determine_k <- function(mydata, kmax = 10, num_iter = 100, diagnostic_plot_file = NULL) {
  .check_packages_honeycomb("cluster")
  
  K <- lapply(1:kmax, function(k) {
    stats::kmeans(mydata, centers = k, iter.max = num_iter, nstart = 10)
  })
  WSS <- sapply(K, function(fit) fit$tot.withinss)
  sil <- sapply(2:kmax, function(k) {
    cluster::pam(mydata, k)$silinfo$avg.width
  })
  k.sil <- which.max(sil) + 1
  
  if (!is.null(diagnostic_plot_file)) {
    grDevices::png(diagnostic_plot_file, width = 1000, height = 500)
    graphics::layout(matrix(c(1, 2), ncol = 2))
    graphics::plot(1:kmax, WSS, type = "b", pch = 19, frame = FALSE,
                   xlab = "Number of clusters K", ylab = "Tot. Within Clusters SS",
                   main = "Cluster Within-Group\nSum of Squares")
    graphics::abline(v = k.sil, lty = 2)
    graphics::plot(c(NA, sil), type = "b", pch = 19, frame = FALSE,
                   main = "Silhouette Analysis", xlab = "Number of clusters K",
                   ylab = "Avg Silhouette Width")
    graphics::abline(v = k.sil, lty = 2)
    grDevices::dev.off()
  }
  
  best_kmeans <- stats::kmeans(mydata, centers = k.sil, iter.max = num_iter, nstart = 10)
  mydata$Cluster <- best_kmeans$cluster
  
  list(K = K, BestK = k.sil, WSS = WSS, Silhouette = sil, ClusteredData = mydata)
}

#' Build the PC1 vs PC2 scatter plot, colored by cluster and shaped by region
plot_pca_clusters <- function(pca, pct_var, region_col = "cont",
                              colors = NULL, legend_position = c(0.75, 0.2)) {
  .check_packages_honeycomb("ggplot2")
  k <- nlevels(pca$Cluster)
  if (is.null(colors)) colors <- .default_honeycomb_colors(k)
  if (length(colors) < k) stop("`colors` must have at least ", k, " entries (one per cluster).")
  
  ggplot2::ggplot(pca, ggplot2::aes(x = PC1, y = PC2, color = Cluster, shape = .data[[region_col]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = paste0("PC1 (", round(pct_var[1], 2), "%)"),
      y = paste0("PC2 (", round(pct_var[2], 2), "%)"),
      shape = "Geography",
      color = "Cluster"
    ) +
    ggplot2::scale_color_manual(values = colors[1:k]) +
    ggplot2::scale_shape_manual(values = seq(0, nlevels(pca[[region_col]]))) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 2),
      color = ggplot2::guide_legend(order = 1)
    )
}

#' Build the honeycomb (circle-packing) plot of cluster composition per region
#'
#' One circle-packed panel per region, each circle a sample colored by its
#' cluster assignment, faceted with `facet_wrap(~region)`.
plot_honeycomb <- function(pca, region_col = "cont", region_order = NULL, colors = NULL) {
  .check_packages_honeycomb(c("dplyr", "ggplot2", "packcircles"))
  
  if (is.null(region_order)) {
    region_order <- sort(unique(as.character(pca[[region_col]])))
  }
  pca[[region_col]] <- factor(pca[[region_col]], levels = region_order)
  
  k <- nlevels(pca$Cluster)
  if (is.null(colors)) colors <- .default_honeycomb_colors(k)
  
  dat_list <- lapply(levels(pca[[region_col]]), function(reg) {
    subset_data <- pca[pca[[region_col]] == reg, c("Cluster", region_col)]
    subset_data <- subset_data[order(subset_data$Cluster), ]
    subset_data$value <- 1
    subset_data$id <- seq_len(nrow(subset_data))
    
    if (nrow(subset_data) == 0) return(NULL)
    
    packing <- packcircles::circleProgressiveLayout(subset_data$value, sizetype = "area")
    packing$theta <- atan2(packing$y, packing$x)
    packing <- packing[order(packing$theta), ]
    
    dat.gg <- packcircles::circleLayoutVertices(packing, idcol = NULL, npoints = 1)
    dat.gg <- merge(dat.gg, subset_data[, c("Cluster", region_col, "id")],
                    by.x = "id", by.y = "id", all.x = TRUE)
    dat.gg
  })
  names(dat_list) <- levels(pca[[region_col]])
  combined <- dplyr::bind_rows(dat_list)
  
  ggplot2::ggplot(combined) +
    ggplot2::geom_point(ggplot2::aes(x, y, colour = factor(Cluster))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 16)
    ) +
    ggplot2::scale_color_manual(values = colors[1:k]) +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", region_col)))
}

#' HoneycombPCA: cluster samples from PCA and plot PC scatter + honeycomb
#'
#' Reads a PLINK2 `.eigenvec`/`.eigenval` pair, chooses the number of
#' k-means clusters automatically (via within-group SS + silhouette
#' width), assigns each sample to a cluster, and produces:
#'   1. a PC1-vs-PC2 scatter plot colored by cluster, shaped by region
#'   2. a circle-packing ("honeycomb") plot of cluster composition per
#'      region
#'   3. a combined side-by-side figure of the two (optionally with a
#'      third panel supplied via `admixture_plot`, e.g. the `$plot`
#'      returned by `LabelFreeADMIX()`)
#'
#' Region/population per sample is the FID column, used as-is (no parsing);
#' supply `region_map` instead if you'd rather provide it explicitly from a
#' separate lookup keyed by sample ID.
#'
#' @param eigenvec_file Path to a PLINK2 `.eigenvec` file (has a header
#'   row starting with `#FID`).
#' @param eigenval_file Path to the matching `.eigenval` file (one
#'   eigenvalue per line, no header).
#' @param pcs Number of leading PCs to use for clustering. Default 5.
#' @param k Optional manual cluster count. If NULL (default), the number
#'   of clusters is chosen automatically via `determine_k()`.
#' @param kmax Max clusters to consider when choosing k automatically.
#' @param num_iter `iter.max` for k-means. Default 100.
#' @param region_map Optional data frame with a sample ID column matching
#'   `iid_col` and a region column, to use instead of the raw FID.
#' @param region_col Name to give the region column internally and in
#'   plots. Default "cont".
#' @param region_order Optional character vector giving the order regions
#'   should appear in the honeycomb facets/legend. If NULL, alphabetical.
#' @param fid_col,iid_col Column names for family/sample ID in the parsed
#'   eigenvec file. Defaults match PLINK2 output ("X.FID", "IID") — PLINK2
#'   writes the header as `#FID`, which R's `read.table` turns into
#'   `X.FID`.
#' @param colors Optional vector of cluster fill/point colors. If NULL, a
#'   default colorblind-friendly palette is used.
#' @param admixture_plot Optional ggplot object (e.g. from
#'   `LabelFreeADMIX()$plot`) to include as a third panel in the combined
#'   figure.
#' @param output_prefix Optional path prefix for saving outputs, e.g.
#'   "~/results/figures/allcohort". If set, saves
#'   `<prefix>_pc1pc2.png`, `<prefix>_honeycomb.png`,
#'   `<prefix>_combined.png`, `<prefix>_silhouette.png`, and
#'   `<prefix>.clusters.tsv`. If NULL, nothing is written to disk.
#' @param plot_width,plot_height Dimensions (inches) for the PC1/PC2 and
#'   honeycomb plots. Default 6x6 and 4x4.
#'
#' @return A list with:
#'   \item{pca}{the eigenvec data, with `region_col` and `Cluster` columns added}
#'   \item{clusters}{data frame of FID, IID, Cluster (as written to `.clusters.tsv`)}
#'   \item{k}{the number of clusters used}
#'   \item{pc_scatter}{ggplot object, PC1 vs PC2}
#'   \item{honeycomb}{ggplot object, circle-packing plot}
#'   \item{combined}{ggplot object combining the above (+ `admixture_plot` if given)}
#'   \item{determine_k}{full output of `determine_k()` (WSS, silhouette, etc.), or NULL if `k` was set manually}
#'
#' @examples
#' \dontrun{
#' result <- HoneycombPCA(
#'   eigenvec_file = "allcohorts.eigenvec",
#'   eigenval_file = "allcohorts.eigenval",
#'   pcs = 5,
#'   output_prefix = "~/results/figures/allcohort"
#' )
#' result$combined
#' }
HoneycombPCA <- function(eigenvec_file, eigenval_file, pcs = 5,
                         k = NULL, kmax = 10, num_iter = 100,
                         region_map = NULL, region_col = "cont", region_order = NULL,
                         fid_col = "X.FID", iid_col = "IID",
                         colors = NULL, admixture_plot = NULL,
                         output_prefix = NULL,
                         plot_width = 6, plot_height = 6) {
  .check_packages_honeycomb(c("dplyr", "ggplot2", "cluster", "packcircles", "cowplot"))
  
  # --- Read PCA output ---
  pca <- utils::read.table(eigenvec_file, header = TRUE, comment.char = "",
                           check.names = TRUE, stringsAsFactors = FALSE)
  if (!fid_col %in% colnames(pca)) {
    stop("`fid_col` ('", fid_col, "') not found in eigenvec file. Columns present: ",
         paste(colnames(pca), collapse = ", "))
  }
  if (!iid_col %in% colnames(pca)) {
    stop("`iid_col` ('", iid_col, "') not found in eigenvec file. Columns present: ",
         paste(colnames(pca), collapse = ", "))
  }
  
  pc_cols <- paste0("PC", 1:pcs)
  missing_pcs <- setdiff(pc_cols, colnames(pca))
  if (length(missing_pcs) > 0) {
    stop("Requested `pcs` = ", pcs, " but these columns are missing from the ",
         "eigenvec file: ", paste(missing_pcs, collapse = ", "))
  }
  
  eigenvalues <- scan(eigenval_file, quiet = TRUE)
  pct_var <- (eigenvalues / sum(eigenvalues)) * 100
  
  # --- Region/population grouping ---
  if (!is.null(region_map)) {
    if (!iid_col %in% colnames(region_map)) {
      stop("`region_map` must contain an ID column matching `iid_col` ('", iid_col, "').")
    }
    map_region_col <- setdiff(colnames(region_map), iid_col)[1]
    pca[[region_col]] <- region_map[[map_region_col]][match(pca[[iid_col]], region_map[[iid_col]])]
    if (any(is.na(pca[[region_col]]))) {
      warning(sum(is.na(pca[[region_col]])), " sample(s) had no match in `region_map`.")
    }
  } else {
    pca[[region_col]] <- as.character(pca[[fid_col]])
  }
  pca[[region_col]] <- as.factor(pca[[region_col]])
  
  # --- Choose k and assign clusters ---
  included_pcs <- pca[, pc_cols, drop = FALSE]
  diag_file <- if (!is.null(output_prefix)) paste0(output_prefix, "_silhouette.png") else NULL
  
  if (is.null(k)) {
    dk <- determine_k(included_pcs, kmax = kmax, num_iter = num_iter, diagnostic_plot_file = diag_file)
    k <- dk$BestK
    pca$Cluster <- as.factor(dk$ClusteredData$Cluster)
  } else {
    dk <- NULL
    km <- stats::kmeans(included_pcs, centers = k, iter.max = num_iter, nstart = 10)
    pca$Cluster <- as.factor(km$cluster)
  }
  
  # --- Cluster assignment table ---
  clusters <- pca[, c(fid_col, iid_col, "Cluster")]
  colnames(clusters) <- c("#FID", "IID", "Cluster")
  
  # --- Plots ---
  pc_scatter <- plot_pca_clusters(pca, pct_var, region_col = region_col, colors = colors)
  honeycomb <- plot_honeycomb(pca, region_col = region_col, region_order = region_order, colors = colors)
  
  if (!is.null(admixture_plot)) {
    combined <- cowplot::plot_grid(
      cowplot::plot_grid(pc_scatter, honeycomb, labels = c("(a)", "(b)"), label_size = 16,
                         ncol = 2, align = "hv", axis = "tb", rel_widths = c(1.3, 1), label_x = -0.035),
      cowplot::plot_grid(admixture_plot, labels = "(c)", label_size = 16),
      ncol = 1, rel_heights = c(1.1, 1)
    )
  } else {
    combined <- cowplot::plot_grid(pc_scatter, honeycomb, labels = c("(a)", "(b)"),
                                   label_size = 16, ncol = 2, align = "hv", axis = "tb",
                                   rel_widths = c(1.3, 1), label_x = -0.035)
  }
  
  # --- Save outputs ---
  if (!is.null(output_prefix)) {
    write.table(clusters, paste0(output_prefix, ".clusters.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    ggplot2::ggsave(paste0(output_prefix, "_pc1pc2.png"), plot = pc_scatter,
                    width = plot_width, height = plot_height)
    ggplot2::ggsave(paste0(output_prefix, "_honeycomb.png"), plot = honeycomb,
                    width = plot_width * 2/3, height = plot_height * 2/3)
    combined_height <- if (!is.null(admixture_plot)) plot_height * 2.2 else plot_height * 1.1
    ggplot2::ggsave(paste0(output_prefix, "_combined.png"), plot = combined,
                    width = plot_width * 1.6, height = combined_height)
  }
  
  list(
    pca = pca,
    clusters = clusters,
    k = k,
    pc_scatter = pc_scatter,
    honeycomb = honeycomb,
    combined = combined,
    determine_k = dk
  )
}