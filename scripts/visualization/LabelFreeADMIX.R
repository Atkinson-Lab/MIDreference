# ============================================================
# LabelFreeADMIX
#
# Order and plot ADMIXTURE Q-matrix results without specifying
# population labels. Region grouping lines/labels are optional
# annotations on top of the plot;
# sample *order* is inferred automatically unless overridden.
#
# Requires: dplyr, reshape2, ggplot2
# ============================================================

.check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      ". Install with install.packages(c(", 
      paste(sprintf('"%s"', missing), collapse = ", "), "))"
    )
  }
}

#' Default colorblind-friendly palette (Okabe-Ito based), recycled/extended as needed
#'
#' @param k Number of colors needed
#' @return Character vector of length k
.default_admix_colors <- function(k) {
  base_palette <- c(
    "#CC79A7", "#882255", "#F0E442", "#009E73",
    "#E69F00", "#56B4E9", "#D55E00", "#332288"
  )
  if (k <= length(base_palette)) {
    base_palette[1:k]
  } else {
    grDevices::colorRampPalette(base_palette)(k)
  }
}

#' Identify the background component and order remaining components
#'
#' Implements the component-ordering logic: the "background" component is
#' the one with the fewest samples above the minimum-max-proportion threshold
#' (i.e. the component that is rarely anyone's dominant ancestry). The
#' remaining components are ordered left-to-right starting from the one with
#' the most highly-admixed samples (proportion > 0.8), then chained forward
#' by shared ancestry with the previous component.
#'
#' @param q Data frame of ADMIXTURE proportions; one row per sample.
#' @param component_cols Character vector of column names in `q` holding the
#'   per-component proportions (should sum to ~1 per row).
#' @return A list with `background` (component name) and `indexes` (ordered
#'   character vector of the remaining k-1 components, left to right).
order_admixture_components <- function(q, component_cols) {
  k <- length(component_cols)
  qk <- q[, component_cols, drop = FALSE]
  
  # 1.1 Minimum of each sample's maximum component (threshold for "admixed")
  min_prop <- min(apply(qk, 1, max))
  remaining_columns <- component_cols
  
  # 1.2 Background = component with fewest samples above min_prop
  background <- names(which.min(sapply(qk, function(x) sum(x < min_prop))))
  remaining_columns <- setdiff(remaining_columns, background)
  
  # 1.3 First plotted component: most samples with proportion > 0.8
  comp1 <- names(which.max(sapply(qk[, remaining_columns, drop = FALSE], function(x) sum(x > 0.8))))
  remaining_columns <- setdiff(remaining_columns, comp1)
  indexes <- comp1
  
  # 1.4 Chain remaining components by shared ancestry with the previous one
  while (length(remaining_columns) > 1) {
    logical_matrix <- qk[, remaining_columns, drop = FALSE][qk[[comp1]] > min_prop, , drop = FALSE]
    comp1 <- names(sort(colSums(logical_matrix), decreasing = TRUE))[1]
    indexes <- c(indexes, comp1)
    remaining_columns <- setdiff(remaining_columns, comp1)
  }
  # 1.5 Add the last remaining component
  indexes <- c(indexes, remaining_columns[1])
  
  list(background = background, indexes = indexes)
}

#' Order samples for plotting based on their dominant (non-background) component
#'
#' For each component (in the order given by `indexes`), samples whose
#' dominant component is that one are grouped together and sorted by their
#' proportion of that component (alternating ascending/descending so
#' neighboring groups blend smoothly in the plot).
#'
#' @param q Data frame of ADMIXTURE proportions, with `id_col` and
#'   `region_col` columns present.
#' @param component_cols Character vector of proportion column names.
#' @param background Name of the background component (excluded as a
#'   sample's "dominant" component; the 2nd-highest is used instead if a
#'   sample's top component is background).
#' @param indexes Ordered character vector of non-background components
#'   (output of `order_admixture_components()`).
#' @param id_col Name of the sample ID column in `q`. Default "id".
#' @param region_col Name of the region/grouping column in `q`. Default "region".
#' @return `q` reordered by group and within-group proportion.
order_admixture_samples <- function(q, component_cols, background, indexes,
                                    id_col = "id", region_col = "region") {
  k <- length(component_cols)
  
  # 2.1 Determine each sample's dominant component (skip background)
  max_one <- matrix(0, nrow = nrow(q), ncol = k,
                    dimnames = list(rownames(q), component_cols))
  for (i in seq_len(nrow(q))) {
    x <- as.numeric(q[i, component_cols])
    ordered_components <- component_cols[order(x, decreasing = TRUE)]
    max_index <- ordered_components[1]
    if (max_index == background) {
      max_index <- ordered_components[2]
    }
    max_one[i, max_index] <- 1
  }
  max_one <- as.data.frame(max_one)
  
  # 2.2 Bin samples by dominant component (in `indexes` order), sort within bin
  sample_order <- data.frame()
  for (i in seq_along(indexes)) {
    comp1 <- indexes[i]
    bin_indices <- which(max_one[[comp1]] == 1)
    bin <- q[bin_indices, ]
    
    if (i == length(indexes)) {
      bin <- bin[order(bin[[comp1]]), ]
    } else if (i %% 2 == 1) {
      bin <- bin[rev(order(bin[[comp1]])), ]
    } else {
      bin <- bin[order(bin[[comp1]]), ]
    }
    sample_order <- rbind(sample_order, bin)
  }
  sample_order
}

#' Infer a display order for regions from each sample's dominant component
#'
#' For each component in `indexes` (left to right), finds which region(s)
#' are most often assigned that component as their dominant ancestry, and
#' appends any not-yet-seen regions to the order. This removes the need to
#' hand-specify a `region_order` vector.
#'
#' @param major_component_df Data frame with a `Major component` column and
#'   a region column.
#' @param indexes Ordered character vector of non-background components.
#' @param region_col Name of the region column. Default "region".
#' @return Character vector giving the inferred region display order.
infer_region_order <- function(major_component_df, indexes, region_col = "region") {
  region_order <- c()
  for (comp in indexes) {
    tab <- major_component_df[major_component_df[["Major component"]] == comp, ]
    counts <- sort(table(tab[[region_col]]), decreasing = TRUE)
    new_regions <- setdiff(names(counts), region_order)
    region_order <- c(region_order, new_regions)
  }
  # append any regions never seen as a majority component (edge case safety net)
  all_regions <- unique(major_component_df[[region_col]])
  unique(c(region_order, all_regions))
}

#' Build the ADMIXTURE bar plot with region grouping lines/labels
#'
#' @param sample_order Data frame of samples in plotting order (output of
#'   `order_admixture_samples()`).
#' @param component_cols Character vector of proportion column names.
#' @param region_order Character vector giving top-to-bottom region label
#'   order. If NULL, inferred automatically via `infer_region_order()`.
#' @param colors Character vector of fill colors, one per component
#'   (length must equal `length(component_cols)`). If NULL, a default
#'   colorblind-friendly palette is used.
#' @param id_col Name of the sample ID column. Default "id".
#' @param region_col Name of the region column. Default "region".
#' @param label_offset Horizontal offset (in sample units) for region text
#'   labels to the left of the plot. Default 8.
#' @param k_label Label for the "k=" annotation on the left of the plot.
#'   If NULL, uses `length(component_cols)`.
#' @return A ggplot object.
plot_labelfree_admixture <- function(sample_order, component_cols,
                                     region_order = NULL, colors = NULL,
                                     id_col = "id", region_col = "region",
                                     label_offset = 8, k_label = NULL) {
  .check_packages(c("dplyr", "reshape2", "ggplot2"))
  k <- length(component_cols)
  if (is.null(colors)) colors <- .default_admix_colors(k)
  if (length(colors) != k) stop("`colors` must have length ", k, " (one per component).")
  if (is.null(k_label)) k_label <- k
  
  admixture_df <- reshape2::melt(
    sample_order,
    id.vars = c(id_col, region_col),
    measure.vars = component_cols,
    variable.name = "Cluster",
    value.name = "Proportion"
  )
  admixture_df$order <- factor(admixture_df[[id_col]], levels = sample_order[[id_col]])
  
  if (is.null(region_order)) {
    region_order <- unique(sample_order[[region_col]])
  }
  
  spacing <- 0.2
  region_line_data <- sample_order
  region_line_data[[region_col]] <- factor(region_line_data[[region_col]], levels = region_order)
  region_line_data <- dplyr::group_by(region_line_data, .data[[region_col]])
  region_line_data <- dplyr::summarise(
    region_line_data,
    x_start = which(admixture_df$order == dplyr::first(.data[[id_col]]))[1],
    x_end   = which(admixture_df$order == dplyr::last(.data[[id_col]]))[1],
    .groups = "drop"
  )
  region_line_data <- dplyr::mutate(region_line_data, region_y = -as.numeric(.data[[region_col]]) * spacing)
  
  sample_lines <- sample_order
  sample_lines[[region_col]] <- factor(sample_lines[[region_col]], levels = region_order)
  sample_lines <- dplyr::left_join(sample_lines, region_line_data, by = region_col)
  
  p <- ggplot2::ggplot(admixture_df, ggplot2::aes(x = order, y = Proportion, fill = Cluster)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::geom_segment(
      data = region_line_data,
      ggplot2::aes(x = x_start, xend = x_end, y = region_y, yend = region_y),
      color = "gray", linetype = "solid", linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = sample_lines,
      ggplot2::aes(x = .data[[id_col]], xend = .data[[id_col]],
                   y = region_y - 0.05, yend = region_y + 0.05),
      color = "black", linetype = "solid", linewidth = 0.6, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = region_line_data,
      ggplot2::aes(x = x_start - label_offset, y = region_y, label = .data[[region_col]]),
      color = "black", size = 5, hjust = 1, vjust = 0.5, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text", x = 0, y = max(admixture_df$Proportion) / 2,
      label = paste0("k=", k_label), size = 5, color = "black", angle = 0, hjust = 1.5
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.1, 0.1)),
      limits = c(min(region_line_data$region_y) - 0.5, max(admixture_df$Proportion) + 0.5)
    ) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 5, 1, 50),
      text = ggplot2::element_text(size = 18),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  p
}

#' Summarize each region's distribution across major (dominant) components
#'
#' @param q Data frame of ADMIXTURE proportions with `region_col`.
#' @param component_cols Character vector of proportion column names.
#' @param region_col Name of the region column. Default "region".
#' @return Data frame: Major component, Region, Sample count, Proportion of
#'   total region samples.
summarize_admixture_components <- function(q, component_cols, region_col = "region") {
  .check_packages("dplyr")
  major <- component_cols[max.col(q[, component_cols, drop = FALSE], ties.method = "first")]
  q$`Major component` <- major
  
  by_comp_region <- dplyr::count(q, `Major component`, .data[[region_col]], name = "count")
  by_region <- dplyr::count(q, .data[[region_col]], name = "total")
  results <- dplyr::left_join(by_comp_region, by_region, by = region_col)
  results <- dplyr::mutate(results, proportion = round(count / total, 2))
  results <- dplyr::select(results, `Major component`, dplyr::all_of(region_col), count, proportion)
  
  colnames(results) <- c("Major component", "Region", "Sample count",
                         "Proportion of total region samples (%)")
  results
}

#' LabelFreeADMIX: order and plot ADMIXTURE results without population labels
#'
#' Given an ADMIXTURE Q-matrix (proportions per sample per component), this
#' automatically identifies the "background" component,
#' determines a left-to-right plotting order for the remaining components
#' based on shared ancestry, orders samples for a clean stacked-bar plot,
#' and (optionally) infers a top-to-bottom order for population/region group labels —
#' all without hand-specifying population labels up front.
#'
#' @param q Data frame with one row per sample: k columns of ADMIXTURE
#'   proportions (should sum to ~1 per row) plus an ID column and a region
#'   column. Component columns are inferred as every column except
#'   `id_col` and `region_col`.
#' @param id_col Name of the sample ID column. Default "id".
#' @param region_col Name of the region/population column. Default "region".
#' @param region_order Optional character vector giving a manual top-to-
#'   bottom order for region labels. If NULL (default), the order is
#'   inferred automatically from which region most often carries each
#'   component as its dominant ancestry.
#' @param colors Optional character vector of fill colors, one per
#'   component. If NULL, a default colorblind-friendly palette is used.
#' @param label_offset Horizontal offset for region text labels. Default 8.
#' @param plot_file Optional file path to save the plot to (via ggsave).
#' @param plot_width,plot_height Dimensions for `ggsave`, if `plot_file` is set.
#' @param table_file Optional file path to save the component/region
#'   summary table to (tab-separated).
#'
#' @return A list with:
#'   \item{plot}{the ggplot object}
#'   \item{sample_order}{`q` reordered for plotting}
#'   \item{component_order}{ordered non-background component names}
#'   \item{background}{the identified background component}
#'   \item{region_order}{the region label order used}
#'   \item{summary}{per-region major-component summary table}
#'
#' @examples
#' \dontrun{
#' result <- LabelFreeADMIX(q, id_col = "id", region_col = "region")
#' result$plot
#' ggsave("admixture.png", result$plot, width = 9, height = 6)
#' }
LabelFreeADMIX <- function(q, id_col = "id", region_col = "region",
                           region_order = NULL, colors = NULL,
                           label_offset = 8,
                           plot_file = NULL, plot_width = 9, plot_height = 6,
                           table_file = NULL) {
  .check_packages(c("dplyr", "reshape2", "ggplot2"))
  
  if (!all(c(id_col, region_col) %in% colnames(q))) {
    stop("`q` must contain both `id_col` ('", id_col, "') and `region_col` ('", region_col, "') columns.")
  }
  component_cols <- setdiff(colnames(q), c(id_col, region_col))
  if (length(component_cols) < 3) {
    stop("Found only ", length(component_cols), " component column(s) after excluding ",
         "id_col/region_col — check that `q` has one column per ADMIXTURE component.")
  }
  
  ordering <- order_admixture_components(q, component_cols)
  sample_order <- order_admixture_samples(
    q, component_cols, ordering$background, ordering$indexes,
    id_col = id_col, region_col = region_col
  )
  
  summary_table <- summarize_admixture_components(q, component_cols, region_col = region_col)
  
  if (is.null(region_order)) {
    region_order <- infer_region_order(
      data.frame(
        `Major component` = summary_table$`Major component`,
        region = summary_table$Region,
        check.names = FALSE
      ),
      ordering$indexes,
      region_col = "region"
    )
  }
  
  p <- plot_labelfree_admixture(
    sample_order, component_cols,
    region_order = region_order, colors = colors,
    id_col = id_col, region_col = region_col,
    label_offset = label_offset, k_label = length(component_cols)
  )
  
  if (!is.null(plot_file)) {
    ggplot2::ggsave(plot_file, plot = p, width = plot_width, height = plot_height)
  }
  if (!is.null(table_file)) {
    write.table(summary_table, table_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = TRUE)
  }
  
  list(
    plot = p,
    sample_order = sample_order,
    component_order = ordering$indexes,
    background = ordering$background,
    region_order = region_order,
    summary = summary_table
  )
}