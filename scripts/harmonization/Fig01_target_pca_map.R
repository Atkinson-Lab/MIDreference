# ==========================================================

# Author: Hatoon AlAli
# Paper: An Expanded Reference for Middle Eastern Populations Improves GWAS
# Generates Figure 1 (PCA of different MID definitions + map)
# R version: 4.4.1

# ==========================================================
# 0. Libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(maps)
library(scatterpie)

# 1. Paths & Parameters
input = "~/midProject/results/raw_outputs/pca/targetcohort"
output = "~/midProject/results/figures/harmonization/pca/"

set.seed(1234)

# ==========================================================
# PART I — PCA
# ==========================================================

# ----------------------------------------------------------
# READ INPUTS
# ----------------------------------------------------------

# Variance file
var <- read.csv(paste0(input, ".eigenval"), header = F, sep = "\t")
# PCA file
pca <- read.csv(paste0(input, ".eigenvec"), header = T, sep = "\t")

pca$id <- as.factor(gsub("^.*]", "", pca$X.FID))
pca$id <- as.factor(gsub("^.*\\|", "", pca$id))
pca$eth <- as.factor(gsub("\\|.*", "", pca$X.FID))
pca$eth <- as.factor(gsub(".*]", "", pca$eth))
pca$cont <- as.factor(gsub("^.*_", "", pca$id))
pca$id <- as.factor(gsub("_.*", "", pca$id))
pca$study <- as.factor(gsub("].*", "", pca$X.FID))
pca$study <- as.factor(gsub("\\[", "", pca$study))

pca$eth <- gsub("^N$", "MoroccoNorth",  pca$eth)
pca$eth <- gsub("^S$", "MoroccoSouth", pca$eth)

# ----------------------------------------------------------
# List definitions
# ----------------------------------------------------------

# definitions
mid1 <- c("Emirates", "Iraq", 
          "Jordan", "Oman", 
          "MID", "Qatar", 
          "SaudiArabia", "Syria", 
          "Yemen", "Kuwait", "Lebanon")
mid2 <- c(mid1, "Egypt", "Libya", "Algeria", "Tunisia", "Morocco", "Sahara")
mid3 <- c(mid1, "Turkey", "Iran", "Cyprus", "Armenia")
mid4 <- c(mid2, "Turkey", "Iran", "Cyprus", "Armenia", "Afghanistan", "CSA") #Greater Middle East
mid5 <- c(mid2, "Sudan", "Comoros", "Somalia") # Arab League countries
#Arab language/ethnicity
eth <- c("EmiratiA", "NA", "Palestinian", "QPRC", "Arab", "Bedouin", "EmiratiB", "Jordan", "Saharawi",
         "Syria", "Yemen.o", "ArabBaggara", "EmiratiC", "MoroccoSouth", "MoroccoNorth", 
         "saudi", "syrian", "yemeni", "ArabKababish", "Druze", "EmiratiD", "kuwait", 
         "Saudi.o", "ArabRashaayda", "egyptian", "lebanese", "Oman", "QBC", "SaudiA", "tunisian", "uae",
         "Arabs", "Emirati.o", "Oman.o", "QIGR", "SaudiB")


definitions <- list(mid1, mid2, mid3, mid4, mid5, eth)

# ----------------------------------------------------------
# Scree plot
# ----------------------------------------------------------

# Save PCA variance 
colnames(var) <- "Variance"
var$Percentage <- round(var$Variance/sum(var$Variance)*100, 2)
var$PC <- 1:nrow(var)
# scree plot
scree_plot <- ggplot(var, aes(x = PC, y = Percentage)) +
  geom_point(size = 3, color = "blue") +
  geom_line(group = 1, color = "blue") +
  geom_hline(yintercept = 4, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = unique(var$PC)) + 
  labs(title = "Scree Plot", x = "PC", y = "Variance Explained (%)") +
  theme_classic()
scree_plot
ggsave(paste0(output, "targetpca_scree.png"), plot = scree_plot, width = 6, height = 4, dpi = 300)

# ----------------------------------------------------------
# Plot PC1 + PC2
# ----------------------------------------------------------

## Colors
colors <- c("#009E73",  # Green mid1
            "#56B4E9",  # Sky Blue mid2
            "#E69F00",  # Orange mid3
            "#CC79A7",  # Reddish Purple mid4
            "#F0E442",  # Yellow mid5
            "mediumpurple1")  # Deep Purple mid6

special_color <- "black"  # HGDP-MID1 and HGDP-MID2

# Assign colors to definitions
definitions <- list(mid1, mid2, mid3, mid4, mid5, eth)

define_shapes <- function(pca, definition, def_color, special_color) {
  pca$color <- ifelse(pca$id %in% c("MID"), special_color,
                      ifelse(pca$id %in% definition | pca$eth %in% definition, def_color, "#C8C8C8"))
  pca$shape <- ifelse(pca$id %in% c("MID"), 43,
                      ifelse(pca$id %in% definition | pca$eth %in% definition, 20, 1))
  pca$layer <- ifelse(pca$id %in% c("MID"), 3,
                      ifelse(pca$id %in% definition | pca$eth %in% definition, 2, 1))
  return(pca)
}

# Generate plots
fir_pc <- 1
sec_pc <- 2

pc_x <- paste0("PC", fir_pc)
pc_y <- paste0("PC", sec_pc)

plots <- list()
for (i in 1:length(definitions)) {
  pca_colored <- define_shapes(pca, definitions[[i]], colors[i], special_color)
  pca_colored <- pca_colored[order(pca_colored$layer), ]
  
  p <- ggplot(pca_colored, 
              aes(x = .data[[pc_x]], 
                  y = .data[[pc_y]], 
                  color = color, 
                  shape = factor(shape))) +
    geom_point() +
    scale_color_identity() +
    scale_shape_manual(values = c(20, 1, 43)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plots[[i]] <- p
}

# Legend data frame
legend_data <- data.frame(
  label = c("HGDP-MID", "MID1", "MID2", "MID3", "MID4", "MID5", "MID6"),
  color = c(special_color, colors[1:6]),
  shape = c(43, 20, 20, 20, 20, 20, 20)
)
legend_data$label <- factor(legend_data$label, 
                            levels = c("HGDP-MID", "MID1", "MID2", "MID3", "MID4", "MID5", "MID6"))

# Dummy data for legend
dummy_data <- data.frame(x = 1, y = 1)

legend_plot <- ggplot(dummy_data, aes(x, y)) +
  geom_point(data = legend_data,
             aes(x = 1, y = 1, color = label, shape = label),
             size = 5) +
  scale_color_manual(values = setNames(legend_data$color, legend_data$label)) +
  scale_shape_manual(values = setNames(legend_data$shape, legend_data$label)) +
  guides(color = guide_legend(title = "Definitions"),
         shape = guide_legend(title = "Definitions")) +
  theme_void() +
  theme(legend.position = "right",
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

legend <- cowplot::get_legend(legend_plot)

# Combine plots + legend
combined_plots <- arrangeGrob(
  do.call(arrangeGrob, c(plots, ncol = 3)),
  legend,
  ncol = 2,
  widths = c(2.7, 0.7)
)

# Axis labels with PC numbers
p <- grid.arrange(
  combined_plots,
  bottom = textGrob(paste0("PC", fir_pc, " (", round(var[fir_pc, 2], 2), "%)"),
                    gp = gpar(fontsize = 18)),
  left = textGrob(paste0("PC", sec_pc, " (", round(var[sec_pc, 2], 2), "%)"),
                  rot = 90, gp = gpar(fontsize = 18))
)

p

# Save figure
ggsave(paste0(output, "definitions_pc", fir_pc, "pc", sec_pc, ".png"),
       p, width = 8, height = 5)
# ==========================================================
# PART II — MAP
# ==========================================================

# add map and save combined plot (only for PC1 and PC2)
mapcolors <- c("#F0E442", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#999999", "purple", "lightgrey")

# Define the map ranges
lon_range <- c(-20, 80)
lat_range <- c(-0, 50)

# Create a base map
base_map <- map_data("world")
base_map <- subset(base_map, long >= lon_range[1] & long <= lon_range[2] &
                     lat >= lat_range[1] & lat <= lat_range[2])

# definitions
mid1 <- c("United Arab Emirates", "Iraq", "Jordan", "Oman", "Palestine", "Qatar", "Saudi Arabia", "Syria", "Yemen", "Israel")
mid2 <- c(mid1, "Algeria", "Egypt", "Morocco", "Western Sahara", "Tunisia", "Libya")
mid3 <- c(mid1, "Iran", "Turkey", "Cyprus")
mid4 <- c(mid2, "Iran", "Turkey", "Cyprus", "Afghanistan", "Pakistan")
mid5 <- c(mid2, "Sudan", "Comoros", "Somalia", "Mauritania", "Djibouti", "Eritrea")

# Define regions and colors
regions_list <- list(mid5, mid4, mid3, mid2, mid1)

# Create an empty ggplot object
map <- ggplot() +
  geom_polygon(data = base_map, aes(x = long, y = lat, group = group), fill = "lightgrey", color = NA) +
  coord_fixed(xlim = lon_range, ylim = lat_range) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), legend.position = "none")

# Increase the distance between layers by adjusting the nudge value
distance_between_layers <- 1.3  # Adjust this value to increase/decrease distance

# Add each map layer
for (i in seq_along(regions_list)) {
  map_layer <- map_data("world", regions = regions_list[[i]])
  map_layer <- subset(map_layer, long >= lon_range[1] & long <= lon_range[2] &
                        lat >= lat_range[1] & lat <= lat_range[2])
  map <- map + geom_polygon(data = map_layer, aes(x = long, y = lat, group = group),
                        fill = mapcolors[i], color = NA,
                        position = position_nudge(y = i * distance_between_layers))  # Adjust `y` for layering effect
}


label_df <- data.frame(
  region = c("SWAsia", "ECEurope", "ArabPen", "NArab", "EAfrica", "SEurope", "NAfrica", "CAfrica", "WAfrica"),
  lon = c(67.195, 45.87, 45.57, 39.60, 36.74, 27, 11.32, 17.11, -3.29),
  lat = c(29.97, 48, 28.11, 38.13, 11.55, 43.81, 35.5, 14.51, 15.34)
)

for (i in seq_along(label_df$region)) {
  lbl <- label_df[i, ]
  map <- map +
    geom_text(
      data = lbl,
      aes(x = lon, y = lat, label = region),
      size = 4
    )
}
# Display the plot
map
ggsave(paste0(output, "map.png"), map, width = 8, height = 3)

# ==========================================================
# PART III — Clinical variant distribtion
# ==========================================================
snps <- read.table("~/midProject/results/tables/clinical_snps_regions.tsv", header = TRUE, sep="\t", quote = "", fill = TRUE)
selected_snp = c("chr1:155313038:A:G")
snps <- snps[snps$SNP %in% selected_snp, ]
snps$id[snps$SNP == "chr1:155313038:A:G"] <- "rs11264359"

# build pie data (long → wide for scatterpie)
pie_df <- snps %>%
  tidyr::pivot_longer(
    cols = matches("\\.(MAF|MAC|NCHROBS)$"),
    names_to = c("region", ".value"),
    names_pattern = "(.*)\\.(MAF|MAC|NCHROBS)"
  ) %>%
  left_join(label_df, by = "region") %>%
  mutate(SNP = factor(SNP, levels = selected_snp))

# spread SNPs into columns for pie slices
pie_wide <- pie_df %>%
  select(region, lon, lat, SNP, MAF) %>%
  tidyr::pivot_wider(names_from = SNP, values_from = MAF)

pie_wide[, selected_snp] <- lapply(pie_wide[, selected_snp, drop = FALSE], as.numeric)
pie_wide[[paste0(selected_snp[1], "_ref")]] <- 1 - pie_wide[[selected_snp[1]]]

pie_wide1 <- pie_wide %>% mutate(x = lon, y = lat+3.5)

map_snps <- map +
  # pies
  geom_scatterpie(data = pie_wide1,
                  aes(x = x, y = y),
                  cols = c(selected_snp[1], paste0(selected_snp[1], "_ref")),
                  pie_scale = 1.5) +
  # add rs ID
  geom_text(data = data.frame(x = min(base_map$lon), 
                              y = max(base_map$lat)+4,
                              label = snps$id),
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 4) +
  
  coord_equal()
map_snps
ggsave(paste0(output, "map_snp.png"), map_snps, width = 8, height = 3)
# ==========================================================
# PART IV — FIGURE 1
# ==========================================================

# Add labels to each plot
p_labeled   <- ggdraw(p)   + draw_plot_label("(a)", x = 0.01, y = 0.98, hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
map_labeled <- ggdraw(map_snps) + draw_plot_label("(b)", x = 0.01, y = 0.98, hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Stack vertically with full width
final_plot <- plot_grid(
  p_labeled,
  map_labeled,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1)
)

final_plot

ggsave(paste0(output, "fig1.png"), final_plot, width = 6, height = 7)
# ============================= END =============================

