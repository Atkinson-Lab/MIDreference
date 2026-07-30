# ==========================================================

# Author: Hatoon AlAli
# Paper: An Expanded Reference for Middle Eastern Populations Improves GWAS
# Generates Supplementary Fig. 3 (PCA labeled by study, country, or ethincity)
# R version: 4.4.1

# ==========================================================
library(paletteer)
library(ggplot2)
library(cowplot)

# Read Variance file
var <- read.csv("~/midProject/results/raw_outputs/pca/targetcohort.eigenval", header = F, sep = "\t")
## Read PCA file
pca <- read.csv("~/midProject/results/raw_outputs/pca/targetcohort.eigenvec", header = T, sep = "\t")
output = "~/midProject/results/figures/harmonization/pca/"

colnames(var) <- "Variance"
var$Percentage <- round(var$Variance/sum(var$Variance)*100, 2)
var$PC <- 1:nrow(var)

pca$nat <- as.factor(gsub("^.*]", "", pca$X.FID))
pca$nat <- as.factor(gsub("^.*\\|", "", pca$nat))
pca$eth <- as.factor(gsub("\\|.*", "", pca$X.FID))
pca$eth <- as.factor(gsub(".*]", "", pca$eth))
pca$cont <- as.factor(gsub("^.*_", "", pca$nat))
pca$nat <- as.factor(gsub("_.*", "", pca$nat))
pca$study <- as.factor(gsub("].*", "", pca$X.FID))
pca$study <- as.factor(gsub("\\[", "", pca$study))

pca$eth <- gsub("^N$", "MoroccoNorth",  pca$eth)
pca$eth <- gsub("^S$", "MoroccoSouth", pca$eth)

# correct eth labels
pca$eth <- gsub("Egypt_NAfrica", "Egyptian", pca$eth)
pca$eth <- gsub("Afghanistan_CWAsia", "Afghan", pca$eth)
pca$eth <- gsub("Algeria_NAfrica", "Algerian", pca$eth)
pca$eth <- gsub("Tunisia_NAfrica", "Tunisian", pca$eth)
pca$eth <- gsub("Yemen_ArabPen", "Yemeni", pca$eth)
pca$eth <- gsub("Yemen.o", "Yemeni", pca$eth)
pca$eth <- gsub("yemeni", "Yemeni", pca$eth)
pca$eth <- gsub("Yemen", "Yemeni", pca$eth)
pca$eth <- gsub("Yemenii", "Yemeni", pca$eth)
pca$eth <- gsub("Oman.o", "Omani", pca$eth)
pca$eth <- gsub("Oman", "Omani", pca$eth)
pca$eth <- gsub("Omanii", "Omani", pca$eth)
pca$eth <- gsub("moroccan", "Moroccan", pca$eth)
pca$eth <- gsub("Morocco_NAfrica", "Moroccan", pca$eth)
pca$eth <- gsub("Libya_NAfrica", "Libyan", pca$eth)
pca$eth <- gsub("Saudi.o", "Saudi", pca$eth)
pca$eth <- gsub("saudi", "Saudi", pca$eth)
pca$eth <- gsub("Syria", "Syrian", pca$eth)
pca$eth <- gsub("syrian", "Syrian", pca$eth)
pca$eth <- gsub("Turkey_SEurope", "Turkish", pca$eth)
pca$eth <- gsub("turkish", "Turkish", pca$eth)
pca$eth <- gsub("iranian", "Iranian", pca$eth)
pca$eth <- gsub("egyptian", "Egyptian", pca$eth)

# Plot PC1 vs PC2, coloring points by cluster
colors <- c(
  "#1F77B4FF",
  "#FF7F0EFF",
  "#2CA02CFF",
  "#D62728FF",
  "#9467BDFF",
  "#8C564BFF",
  "#E377C2FF",
  "#7F7F7FFF",
  "#BCBD22FF",
  "#17BECFFF",
  "#F9CA24FF",
  "#C7ECEEFF",
  "#CD534CFF",
  "#8F7700FF",
  "#003C67FF",
  "#A20056FF"
)

# plot by nationality
p1 <- ggplot(data = pca, aes(x = PC1, y = PC2, label = nat, colour = nat, shape = nat)) +
  geom_point(size = 2, stroke = 1, alpha = 0.8, aes(color = nat)) +
  scale_color_manual(
    values = rep(colors, length.out = nlevels(factor(pca$nat)))  # Repeat colors as needed
  ) +
  scale_shape_manual(values = rep(0:20, length.out = nlevels(pca$nat))) +  # Shape mapping based on levels of 'nat'
  labs(
    x = paste("PC1", " (", round(var[1, 2], 2), "%)", sep = ""),
    y = paste("PC2", " (", round(var[2, 2], 2), "%)", sep = ""),
    colour = "Country",  
    shape = "Country"    
  ) +
  theme_classic() +
  guides(
    colour = guide_legend(ncol = 3),  # Set color legend to 3 columns
    shape = guide_legend(ncol = 3)    # Set shape legend to 3 columns
  ) +
  theme(
    legend.position = c(.76, 0.22),
    legend.key = element_rect(color = "grey", size = 0.3),  # Border around legend keys (individual legend items)
    legend.background = element_rect(color = "grey", size = 0.3)  # Border around the entire legend box
  )

p1

ggsave(paste0(output, "country_targetPCA.png"), p1, width = 7, height = 7)

# plot by study
p2 <- ggplot(data = pca, aes(x = PC1, y = PC2, label = study, colour = study, shape = study)) +
  geom_point(size = 2, stroke = 1, alpha = 0.8, aes(color = study)) +
  scale_color_manual(
    values = rep(colors, length.out = nlevels(factor(pca$study)))  # Repeat colors as needed
  ) +
  scale_shape_manual(values = rep(1:25, length.out = nlevels(factor(pca$study)))) +  # Repeat shapes from 1 to 25
  labs(
    x = paste("PC1", " (", round(var[1, 2], 2), "%)", sep = ""),
    y = paste("PC2", " (", round(var[2, 2], 2), "%)", sep = ""),
    colour = "Study",  
    shape = "Study"    
  ) +
  theme_classic() +
  guides(
    colour = guide_legend(ncol = 3),  # Set color legend to 3 columns
    shape = guide_legend(ncol = 3)    # Set shape legend to 3 columns
  ) +
  theme(
    legend.position = c(.76, .25),
    legend.key = element_rect(color = "grey", size = 0.3),  # Border around legend keys (individual legend items)
    legend.background = element_rect(color = "grey", size = 0.3)  # Border around the entire legend box
  )


p2
ggsave(paste0(output, "study_targetPCA.png"), p2, width = 7)

# plot by ethnicity
p3 <- ggplot(data = pca, aes(x = PC1, y = PC2, label = eth, colour = eth, shape = eth)) +
  geom_point(size = 2, stroke = 1, alpha = 0.8, aes(color = eth)) +
  scale_color_manual(
    values = rep(colors, length.out = nlevels(factor(pca$eth)))  # Repeat colors as needed
  ) +
  scale_shape_manual(values = rep(1:25, length.out = nlevels(factor(pca$eth)))) +  # Repeat shapes from 1 to 25
  labs(
    x = paste("PC1", " (", round(var[1, 2], 2), "%)", sep = ""),
    y = paste("PC2", " (", round(var[2, 2], 2), "%)", sep = ""),
    colour = "Ethnicity/nationality",  
    shape = "Ethnicity/nationality"    
  ) +
  theme_classic() +
  theme(
    legend.position = "right",  # Position the legend to the right
    legend.key = element_rect(color = "grey", size = 0.3),  # Border around legend keys (individual legend items)
    legend.background = element_rect(color = "grey", size = 0.3)  # Border around the entire legend box
  )

p3
ggsave(paste0(output, "ethnicity_targetPCA.png"), p3, width = 13)


## combine the 3 plots
final_plot <- plot_grid(
  plot_grid(p1, p2, labels = c("(a)", "(b)"), label_size = 16, ncol = 2, align = "hv"),
  plot_grid(p3, labels = "(c)", label_size = 16),
  ncol = 1,
  rel_heights = c(1, 1)
)

final_plot
ggsave(paste0(output, "figsup1.png"), final_plot, width = 14, height = 16)

