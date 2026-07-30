# ==========================================================

# Author: Hatoon AlAli
# Paper: 
# Generates Figure 2 (ADMIXTURE + PCA for target nd neighboring cohorts)
# R version: 4.4.1

# ==========================================================
# 0. Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(gridExtra)
  library(grid)
  library(packcircles)
  library(cluster)
  library(patchwork)
  library(cowplot)
})

# 1. Paths & Parameters
ADMIX_IN  <- "~/midProject/results/raw_outputs/admixture/allcohorts_random"
PCA_IN    <- "~/midProject/results/raw_outputs/pca/allcohorts_random"
ADMIX_OUT <- "~/midProject/results/figures/harmonization/admixture/"
PCA_OUT   <- "~/midProject/results/figures/harmonization/pca/"

k <- 6       # chosen K base on K error plot
pcs <- 5       # number of PCs used for clustering based on scree plot

set.seed(1234)

# ==========================================================
# PART I — ADMIXTURE
# ==========================================================


# ----------------------------------------------------------
# READ INPUTS
# ----------------------------------------------------------

data <- read.delim(paste0(ADMIX_IN, "_kerror.tsv"), header = FALSE, col.names = c("K", "Cross_Validation"))
# read table with lowest cross validation error (or the selected K)
q <- read.table(paste0(ADMIX_IN, ".",k,".Q"))
id <- read.table(paste0(ADMIX_IN, ".fam"))
q$region <- sub(".*_", "", id$V1)
q$id <- paste(id$V1, "_", id$V2, sep = "")

# ----------------------------------------------------------
# STEP 0:  ADMIXTURE K‑error (Cross‑Validation)
# ----------------------------------------------------------

png(filename=paste0(ADMIX_OUT, "kerror_admixture.png"))
plot(data$K, data$Cross_Validation, type = "o", # 'o' for both points and lines
     xlab = "K", ylab = "Cross Validation",
     col = "blue", pch = 16) # Blue color and solid circle for points
dev.off()

# ----------------------------------------------------------
# STEP 1:  Order ADMIXTURE Components
# ----------------------------------------------------------

# 1.1. Find the sample with highest admixture (smallest maximum component)
min_prop <- min(apply(q[, 1:(ncol(q) - 2)], 1, max))
remaining_columns <- names(q[, 1:(ncol(q) - 2)]) #define remaining components not added to indexes list 

# 1.2. exclude the component with the largest admixture (lowest median) from ordering (background color)
background = names(which.min(sapply(q[, 1:6], function(x) sum(x < min_prop))))
remaining_columns <- setdiff(remaining_columns, background)

# 1.3. Find component 1 (comp1) with the most values greater than 0.8
comp1 <- names(which.max(sapply(q[, 1:6], function(x) sum(x > 0.8))))
remaining_columns <- setdiff(remaining_columns, comp1)
indexes <- comp1 #define order of components plotted (left to right)

# 1.4. Find the other components, order by comp1: largest component -> comp2: component mostly shared by comp1
while(length(remaining_columns) > 1) {
  # Find top component with comp1 and move forward
  logical_matrix <- q[, remaining_columns][q[[comp1]] > min_prop, ]
  comp1 <- names(colSums(logical_matrix)[order(colSums(logical_matrix), decreasing = TRUE)])[1]
  indexes <- c(indexes, comp1)
  remaining_columns <- setdiff(remaining_columns, comp1)
}
# 1.5. Add last component
indexes <- c(indexes, remaining_columns[1])

# ----------------------------------------------------------
# STEP 2:  Order Samples
# ----------------------------------------------------------

# 2.1 Annotate the largest proportion in each sample
max_one <- matrix(0, nrow = nrow(q), ncol = k,
                  dimnames = list(rownames(q), colnames(q)[1:k]))
for (i in 1:nrow(q)) {
  x <- as.numeric(q[i, 1:k])
  # Order indices from largest to smallest
  ordered_indices <- colnames(q)[1:k][order(x, decreasing = TRUE)]
  # Take the max index
  max_index <- ordered_indices[1]
  # If the max is background, use the 2nd max
  if (max_index == background) {
    max_index <- ordered_indices[2]
  }
  # Mark the chosen cluster with 1
  max_one[i, max_index] <- 1
}

max_one <- as.data.frame(max_one)

# 2.2. bin samples and order them based on their largest components 
sample_order <- data.frame()
results <- data.frame()  # Initialize an empty data frame to store results

for (i in 1:length(indexes)) {
  # Get the current index and print it
  comp1 <- indexes[i]
  #print(comp1)
  
  # Find rows where max_one[[comp1]] equals 1
  bin_indices <- which(max_one[[comp1]] == 1)
  bin <- q[bin_indices, ]
  
  # Conditionally arrange based on whether i is odd, even, or the last index
  if (i == length(indexes)) {
    # If i is the last index, arrange in ascending order
    sample_order <- rbind(sample_order, bin[order(bin[[comp1]]), ])
  } else if (i %% 2 == 1) {
    # If i is odd, arrange in descending order
    sample_order <- rbind(sample_order, bin[rev(order(bin[[comp1]])), ])
  } else {
    # If i is even, arrange in ascending order
    sample_order <- rbind(sample_order, bin[order(bin[[comp1]]), ])
  }
}

# ----------------------------------------------------------
# STEP 3:  ADMIXTURE Plot
# ----------------------------------------------------------

# Reshape data for ggplot
admixture_df <- reshape2::melt(
  sample_order,
  id.vars = c("id", "region"),
  measure.vars = names(q[, 1:(ncol(q) - 2)]),
  variable.name = "Cluster",
  value.name = "Proportion"
)

# Ensure sample order matches `sample_order$id`
admixture_df$order <- factor(admixture_df$id, levels = sample_order$id)

# Desired order of region lines under the admixture plot (change manually, can be coded by taking most common region in comp1 ..etc)
region_order <- c("SWAsia", "ECEurope", "SEurope", "NArab", 
                  "ArabPen", "NAfrica", "CAfrica", "EAfrica", "WAfrica")

# Assign regions to unique vertical positions (bring them closer together)
SPACING <- 0.2  # Smaller value for reduced distance between lines

# Compute the start and end positions of each region line
region_line_data <- sample_order %>%
  mutate(region = factor(region, levels = region_order)) %>%
  group_by(region) %>%
  summarise(
    x_start = which(admixture_df$order == first(id))[1],  # Start position
    x_end = which(admixture_df$order == last(id))[1],     # End position
    .groups = "drop"
  ) %>%
  mutate(region_y = -as.numeric(region) * SPACING)  # Maintain consistent vertical spacing

# Define colors for clusters
colors <- c("#CC79A7", 
            "#882255", 
            "#F0E442",  
            "#009E73", 
            "#E69F00",
            "#56B4E9") 


# Create the plot
p.admix <- ggplot(admixture_df, aes(x = order, y = Proportion, fill = Cluster)) +
  # Admixture bars
  geom_bar(stat = "identity", width = 1) +
  
  # Horizontal region lines under the admixture plot
  geom_segment(data = region_line_data,
               aes(x = x_start, xend = x_end, y = region_y, yend = region_y),
               color = "gray", linetype = "solid", size = 0.8, inherit.aes = FALSE) +
  
  # Smaller vertical sample intercept lines under the admixture plot
  geom_segment(data = sample_order %>%
                 mutate(region = factor(region, levels = region_order)) %>%
                 left_join(region_line_data, by = "region"),
               aes(x = id, xend = id, y = region_y - 0.05, yend = region_y + 0.05),
               color = "black", linetype = "solid", size = 0.6, inherit.aes = FALSE) +
  
  # Add region labels on the y-axis, same level as the region line
  geom_text(data = region_line_data,
            aes(x = x_start - 8, y = region_y, label = region),  # Offset label horizontally
            color = "black", size = 5, hjust = 1, vjust = 0.5, inherit.aes = FALSE) +
  
  # Annotate "k=6" on the left of the admixture plot
  annotate("text", x = 0, y = max(admixture_df$Proportion) / 2, label = paste0("k=", k), 
           size = 5, color = "black", angle = 0, hjust = 1.5) +
  
  # Adjust y-axis to ensure space for the first label
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)),
                     limits = c(min(region_line_data$region_y) - 0.5,
                                max(admixture_df$Proportion) + 0.5)) +
  
  # Customize colors and labels
  scale_fill_manual(values = colors, guide = "none") +  # Remove the legend
  
  # Minimal theme adjustments
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.text.y = element_blank(),  # Hide y-axis text
    axis.ticks.y = element_blank(),  # Hide y-axis ticks
    panel.grid = element_blank(),    # Remove gridlines
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    plot.title = element_blank(),     # Remove the title
    plot.margin = margin(1, 5, 1, 50), # Add left margin for lables
    text = element_text(size = 18),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(clip = "off")


p.admix
ggsave(paste0(ADMIX_OUT, "allcohort_admixture.png"), plot = p.admix, width = 9, height = 6)

# save a table with Major component, region, number of samples
results <- q %>%
  mutate(`Major component` = names(.)[max.col(select(., V1:V6), ties.method = "first")]) %>%
  count(`Major component`, region, name = "count") %>%
  left_join(q %>% count(region, name = "total"), by = "region") %>%
  mutate(proportion = round(count / total, 2)) %>%
  select(`Major component`, region, count, proportion)

colnames(results) <- c("Major component", "Region", "Sample count", "Proportion of total region samples (%)")

output_file <- paste0(output, "samples_per_max_comp.tsv")
write.table(results, output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ==========================================================
# PART II — PCA
# ==========================================================

# ----------------------------------------------------------
# II.A  Read & Prepare PCA
# ----------------------------------------------------------

pca <- read.csv(paste0(PCA_IN, ".eigenvec"), header = T, sep = "\t")
var <- read.csv(paste0(PCA_IN, ".eigenval"), header = F, sep = "\t")

pca$id <- as.factor(gsub("^.*]", "", pca$X.FID))
pca$id <- as.factor(gsub("^.*\\|", "", pca$id))
pca$eth <- as.factor(gsub("\\|.*", "", pca$X.FID))
pca$eth <- as.factor(gsub(".*]", "", pca$eth))
pca$cont <- as.factor(gsub("^.*_", "", pca$id))
pca$id <- as.factor(gsub("_.*", "", pca$id))
pca$study <- as.factor(gsub("].*", "", pca$X.FID))
pca$study <- as.factor(gsub("\\[", "", pca$study))

# ----------------------------------------------------------
# STEP 1:  Determine the number of included PCS (sscree plot)
# ----------------------------------------------------------

colnames(var) <- "Variance"
var$Percentage <- round(var$Variance/sum(var$Variance)*100, 2)
var$PC <- 1:nrow(var)

p.scree <- ggplot(var, aes(x = PC, y = Percentage)) +
  geom_point(size = 3, color = "blue") +
  geom_line(group = 1, color = "blue") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = unique(var$PC)) + 
  labs(title = "Scree Plot", x = "PC", y = "Variance Explained (%)") +
  theme_classic()

p.scree
ggsave(file.path(PCA_OUT, "allcohort_scree.png"), p.scree, width = 6, height = 4, dpi = 300)

# Review scree plot and PC GWAS to determine the number of PCs to include
included_pcs <- pca[, paste0("PC", 1:pcs)]

# ----------------------------------------------------------
# STEP 2:  Perform K-means Clustering
# ----------------------------------------------------------

# Step 2: Get the optimum number of clusters for k-means clustering (using 00-DetermineK.R by  Constantinescu et al. 2022)
DetermineK = function(mydata, kmax = 10, num_iter = 100, plotrepress=FALSE){
  ## 1. Kmeans
  K = lapply(1:kmax, function(k){
    kmeans(mydata, centers = k, iter.max = num_iter, nstart = 10 )
  })  
  ## 2. Within Group SumOfSquares
  WSS <- sapply(1:kmax,function(i){
    K[[i]]$tot.withinss
  })
  ## 3. Silhouette analysis
  sil = sapply(2:kmax, function(k){
    cluster::pam(mydata, k)$silinfo$avg.width
  })
  k.sil = which.max(sil) + 1
  ## 4. Silhouette analysis 2
  #sil.est = fpc::pamk(mydata, krange=1:kmax)$nc
  ## 5. Plot
  if(plotrepress == 0){
    ## PLOT LAYOUT
    mat = matrix(c(1,2,3,3), ncol = 2, byrow = FALSE)
    layout(mat)
    ## Plot WSS
    plot(1:kmax, WSS, type="b", pch = 19, frame = FALSE, 
         xlab="Number of clusters K", 
         ylab="Tot. Within Clusters SS", 
         main = "Cluster Within Group\nSum of Squares")
    abline(v = k.sil , lty =2)
    #### PLOT Silhouette & save
    png(paste0(PCA_OUT, "silhouette_analysis_allcohort.png"), width = 800, height = 600)
    plot(c(0,sil), type = "b", pch = 19, frame = FALSE, 
         main = "Silhouette Analysis", 
         xlab="Number of clusters K", 
         ylab  = "Avg Silhouette Width")
    abline(v = k.sil , lty =2)
    dev.off()
    ### Distance Matrix
    d <- dist(mydata, method = "euclidean")
    plot(cluster::pam(d, k.sil))
  }
  ## 6. Assign clusters using BestK
  best_kmeans <- kmeans(mydata, centers = k.sil, iter.max = num_iter, nstart = 10)
  mydata$Cluster <- best_kmeans$cluster  # Add cluster labels to mydata
  return(list( K = K, BestK = k.sil, WSS = WSS, Silhouette = sil, ClusteredData = mydata))
}

result <- DetermineK(included_pcs)
k <- result$BestK # Define the number of clusters
kmeans_result <- result$ClusteredData
# Add cluster assignments to the original data
pca$Cluster <- as.factor(kmeans_result$Cluster)

# Save samples with assigned cluster (FID, IID, cluster) to be used for Fst
df <- subset(pca, select = c('X.FID', 'IID', 'Cluster'))
colnames(df)[1] <- '#FID'
write.table(df, "~/midProject/results/raw_outputs/pca/allcohorts_random.clusters.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------
# STEP 3:   PCA Plot (PC1 vs PC2)
# ----------------------------------------------------------

colors <- c("#882255",  
            "#56B4E9", 
            "#E69F00",  
            "#F0E442",  
            "#CC79A7",
            "#009E73")

# Plot PC1 vs PC2, coloring points by cluster and using letters as shapes
p1 <- ggplot(pca, aes(x = PC1, y = PC2, color = Cluster, shape = cont)) +
  geom_point(size = 3) + # Points with variable shapes and colors
  labs(x = paste("PC1", " (", round(var[1, 2], 2), "%)", sep = ""),
       y = paste("PC2", " (", round(var[2, 2], 2), "%)", sep = ""),
       shape = "Geography",
       color = "Cluster") +
  scale_color_manual(values = colors) + # Cluster color palette
  scale_shape_manual(values = seq(0, nlevels(pca$cont))) + # Map shapes to levels of pca$cont
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position = c(.75, 0.2),
    legend.box = "horizontal",
    legend.title=element_text(size=10),
    legend.text=element_text(size=8)) +
  guides(
    shape = guide_legend(order = 2),  # Shape first (left)
    color = guide_legend(order = 1)   # Color second (right)
  )

p1
ggsave(file.path(PCA_OUT, "allcohort_pc1pc2.png"), p1, width = 6, height = 6)

# ----------------------------------------------------------
# STEP 4:  Honeycomb Plot
# ----------------------------------------------------------

# Create a new dataframe with Cluster, cont, and value columns
data <- pca %>%
  select(Cluster, cont) %>%      # Select the Cluster and cont columns
  dplyr::mutate(value = 1,              # Add a value column with 1 for each row
                id = row_number())      # Add an id column with row numbers

# Create an empty list to store the results for each category
dat.gg.list <- list()

# Define the desired order of categories
cont_order <- c("ArabPen", "NArab", "NAfrica", 
                "EAfrica", "WAfrica", "CAfrica", 
                "ECEurope", "SEurope", "SWAsia")

# Convert the 'cont' column to a factor with the specified order
pca$cont <- factor(pca$cont, levels = cont_order)

# Iterate through each unique category in data$cont
for (category in unique(pca$cont)) {
  # Subset the data for the current category
  subset_data <- pca[pca$cont == category, ] %>%
    select(Cluster, cont) %>%      # Select the Cluster and cont columns
    arrange(Cluster) %>%
    dplyr::mutate(value = 1,              # Add a value column with 1 for each row
                  id = row_number())  #Add row number to order the dots
  
  # Perform the circle packing layout
  packing <- circleProgressiveLayout(subset_data$value, sizetype = "area")
  packing$theta <- atan2(packing$y, packing$x)
  packing <- packing[order(packing$theta), ]
  
  # Generate the layout vertices
  dat.gg <- circleLayoutVertices(packing, idcol = NULL, npoints = 1)
  dat.gg <- merge(dat.gg, subset_data[, c("Cluster", "cont", "id")], by.x = "id", by.y = "id", all.x = TRUE)
  
  # Store the result in the list
  dat.gg.list[[as.character(category)]] <- dat.gg
}

# Optionally, you can name the list by the categories
names(dat.gg.list) <- unique(pca$cont)

# Combine all the data frames in dat.gg.list into a single data frame
combined_dat.gg <- bind_rows(lapply(names(dat.gg.list), function(category) {
  data <- dat.gg.list[[category]]
  pca$cont <- category  # Add the cont category as a new column
  return(data)
}))

# Plot all categories in the same panel using facet_wrap
p2 <- ggplot(combined_dat.gg) + 
  geom_point(aes(x, y, colour = factor(Cluster))) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    text = element_text(size = 16)
  ) +
  scale_color_manual(values = colors) +  # Cluster color palette
  facet_wrap(~ cont)  # Create a facet for each cont category

p2
ggsave(file.path(PCA_OUT, "allcohort_honeycomb.png"), p2, width = 4, height = 4)

# ----------------------------------------------------------
# STEP 5:  Plot other PCs
# ----------------------------------------------------------

# Define principal components to plot
pc_pairs <- list(c("PC1", "PC2"), c("PC3", "PC4"), c("PC5", "PC6"))

# Generate plots for each pair
plots <- lapply(pc_pairs, function(pair) {
  ggplot(pca, aes_string(x = pair[1], y = pair[2], color = "Cluster", shape = "cont")) +
    geom_point(size = 3) +
    labs(x = paste0(pair[1], " (", round(var[as.numeric(substr(pair[1], 3, 3)), 2], 2), "%)"),
         y = paste0(pair[2], " (", round(var[as.numeric(substr(pair[2], 3, 3)), 2], 2), "%)"),
         color = "Cluster",
         shape = "Geography") +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = seq(0, nlevels(pca$cont))) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          legend.position = "right",   # Move legend to the right
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.box = "vertical")     # Make legend vertical
})

# Combine all plots into a grid with a shared vertical legend
final_plot <- wrap_plots(plots) + plot_layout(guides = "collect") & theme(legend.position = "right")
final_plot
# Save the final plot
ggsave(paste0(PCA_OUT, "allcohort_pcs1to6.png"), final_plot, width = 15, height = 5)

# ==========================================================
# PART III — FIGURE 2
# ==========================================================

fig2 <- plot_grid(
  plot_grid(p1, p2, labels = c("(a)", "(b)"), label_size = 16, ncol = 2, align = "hv", rel_widths = c(1.3, 1), label_x = -0.035),
  plot_grid(p.admix, labels = "(c)", label_size = 16),
  ncol = 1,
  rel_heights = c(1.1, 1)
)

ggsave(file.path(PCA_OUT, "fig2.png"), fig2, width = 9.5, height = 8)

# ============================= END =============================

