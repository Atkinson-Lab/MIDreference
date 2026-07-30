# ==========================================================

# Author: Hatoon AlAli
# Paper: An Expanded Reference for Middle Eastern Populations Improves GWAS
# Generates Supplementary Fig. 5 (Fst heat map)
# R version: 4.4.1

# ==========================================================

# Fst
library(ggplot2)

input='/storage/atkinson/home/u250489/midProject/results/raw_outputs/fst/allcohorts_random.clusters.fst.summary'
output='~/midProject/results/figures/harmonization/fst/allcohorts_random.clusters.fst.png'

fst <- read.csv(input, sep="\t")

colors <- c("#882255",  #C1
            "#56B4E9",  #C2
            "#E69F00",  #C3
            "#F0E442",  #C4
            "#CC79A7",  #C5
            "#009E73")  #C6

#plot
fst$X.POP1 <- factor(fst$X.POP1)
fst$POP2 <- factor(fst$POP2)
fst$HUDSON_FST <- as.numeric(fst$HUDSON_FST)
fst$f <- as.numeric(format(round(fst$HUDSON_FST, 3), nsmall = 3))


cluster_colors <- setNames(colors, sort(unique(c(fst$X.POP1, fst$POP2))))
colnames(fst)[1:2] <- c("Clusters", "clusters2")

p <- ggplot(fst, aes(Clusters, clusters2)) +
  geom_tile(aes(fill = f), colour = "white") +
  geom_text(aes(label = f), size = 3) +
  scale_fill_gradient("Fst", low = "#e9e6e6", high = "#908484") +
  coord_equal(expand = FALSE) +
  labs(y = NULL) +  # remove Y-axis label
  theme(
    axis.text.x = element_text(color = cluster_colors[1:5], face = "bold"),
    axis.text.y = element_text(color = cluster_colors[2:6], face = "bold")
  )

p
# save
ggsave(output, plot = p, width = 8, height = 6, dpi = 300)
