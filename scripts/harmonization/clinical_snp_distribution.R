# ==========================================================

# Author: Hatoon AlAli
# Paper: An Expanded Reference for Middle Eastern Populations Improves GWAS
# Compares the clinical distribution of a clinpgx variant by region
# R version: 4.4.1

# ==========================================================

# Clinically relevant snp distribution
library(tidyr)
library(dplyr)

freq=read.table("~/midProject/results/raw_outputs/annotation/allcohorts_regions.frq.strat", header = TRUE)
regions=read.table("~/midProject/results/raw_outputs/pca/allcohorts_random.regions.tsv")
clinpgx=read.csv("~/midProject/data/reference/annotations/clinpgx/allcohorts_random_strictqc_clinpgx.intersect.txt", sep="\t", header = FALSE)
output="~/midProject/data/reference/annotations/clinical_snps_regions.tsv"

### count number of samples per region
### Get freq per snp
# Exclude reads with high missing rate (each pop has 200 samples, keep reads > 50 (>10%))
#freq <- freq[freq$NCHROBS > 50, ]

# Transform freq df
freq_wide <- freq %>%
  pivot_wider(
    id_cols = c(CHR, SNP, A1, A2),
    names_from = CLST,
    values_from = c(MAF, MAC, NCHROBS),
    names_glue = "{CLST}.{.value}"
  )
#exclude NA reads
freq_wide <- na.omit(freq_wide)

### Intersect snps with annotation
colnames(clinpgx) <- c('CHR', 'pos', 'ref', 'alt', 'clinpgx')
clinpgx$CHR <- as.numeric(sub("^chr", "", clinpgx$CHR))
freq_wide$pos <- sapply(strsplit(freq_wide$SNP, ":"), `[`, 2)
freq_wide$ref <- sapply(strsplit(freq_wide$SNP, ":"), `[`, 3)
freq_wide$alt <- sapply(strsplit(freq_wide$SNP, ":"), `[`, 4)
freq_wide$CHR <- as.integer(freq_wide$CHR )
freq_wide$pos <- as.integer(freq_wide$pos )
df <- left_join(df, clinpgx, by=c('CHR', 'pos', 'ref', 'alt'))
df <- df[!(is.na(df$clinpgx)), ]
df <- df[!duplicated(df), ]
df <- df %>%
  dplyr::group_by(SNP) %>%
  dplyr::summarise(dplyr::across(everything(), ~ paste(unique(.), collapse = ",")))

# save
write.table(df, output, quote = FALSE, row.names = FALSE, sep="\t", na = "NA")

# Compute AF by definition
tmp <- df[df$SNP =="chr1:155313038:A:G", ]
attach(tmp)
mid1 <- (ArabPen.MAC + NArab.MAC)/(ArabPen.NCHROBS + NArab.NCHROBS)
mid2 <- (ArabPen.MAC + NArab.MAC + NAfrica.MAC)/(ArabPen.NCHROBS + NArab.NCHROBS + NAfrica.NCHROBS)
mid3 <- (ArabPen.MAC + NArab.MAC + NAfrica.MAC + EAfrica.MAC)/(ArabPen.NCHROBS + NArab.NCHROBS + NAfrica.NCHROBS + EAfrica.NCHROBS)
mid4 <- (ArabPen.MAC + NArab.MAC + NAfrica.MAC + SWAsia.MAC)/(ArabPen.NCHROBS + NArab.NCHROBS + NAfrica.NCHROBS + SWAsia.NCHROBS)
