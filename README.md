# MIDreference

Reference code for the manuscript "An Expanded Reference for Middle Eastern Populations Improves GWAS"
prediction.

## Structure

```
harmonization/
├── harmonization.Rmd              Description of code used to generate data for figures 1 and 2: 
│                                  Harmonize multiple MID cohorts into one dataset:
│                                  liftover, QC, merge, downsampling, PCA, ADMIXTURE, Fst,
│                                  and clinically relevant SNP annotation
├── clinical_snps_regions.R        Intersect region-level allele
│                                  frequencies with ClinPGx/GWAS Catalog annotations
├── Fig01_target_pca_map.R         Code used to generate Figure 1
├── Fig02_neighboring_pca_admix.R  Code used to generate Figure 2
├── pca_country_ethnicity_study.R  Code used to generate Figure S3
└── plot_fst_heatMap.R             Code used to generate Figure S7

ancestry_inference/
├── reference_qc.Rmd               Description of code used to build the 
│                                  expanded reference panel (HGDP-1kGP-Almarri):
│                                  per-cohort QC, merge, post-merge QC, liftover to hg37
├── plink_ukbb/                    Code used for genetic ancestry inference and GWAS in UK Biobank:
│                                  Intersect UKBB with the reference panel, train a
│                                  random forest ancestry classifier, predict and refine
│                                  MID ancestry calls, and validate against place-of-birth
│                                  and GWAS
└── hail_aou/                      Code used for genetic ancestry inference and GWAS in All of us    

visualization/                     
├── HoneycombPCA.R                 PCA clustering and Honeycomb plot function
└── LabelFreeADMIX.R               ADMIXTURE ordering and plotting function

```

## Citation
