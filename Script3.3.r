# Script3.3 for ranking enhancer peaks
# 20251112

library("tidyverse")
library("dplyr")


dir <- "/Users/k2477939/Library/CloudStorage/OneDrive-King'sCollegeLondon/General - Laboratory of Molecular Immunotherapy and Antibiotics/Data/Georgia Miller/2025/CUT&Tag_analysis/"

# Classify ATAC-seq peaks into promoters vs enhancers using CUT&Tag consensus peaks -------------------------

# these files have ATAC-seq SteE-dependent accessible regions classifified into promoters or enhancers using CUT&Tag consensus peaks
promoters <- read.table(paste0(dir, "with_3_reps/SteE_vs_WT_peaks_promoters.bed"), header = TRUE) 
  mutate(PeakID = paste0(chrom, "_", start, "_", end)) %>% 
  pull(PeakID)

enhancers <- read.table(paste0(dir, "with_3_reps/SteE_vs_WT_peaks_enhancers.bed"), header = TRUE) 
  mutate(PeakID = paste0(chrom, "_", start, "_", end)) %>% 
  pull(PeakID)

# load ATAC file of top 300 regions
atac <- read_csv(paste0(dir, "peak_data/ATAC_seq_top300_steE_vs_WT_peaks_iBMDMs.csv")) %>% 
  separate_wider_delim(cols = PeakID, delim = "_", names = c("chrom", "start", "end"))

atac_chr2 <- atac %>%  filter(chrom == "chr2") %>% 
  arrange(start)

# add the t-statistic form the limma analysis on ATAC peaks to the enhancer and promoter lists
intersect(atac, promoters)

promotersX <- atac %>% 
  filter(PeakID == "chr1_36569922_36570112")

head(promotersX)
