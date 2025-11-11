library("tidyverse")
library("readr")

# convert csv file of top 300 peaks SteE vs WT STm infected iBMDMs to a bed file that can be used in bedtools multiinter

ATAC_peaks <- read_csv("/Users/k2477939/Library/CloudStorage/OneDrive-King'sCollegeLondon/General - Laboratory of Molecular Immunotherapy and Antibiotics/Data/Georgia Miller/2025/CUT&Tag_analysis/peak_data/ATAC_seq_top300_steE_vs_WT_peaks_iBMDMs.csv")
head(ATAC_peaks)

# bed file is made up of: chromosome number, start and end, must extract that from PeakID
ATAC_peaks2 <- separate_wider_delim(ATAC_peaks, cols = PeakID, delim = "_", names = c("chrom", "chromStart", "chromEnd")) %>% 
  select(chrom, chromStart, chromEnd)
head(ATAC_peaks2)

# remove chr for numeric sorting (altohugh one X chromosome) and make cols the correct types
chrom_order <- c(paste0("chr", 1:22), "chrX")

ATAC_peaks2 <- ATAC_peaks2 %>% 
  mutate(chrom = factor(chrom, levels = chrom_order),
         chromStart = as.numeric(chromStart),
         chromEnd = as.numeric(chromEnd)) %>% 
  arrange(chrom, chromStart)
head(ATAC_peaks2)


# sort by chrom then chromStart as needed for bedtools multiinter
ATAC_peaks2 <- ATAC_peaks2 %>% 
  arrange(chrom, chromStart)
head(ATAC_peaks2)

# save as a bed file
write.table(ATAC_peaks2,"/Users/k2477939/Library/CloudStorage/OneDrive-King'sCollegeLondon/General - Laboratory of Molecular Immunotherapy and Antibiotics/Data/Georgia Miller/2025/CUT&Tag_analysis/peak_data/ATAC_seq_top300_steE_vs_WT_peaks_iBMDMs.bed",
  row.names = F,col.names = F,quote = FALSE,sep="\t")
