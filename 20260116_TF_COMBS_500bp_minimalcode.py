# TF-COMBs to find what TFs co-occur with STAT3 on the top 300 differentially accessible peaks between WT and SteE KO Salmonella infected iBMDMs

# to install TF-COMB:
#%pip install git+https://github.com/loosolab/TF-COMB

# NB: need to comment out the line in the utils.py file starting: turtle import end_fill ...
# NB: also need to replace in objects.py: var_name=['TF2'] with var_name='TF2' so market_basket works

import os

import pandas as pd

from tfcomb import CombObj, DiffCombObj
import tfcomb


##### do for all peaks (background) #####

# set up a CombObj
background = CombObj()

# search for motifs
dir_path = "/scratch/SteE_revisions/"
background.TFBS_from_motifs(regions = os.path.join(dir_path, "TF-COMBs/ATACSeq.OCRs.DESeq2.ssaV_vs_WT.all_regions.bed"),
                   motifs = os.path.join(dir_path, "TF-COMBs/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"),
                   genome = os.path.join(dir_path, "genome.fa"),
                   #motif_pvalue = 0.01, # can use to be less strict than default 0.00001
                   threads = 10)

# overwrite counts with larger window
background.count_within(max_dist=500,
                        binarize = False, # if True then each cooccurence is only counted once per window
                        threads=10)

background.market_basket(threads=10)



##### do for WT vs SPI2 opened iBMDM peaks #####

# set up a CombObj
opened = CombObj()

# search for motifs
opened.TFBS_from_motifs(regions = os.path.join(dir_path, "TF-COMBs/ATACSeq.OCRs.DESeq2.ssaV_vs_WT.SPI2_opened_regions.bed"),
                   motifs = os.path.join(dir_path, "TF-COMBs/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"),
                   genome = os.path.join(dir_path, "genome.fa"),
                   #motif_pvalue = 0.01, # # can use to be less strict than default 0.00001
                   threads = 4)

# overwrite counts with larger window
opened.count_within(max_dist=500,
                    #binarize = True, # if True, then each cooccurence is only counted once per window
                    threads=10) 

opened.market_basket(threads=10)



##### do differential analysis between background and differentially expressed peaks #####

# combine the 2 objects, positive log2fc will show TF pairs enriched in the opened regions (first object)
opened.set_prefix("opened")
background.set_prefix("background")

diff = DiffCombObj([opened, background],
                   join = 'inner' # only TFs present on both objects are kept
                       )

# normalise
diff.normalize()

# calculate fold changes
diff.calculate_foldchanges() # use default pseudo = 0.01, defaults to cosine



#### generate table for enriched stat3 pairs containing metrics used in market basket association rules ####

# get metrics from both objects (excluding cosine as it's already in diff.rules)
opened_metrics = opened.rules[['TF1_TF2_count', 'TF1_count', 'TF2_count', 'zscore']]
opened_metrics.columns = ['opened_cooccur_count', 'opened_TF1_count', 'opened_TF2_count', 'opened_zscore']

background_metrics = background.rules[['TF1_TF2_count', 'TF1_count', 'TF2_count', 'zscore']]
background_metrics.columns = ['background_cooccur_count', 'background_TF1_count', 'background_TF2_count', 'background_zscore']

# merge everything
combined = diff.rules.copy()
combined = combined.join(opened_metrics, how='left')
combined = combined.join(background_metrics, how='left')

# calculate key significance metrics
n_opened_regions = len(opened.regions) #1385
n_background_regions = len(background.regions) #27041

combined['opened_pct_regions'] = (combined['opened_cooccur_count'] / n_opened_regions) * 100
combined['background_pct_regions'] = (combined['background_cooccur_count'] / n_background_regions) * 100

# enrichment ratio (raw fold change, not log2)
combined['enrichment_ratio'] = combined['opened_cosine'] / combined['background_cosine']

# difference in z-scores (higher = more statistically different)
combined['zscore_diff'] = combined['opened_zscore'] - combined['background_zscore']

# absolute counts for context
combined['count_difference'] = combined['opened_cooccur_count'] - combined['background_cooccur_count']

# final table with most relevant columns
significance_table = combined[[
    'TF1', 'TF2',
    # opened metrics
    'opened_cooccur_count', 'opened_pct_regions', 'opened_cosine', 'opened_zscore',
    # background metrics  
    'background_cooccur_count', 'background_pct_regions', 'background_cosine', 'background_zscore',
    # combined comparison metrics
    'opened/background_cosine_log2fc',  # magnitude of enrichment
    'enrichment_ratio',                  # fold change (non-log)
    'zscore_diff',                       # statistical difference
    'count_difference'                   # absolute difference in co-occurrences
]]

# filter for STAT3 enriched pairs and sort
stat3_enriched = significance_table[
    (significance_table.index.str.contains('STAT3')) &
    (significance_table['opened/background_cosine_log2fc'] > 0)
].sort_values('opened/background_cosine_log2fc', ascending=False)

# remove duplicates
stat3_enriched_unique = stat3_enriched[stat3_enriched['TF1'] <= stat3_enriched['TF2']]

# save table
stat3_enriched_unique.to_csv(os.path.join(dir_path, "TF-COMBs/STAT3_TF_pairs_enrichedInDiffAccessOCRs_500bpwindow_allmetrics.csv"))


