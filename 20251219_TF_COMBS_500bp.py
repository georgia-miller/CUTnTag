#!/usr/bin/env python
# coding: utf-8

# In[1]:


# 18/12/2025 Running TF-COMBs to find what TFs co-occur with STAT3 on the top 300 differentially accessible peaks between WT and SteE KO Salmonella infected iBMDMs

# Using new tables Peter sent me for iBMDMs
#%pip install --upgrade pip
#%pip install git+https://github.com/loosolab/TF-COMB

# NB: need to comment out this line from turtle import end_fill from the utils.py file
# NB: also need to replace var_name=['TF2'] with var_name='TF2' in objects.py so market_basket works

# needed for os.path.join
import os

# use to show all outputs from a cell, not just the last one
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

import pandas as pd

from tfcomb import CombObj, DiffCombObj

import tfcomb
print("tfcomb version: ", tfcomb.__version__)


# In[2]:


# set up a CombObj
background = CombObj()

# search for motifs
dir_path = "/scratch/prj/id_hill_sims_wellcda/SteE_revisions/"
background.TFBS_from_motifs(regions = os.path.join(dir_path, "TF-COMBs/ATACSeq.OCRs.DESeq2.ssaV_vs_WT.all_regions.bed"),
                   motifs = os.path.join(dir_path, "TF-COMBs/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"),
                   genome = os.path.join(dir_path, "genome.fa"),
                   #motif_pvalue = 0.01, # less strict than default 0.00001
                   threads = 10)

# overwrite counts with larger window
background.count_within(max_dist=500,
                        #binarize = True, # if True then each cooccurence is only counted once per window
                        threads=10)

background.market_basket(threads=10)


# In[3]:


background.rules


# In[4]:


_ = background.plot_heatmap()
_ = background.plot_bubble()


# In[5]:


# filter for STAT3 containing pairs
background.TFBS[0].__dict__

# look at all TFBS sites containing STAT3
# this way returns a list of TFBS object
stat3_tfbs = [x for x in background.TFBS if "STAT3" in x.name]
stat3_tfbs[:5]

len(stat3_tfbs)
## 1,773 stat3 sites

# get all STAT3 motif variants
stat3_variants = sorted(set(x.name for x in background.TFBS if "STAT3" in x.name))
print(stat3_variants)

# filter the significant co-occurrence pairs to those with STAT3 co-occurrence
sorted(set(x.name for x in background.TFBS if "STAT3" in x.name))

stat3_pairs = background.select_TF_rules(["STAT3_MOUSE.H11MO.0.A"])

stat3_pairs.rules.head()


# In[6]:


background.__dict__.keys()


# In[7]:


_ = stat3_pairs.plot_heatmap()
_ = stat3_pairs.plot_bubble()


# In[8]:


# do it for WT vs SPI2 opened iBMDM peaks
# set up a CombObj
opened = CombObj()

# search for motifs
dir_path = "/scratch/prj/id_hill_sims_wellcda/SteE_revisions/"
opened.TFBS_from_motifs(regions = os.path.join(dir_path, "TF-COMBs/ATACSeq.OCRs.DESeq2.ssaV_vs_WT.SPI2_opened_regions.bed"),
                   motifs = os.path.join(dir_path, "TF-COMBs/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"),
                   genome = os.path.join(dir_path, "genome.fa"),
                   #motif_pvalue = 0.01, # default is 0.00001
                   threads = 4)

# overwrite counts with larger window
opened.count_within(max_dist=500,
                    #binarize = True, # if True, then each cooccurence is only counted once per window
                    threads=10) 

opened.market_basket(threads=10)


opened.rules
_ = opened.plot_heatmap()
_ = opened.plot_bubble()


# In[9]:


# look at all TFBS sites containing STAT3
# this way returns a list of TFBS object
stat3_tfbs1 = [x for x in opened.TFBS if "STAT3" in x.name]
stat3_tfbs1[:5]

len(stat3_tfbs1)
## 47 stat3 sites

# get all STAT3 motif variants
stat3_variants1 = sorted(set(x.name for x in opened.TFBS if "STAT3" in x.name))
print(stat3_variants1)

# filter the significant co-occurrence pairs to those with STAT3 co-occurrence
sorted(set(x.name for x in opened.TFBS if "STAT3" in x.name))

stat3_pairs1 = opened.select_TF_rules(["STAT3_MOUSE.H11MO.0.A"])

stat3_pairs1.rules.head()


# In[10]:


# do differential analysis between background and differentially expressed peaks

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

print(diff.rules.columns)

# print table
print("All comparisons:")
pd.set_option('display.max_columns', None)  # show all columns
pd.set_option('display.width', None)        # auto-detect width
pd.set_option('display.max_colwidth', None) # don't truncate column content
diff.rules.sort_values('opened/background_cosine_log2fc', ascending=False)

print("Enriched in opened regions:")
diff.rules[diff.rules['opened/background_cosine_log2fc'] > 0
    ].sort_values('opened/background_cosine_log2fc', ascending = False
                 ).head(15)


print("STAT3 is the 13th most enriched TF pair out of all TF pairs")

# plot
diff.plot_heatmap()


# In[11]:


# filter for stat3 pairs enriched in opened regions
diff_stat3 = diff.select_TF_rules(["STAT3_MOUSE.H11MO.0.A"])
print(diff_stat3.rules.columns)

# inspect positive log2fc (enriched in opened vs all)
# print table
print("Enriched in opened regions:")
pd.set_option('display.max_columns', None)  # show all columns
pd.set_option('display.width', None)        # auto-detect width
pd.set_option('display.max_colwidth', None) # don't truncate column content
diff_stat3.rules.sort_values('opened/background_cosine_log2fc', ascending=False)

# plot
diff_stat3.plot_heatmap()


# In[12]:


# realised that if filter in a different way, get non-sf tf pairs as well
# select rules only keeps signficant ones - check this is why others aren't showing
diff_stat3_all = diff.rules[
    diff.rules.index.str.contains('STAT3')
].sort_values('opened/background_cosine_log2fc', ascending=False)

print(f"Total STAT3 rows: {len(diff_stat3_all)}")

# Remove duplicates (keep only one direction per pair)
diff_stat3_unique = diff_stat3_all[diff_stat3_all['TF1'] <= diff_stat3_all['TF2']]
print(f"Unique STAT3 pairs: {len(diff_stat3_unique)}")

# Look at top pairs
print("\nTop 10 enriched STAT3 pairs:")
diff_stat3_unique.head(10)

# Check how many are enriched vs depleted
enriched = diff_stat3_unique[diff_stat3_unique['opened/background_cosine_log2fc'] > 0]
depleted = diff_stat3_unique[diff_stat3_unique['opened/background_cosine_log2fc'] < 0]
print(f"\nEnriched in opened: {len(enriched)}")
print(f"Depleted in opened: {len(depleted)}")

# save csv
#enriched.to_csv(os.path.join(dir_path, "TF-COMBs/STAT3_TF_pairs_enrichedInDiffAccessOCRs_500bpwindow.csv"))


# In[13]:


# look at structure of the objects
print("Opened object rules shape:", opened.rules.shape)
print("\nColumns in opened.rules:")
print(opened.rules.columns)

print("Diff object rules shape:", diff.rules.shape)
print("\nColumns in diff.rules:")
print(diff.rules.columns)


# In[17]:


# generate table for enriched stat3 pairs that has the association rule metrics in it

# get metrics from both objects (excluding cosine since it's already in diff.rules)
opened_metrics = opened.rules[['TF1_TF2_count', 'TF1_count', 'TF2_count', 'zscore']]
opened_metrics.columns = ['opened_cooccur_count', 'opened_TF1_count', 'opened_TF2_count', 'opened_zscore']

background_metrics = background.rules[['TF1_TF2_count', 'TF1_count', 'TF2_count', 'zscore']]
background_metrics.columns = ['background_cooccur_count', 'background_TF1_count', 'background_TF2_count', 'background_zscore']

# merge everything
combined = diff.rules.copy()
combined = combined.join(opened_metrics, how='left')
combined = combined.join(background_metrics, how='left')

# calculate key significance metrics
n_opened_regions = 1385
n_background_regions = 27041

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

print(f"STAT3 enriched pairs: {len(stat3_enriched_unique)}")
stat3_enriched_unique.head(20)

# save
#stat3_enriched_unique.to_csv(os.path.join(dir_path, "TF-COMBs/STAT3_TF_pairs_enrichedInDiffAccessOCRs_500bpwindow_allmetrics.csv"))


# In[18]:


# get the metrics used by market basket to determine what is an association or not (only the ones that pass the thresholding are kept)
# want:
    # motif pvalue = statistical significance/motif match
    # support = frequency of tf pair. count(A&B)/total number of windows i.e. baskets
    # confidence P(TF2|TF1). count(A&B)/count(A)
    # lift = enrihcment over random. support(A&B)/ [support(A)*support(B)]

# inside the opened or background object there is:
    # rules = filtered association table with TF1, TF2, support, confidence, lift, count
    # pair_counts = square matrix with TF1, TF2, values = number of regions where they co-occur
    # TFBS = list of individual TFBSs

import pandas as pd
import numpy as np

STAT3 = "STAT3_MOUSE.H11MO.0.A"

# core data
tf_names = opened.count_names
pair_counts = opened.pair_counts      # numpy array
tf_counts = opened.TF_counts           # numpy array
N = max(tf_counts)                     # total windows

i = tf_names.index(STAT3)

rows = []
for j, tf2 in enumerate(tf_names):
    c_xy = pair_counts[i, j]
    c_x = tf_counts[i]
    c_y = tf_counts[j]

    support = c_xy / N if N > 0 else 0
    confidence = c_xy / c_x if c_x > 0 else 0
    lift = (c_xy / N) / ((c_x / N) * (c_y / N)) if c_x > 0 and c_y > 0 else 0

    rows.append((STAT3, tf2, c_xy, support, confidence, lift))

stat3_metrics = pd.DataFrame(
    rows,
    columns=["TF1", "TF2", "cooccurrence_count", "support", "confidence", "lift"]
)

# remove STAT3–STAT3
stat3_metrics = stat3_metrics[stat3_metrics["TF2"] != STAT3]

# sort by lift (or confidence)
stat3_metrics.sort_values(
    ["lift", "confidence", "support"],
    ascending=[False, False, False]
).head(20)

# STAT3-STAT3 homotypic interactions were probably excluded by default from this


# In[71]:


import matplotlib.pyplot as plt
import numpy as np

# Get all STAT3 pairs
stat3_opened = opened.rules[opened.rules.index.str.contains('STAT3')]
stat3_background = background.rules[background.rules.index.str.contains('STAT3')]

# Remove duplicates
stat3_opened_unique = stat3_opened[stat3_opened['TF1'] <= stat3_opened['TF2']]
stat3_background_unique = stat3_background[stat3_background['TF1'] <= stat3_background['TF2']]

# Identify significant pairs
opened_sig = (stat3_opened_unique['zscore'] > 2) & (stat3_opened_unique['cosine'] > 0.1)
background_sig = (stat3_background_unique['zscore'] > 2) & (stat3_background_unique['cosine'] > 0.1)

# Shared axis limits
all_cosine = np.concatenate([
    stat3_opened_unique['cosine'],
    stat3_background_unique['cosine']
])
all_zscore = np.concatenate([
    stat3_opened_unique['zscore'],
    stat3_background_unique['zscore']
])

xlim = (all_cosine.min() - 0.05, all_cosine.max() + 0.05)
ylim = (all_zscore.min() - 2, all_zscore.max() + 5)

# Create figure (BACKGROUND first)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# =========================
# Plot 1: BACKGROUND
# =========================
ax_bg = axes[0]

ax_bg.scatter(
    stat3_background_unique.loc[~background_sig, 'cosine'],
    stat3_background_unique.loc[~background_sig, 'zscore'],
    alpha=0.6, s=50, c='darkgray', label='Non-significant', zorder=1
)

ax_bg.scatter(
    stat3_background_unique.loc[background_sig, 'cosine'],
    stat3_background_unique.loc[background_sig, 'zscore'],
    alpha=0.9, s=150, c='red',
    edgecolors='black', linewidths=2,
    label='Significant', zorder=2
)

ax_bg.axhline(2, linestyle='--', color='red', alpha=0.4, linewidth=1.5)
ax_bg.axvline(0.1, linestyle='--', color='red', alpha=0.4, linewidth=1.5)

for _, row in stat3_background_unique[background_sig].iterrows():
    label = f"{row['TF1'].split('_')[0]}-{row['TF2'].split('_')[0]}"
    ax_bg.annotate(
        label,
        (row['cosine'], row['zscore']),
        xytext=(8, 15),
        textcoords='offset points',
        fontsize=11,
        arrowprops=dict(
            arrowstyle='-',
            connectionstyle='arc3,rad=0',
            lw=1.5,
            color='black'
        ),
        fontweight='bold',
        zorder=3
    )

ax_bg.set_title('STAT3 Pairs – Background Regions', fontsize=15, fontweight='bold')
ax_bg.set_xlabel('Cosine Similarity', fontsize=13, fontweight='bold')
ax_bg.set_ylabel('Z-score', fontsize=13, fontweight='bold')
ax_bg.legend(loc='center right', fontsize=10)
ax_bg.grid(alpha=0.3, linestyle=':', linewidth=0.5)
ax_bg.spines['top'].set_visible(False)
ax_bg.spines['right'].set_visible(False)
ax_bg.set_xlim(xlim)
ax_bg.set_ylim(ylim)

# =========================
# Plot 2: OPENED
# =========================
ax_op = axes[1]

ax_op.scatter(
    stat3_opened_unique.loc[~opened_sig, 'cosine'],
    stat3_opened_unique.loc[~opened_sig, 'zscore'],
    alpha=0.6, s=50, c='darkgray', label='Non-significant', zorder=1
)

ax_op.scatter(
    stat3_opened_unique.loc[opened_sig, 'cosine'],
    stat3_opened_unique.loc[opened_sig, 'zscore'],
    alpha=0.9, s=150, c='red',
    edgecolors='black', linewidths=2,
    label='Significant', zorder=2
)

ax_op.axhline(2, linestyle='--', color='red', alpha=0.4, linewidth=1.5)
ax_op.axvline(0.1, linestyle='--', color='red', alpha=0.4, linewidth=1.5)

for _, row in stat3_opened_unique[opened_sig].iterrows():
    label = f"{row['TF1'].split('_')[0]}-{row['TF2'].split('_')[0]}"
    ax_op.annotate(
        label,
        (row['cosine'], row['zscore']),
        xytext=(8, 15),
        textcoords='offset points',
        fontsize=11,
        arrowprops=dict(
            arrowstyle='-',
            connectionstyle='arc3,rad=0',
            lw=1.5,
            color='black'
        ),
        fontweight='bold',
        zorder=3
    )

ax_op.set_title('STAT3 Pairs – Opened Regions', fontsize=15, fontweight='bold')
ax_op.set_xlabel('Cosine Similarity', fontsize=13, fontweight='bold')
ax_op.set_ylabel('Z-score', fontsize=13, fontweight='bold')
ax_op.legend(loc='center right', bbox_to_anchor=(1.02, 0.5), fontsize=10)
ax_op.grid(alpha=0.3, linestyle=':', linewidth=0.5)
ax_op.spines['top'].set_visible(False)
ax_op.spines['right'].set_visible(False)
ax_op.set_xlim(xlim)
ax_op.set_ylim(ylim)

plt.tight_layout()
#plt.savefig(os.path.join(dir_path, "TF-COMBs/STAT3_TF_pairs_scatterplots.pdf"), dpi=300, bbox_inches="tight") 
plt.show()

# Summary
print(f"\nSignificant STAT3 pairs in background: {background_sig.sum()}")
print(stat3_background_unique[background_sig][['TF1', 'TF2', 'cosine', 'zscore']])

print(f"\nSignificant STAT3 pairs in opened: {opened_sig.sum()}")
print(stat3_opened_unique[opened_sig][['TF1', 'TF2', 'cosine', 'zscore']])


# In[70]:


# plot pseudo-volcano plot (as zscore not pvalue)
import matplotlib.pyplot as plt

# Columns for volcano plot
x_col = "opened/background_cosine_log2fc"
y_col = "zscore_diff"

# Define thresholds
x_thresh = 0.5   # example log2FC threshold
y_thresh = 2     # example z-score difference threshold

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Plot all points
ax.scatter(
    df[x_col],
    df[y_col],
    s=40,
    c="lightgray",
    alpha=0.7,
    label="All TF pairs",
    zorder=1
)

# Highlight significant points
sig = (df[x_col].abs() >= x_thresh) & (df[y_col].abs() >= y_thresh)
ax.scatter(
    df.loc[sig, x_col],
    df.loc[sig, y_col],
    s=80,
    c="red",
    alpha=0.8,
    label="Significant",
    zorder=2
)

# Label STAT3–STAT3 with a line (arrow) only
stat3_stat3 = (
    df["TF1"].str.contains("STAT3", case=False, na=False) &
    df["TF2"].str.contains("STAT3", case=False, na=False)
)
for _, row in df.loc[stat3_stat3].iterrows():
    ax.annotate(
        "STAT3–STAT3",
        (row[x_col], row[y_col]),
        xytext=(25, 3),        # offset for text
        textcoords="offset points",
        fontsize=11,
        fontweight="bold",
        color="black",
        arrowprops=dict(arrowstyle="-", lw=1.5, color="black")
    )

# Add threshold lines
ax.axhline(y_thresh, color="black", linestyle="--", alpha=0.5)
ax.axhline(-y_thresh, color="black", linestyle="--", alpha=0.5)
ax.axvline(x_thresh, color="black", linestyle="--", alpha=0.5)
ax.axvline(-x_thresh, color="black", linestyle="--", alpha=0.5)

# Labels and styling
ax.set_xlabel("log2FC (opened / background)", fontsize=12)
ax.set_ylabel("Z-score difference", fontsize=12)
ax.set_title("STAT3-containing TF-pair enrichment in opened regions", fontsize=14, fontweight="bold")
ax.legend()
ax.grid(alpha=0.3, linestyle=":")

plt.tight_layout()
#plt.savefig(os.path.join(dir_path, "TF-COMBs/STAT3_TF_pairs_volcanoplot.pdf"), dpi=300, bbox_inches="tight") 
plt.show()


# In[ ]:




