import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import json

counts_df = pd.read_csv('chernobyl_raw_counts.tsv', sep = '\t')
counts_df = counts_df.set_index('Gene Name')
counts_df = counts_df[counts_df.index.notnull()]
counts_df.rename(columns = {'SRR934360':'control_SRR934360', 'SRR934361':'control_SRR934361', 'SRR934362':'control_SRR934362', 'SRR934363':'control_SRR934363', 'SRR934364':'control_SRR934364', 'SRR934365':'carc_SRR934365', 'SRR934366':'carc_SRR934366', 'SRR934367':'carc_SRR934367', 'SRR934368':'carc_SRR934368', 'SRR934369':'carc_SRR934369', 'SRR934370':'carc_SRR934370'}, inplace = True)
counts_df = counts_df.drop(columns = ['Gene ID', 'control_SRR934360', 'carc_SRR934365', 'carc_SRR934368'])
sample_names = counts_df.columns

# label our data
labels = np.array(['control']*4 + ['carc']*4)
# take our control and carcinoma data
control_vals = counts_df[sample_names[labels == 'control']].values
treatment_vals = counts_df[sample_names[labels == 'carc']].values
# perform a t-test
counts_df['ttest_p'] = scipy.stats.ttest_ind(control_vals.T,treatment_vals.T, equal_var=False).pvalue
# calculate log2 fold change
counts_df['log2FC'] = np.log2(np.mean(treatment_vals, axis=1)/np.mean(control_vals, axis=1))

counts_df = counts_df.replace([np.inf, -np.inf], np.nan)
counts_df = counts_df[counts_df["log2FC"].notnull()]
counts_df