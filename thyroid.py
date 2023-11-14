import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import json


# Set genes of interest
thyr_json = json.load(open('pathways/KEGG_THYROID_CANCER.v2023.1.Hs.json'))
thyr_genes = []

for gene in thyr_json['KEGG_THYROID_CANCER']['geneSymbols']:
    thyr_genes.append(gene)
    
thyr_genes = ['NTRK1', 'RET', 'TPM3', 'NCOA4', 'BRAF', 'NRAS', 'TFG',
              'CCND1', 'PAX8', 'MYC', 'TPR', 'CCDC6', 'CTNNB1', 'KRAS',
              'PPARG', 'HRAS']


# Enrichment Score Analysis (Up or Downregulated)
counts_df_sorted = counts_df.sort_values('log2FC', ascending=False)
barcode = counts_df_sorted.index.isin(thyr_genes)

fig = plt.figure(figsize=(5,5))
plt.scatter(counts_df['log2FC'], -np.log10(counts_df['ttest_p']))
counts_df_slice = counts_df[counts_df.index.isin(thyr_genes)]
plt.scatter(counts_df_slice['log2FC'], -np.log10(counts_df_slice['ttest_p']), label='Thyroid Cancer Genes')

ax = plt.gca()
ax.set_xlim([-6,6])
ax.set_ylim([0,8])
ax.legend(fontsize=16)
ax.axvline(0, color='k', linewidth=0.5)
ax.tick_params(labelsize=14) # set the fontsize of the tick labels
ax.set_xlabel(r'log$_2$(FC)', fontsize=16)
ax.set_ylabel(r'-log$_{10}$(pval)', fontsize=16)

# enrichment score
running_sum_storage = []
running_sum = 0
running_sum_storage.append(running_sum)
n_in_pathway = np.sum(barcode)
n_not_in_pathway = len(barcode) - np.sum(barcode)
for ele in barcode:
    if ele: # gene is in gene set
        running_sum += 1/n_in_pathway
    else: # gene is not in gene set
        running_sum -= 1/n_not_in_pathway
    running_sum_storage.append(running_sum)
running_sum_storage = np.array(running_sum_storage)

es_true = running_sum_storage[np.argmax(np.abs(running_sum_storage))]
print("Since the enrichment score", es_true, "is positive, we are studying the upregulated genes from the general thyroid cancer pathway.")

# Distribution Plot
n_permutations = 1000
es_null_storage_array = []
for permute_n in range(n_permutations):
    new_labels = np.random.permutation(labels)

    control_vals = counts_df[sample_names[new_labels == 'control']].values
    treatment_vals = counts_df[sample_names[new_labels == 'carc']].values

    # calculate log2 fold change
    counts_df_4 = counts_df.copy()
    counts_df_4['log2FC'] = np.log2(np.mean(treatment_vals, axis=1)/np.mean(control_vals, axis=1))
    counts_df_4 = counts_df_4.replace([np.inf, -np.inf], np.nan)
    counts_df_4 = counts_df_4[counts_df_4["log2FC"].notnull()]
    counts_df_4_sorted = counts_df_4.sort_values('log2FC', ascending=False)

    barcode = counts_df_4_sorted.index.isin(thyr_genes)

    running_sum_storage = []
    running_sum = 0
    running_sum_storage.append(running_sum)
    n_in_pathway = np.sum(barcode)
    n_not_in_pathway = len(barcode) - np.sum(barcode)
    for ele in barcode:
        if ele: # gene is in gene set
            running_sum += 1/n_in_pathway
        else: # gene is not in gene set
            running_sum -= 1/n_not_in_pathway
        running_sum_storage.append(running_sum)
    running_sum_storage = np.array(running_sum_storage)
    # get value of max deviation from zero
    es_null_storage_array.append(running_sum_storage[np.argmax(np.abs(running_sum_storage))])

es_null_storage_array = np.array(es_null_storage_array)

plt.figure(figsize=(8,5))
sns.distplot(es_null_storage_array, bins=np.arange(-1,1.01,0.1))

ax = plt.gca()
ax.axvline(es_true, color='r', label='true ES')
ax.legend(fontsize=14)
ax.set_xlabel('ES (null)', fontsize=16)
ax.set_ylabel('density', fontsize=16)
ax.tick_params(labelsize=14)

pos_es = es_null_storage_array[es_null_storage_array > 0] # get array of positive enrichment scores
print("The empirical p-value of", np.sum(pos_es >= es_true)/len(pos_es), "is above the 0.05 threshold, meaning that the set of general thyroid cancer genes are not all significantly upregulated in post-Chernobyl radiation-induced pediatric thyroid cancers.")