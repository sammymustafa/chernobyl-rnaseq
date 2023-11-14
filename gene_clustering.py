import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import json



def get_tmm(counts, ref_col, test_col, trim_m=0.30, trim_a=0.05):

    '''
        Calculates TMM value between ref_col and test_col
        
        counts = count matrix (pd.DataFrame)
        ref_col = reference column (str)
        test_col = test columns (str)
        trim_m = extremes to trim when taking trimmed mean of M values (float)
        trim_a = extremes to trim when taking trimmed mean of M values (float)
        
    '''

    if ref_col == test_col:
        counts_slim = counts[(counts[ref_col].values > 0) & (counts[test_col].values > 0)][[ref_col]]
    else:
        counts_slim = counts[(counts[ref_col].values > 0) & (counts[test_col].values > 0)][[ref_col, test_col]]

    n_k = counts_slim[test_col].sum()
    n_r = counts_slim[ref_col].sum()

    m_k = np.log2(counts_slim[test_col].values/n_k)-np.log2(counts_slim[ref_col].values/n_r)
    a_k = 0.5*np.log2((counts_slim[test_col].values/n_k)*(counts_slim[ref_col].values/n_r))
    w_k = (n_k - counts_slim[test_col].values)/(n_k*counts_slim[test_col].values) + (n_r - counts_slim[ref_col].values)/(n_r*counts_slim[ref_col].values)

    trim_array_m = (m_k <= np.percentile(m_k, 100*(1-(trim_m/2)))) & ((m_k >= np.percentile(m_k, 100*((trim_m/2)))))
    trim_array_a = (a_k <= np.percentile(a_k, 100*(1-(trim_a/2)))) & ((a_k >= np.percentile(a_k, 100*((trim_a/2)))))

    m_k = m_k[trim_array_m & trim_array_a]
    a_k = a_k[trim_array_m & trim_array_a]
    w_k = w_k[trim_array_m & trim_array_a]
    tmm_k = 2**(np.sum(w_k*m_k)/np.sum(w_k))
    return tmm_k

def get_m_a(counts, ref_col, test_col, trim_m=0.30, trim_a=0.05):

    '''
        Calculates TMM value between ref_col and test_col
        
        counts = count matrix (pd.DataFrame)
        ref_col = reference column (str)
        test_col = test columns (str)
        trim_m = extremes to trim when taking trimmed mean of M values (float)
        trim_a = extremes to trim when taking trimmed mean of M values (float)
        
    '''

    if ref_col == test_col:
        counts_slim = counts[(counts[ref_col].values > 0) & (counts[test_col].values > 0)][[ref_col]].copy()
    else:
        counts_slim = counts[(counts[ref_col].values > 0) & (counts[test_col].values > 0)][[ref_col, test_col]]

    n_k = 1 #counts_slim[test_col].sum()
    n_r = 1 #counts_slim[ref_col].sum()

    m_k = np.log2(counts_slim[test_col].values/n_k)-np.log2(counts_slim[ref_col].values/n_r)
    a_k = 0.5*np.log2((counts_slim[test_col].values/n_k)*(counts_slim[ref_col].values/n_r))
    return m_k, a_k

def norm_tmm(counts, columns_to_norm, ref_col=None):
    '''
        Normalizes count matrix by TMM.
        
        counts = count matrix (pd.DataFrame)
        ref_col = reference column (str)
        columns_to_norm = columns to be normed (str)
        
    '''
    if ref_col == None:
        ref_col = columns_to_norm[0]
    counts[counts.columns[~counts.columns.isin(cols_to_norm)]]
        
    tmm_array = []
    for col in columns_to_norm:
        tmm_array.append(get_tmm(counts = counts, ref_col=ref_col, test_col=col))
    tmm_array = np.array(tmm_array)
    #print(tmm_array)

    norm_counts = counts[columns_to_norm]/(np.sqrt(tmm_array)*counts[columns_to_norm].sum(axis=0))
    norm_counts[counts.columns[~counts.columns.isin(columns_to_norm)]] = counts[counts.columns[~counts.columns.isin(columns_to_norm)]]
    return norm_counts

counts = pd.read_csv('chernobyl_raw_counts.tsv', sep='\t')
counts = counts[counts["Gene Name"].notnull()]
design = pd.read_csv('experiment_design.tsv', sep='\t')
design = design[design['Analysed'] == 'Yes'].reset_index(drop=True)

design['labels'] = ['']*len(design)

design['labels'].loc[(design['Factor Value[disease]'].values == "papillary thyroid carcinoma")] = 'PTC'
design['labels'].loc[(design['Factor Value[disease]'].values == "normal")] = 'NORM'

design = design.sort_values('labels')
cols_to_norm = design['Run'].values
labels = design['labels'].values

counts_norm = counts.copy()

counts_norm = norm_tmm(counts_norm, columns_to_norm=cols_to_norm, ref_col=cols_to_norm[0])

a = counts_norm[cols_to_norm].values[:,labels == 'PTC']
b = counts_norm[cols_to_norm].values[:,labels == 'NORM']
counts_norm['ttest_p'] = scipy.stats.f_oneway(a.T, b.T).pvalue
counts_norm = counts_norm.dropna(subset=['ttest_p'])
counts_norm['bonferroni_fwer'] = counts_norm['ttest_p']*len(counts_norm['ttest_p'])
counts_norm['bonferroni_fwer'][counts_norm['bonferroni_fwer'] > 1] = 1.0
counts_norm['bh_fdr'] = counts_norm['ttest_p']*len(counts_norm['ttest_p'])/counts_norm['ttest_p'].rank()
counts_norm['bh_fdr'][counts_norm['bh_fdr'] > 1] = 1.0

counts_norm_slim = counts_norm[counts_norm['bonferroni_fwer'] < 0.5]
counts_norm_slim = counts_norm_slim.set_index('Gene Name')[cols_to_norm]
counts_norm_slim.columns = labels

# Define colors to label our columns to make the diagnosis categories of our clustermap easier to read
lut = dict(zip(set(labels), ['r', 'b'])) # red, blue
col_colors = pd.DataFrame(labels)[0].map(lut)

g = sns.clustermap(counts_norm_slim, # data to cluster
                   method='ward', # method for heiarchical clustering
                   col_cluster=True, # cluster columns
                   z_score=0, # axis to perform z standardization
                   vmin=-2.5, vmax=2.5, # upper and lower limits of our colorbar
                   cmap='seismic_r', # palette to use for our colorbar
                   col_colors=[col_colors]) # what colors to label our columns on the margin

# define legend for col colors
ax = g.ax_heatmap
handles1 = [mpl.patches.Patch(facecolor=lut[name]) for name in lut]

leg1 = g.ax_row_dendrogram.legend(handles1, lut, title='diagnosis',
           bbox_to_anchor=(0.97, 0.70), bbox_transform=plt.gcf().transFigure, loc='upper left', frameon=False, fontsize=16)
plt.setp(leg1.get_title(),fontsize=16)
leg1._legend_box.align = "left"