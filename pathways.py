import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import json


fold_df = pd.read_csv('chernobyl_analytics.tsv', sep='\t')
fold_df = fold_df[fold_df["Gene Name"].notnull()]
#fold_df = fold_df[fold_df["'papillary thyroid carcinoma' vs 'normal'.p-value"].notnull()]
fold_df = fold_df[fold_df["'papillary thyroid carcinoma' vs 'normal'.log2foldchange"].notnull()]
fold_df


fold_df_slice = fold_df[["Gene Name", "'papillary thyroid carcinoma' vs 'normal'.log2foldchange"]].dropna()
fold_df_slice.to_csv('chern_tc_fc.rnk', index=False, header=False, sep='\t')


# Most Upregulated
# 1) hsa04650: Natural killer cell mediated cytotoxicity (NES = 1.8662)
# 2) hsa05150: Staphylococcus aureus infection (NES = 1.7967)
# 3) hsa05320: Autoimmune thyroid disease (NES = 1.7433)

# Most Downregulated
# 1) hsa00053: Ascorbate and aldarate metabolism (NES = -1.2084)
# 2) hsa00750: Vitamin B6 metabolism (NES = -1.0648)
# 3) hsa00120: Primary bile acid biosynthesis (NES = -1.0309)

