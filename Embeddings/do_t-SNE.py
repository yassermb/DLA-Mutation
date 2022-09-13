# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 02:56:41 2021

@author: awadmin
"""

from sklearn.manifold import TSNE
import pandas as pd
import numpy as np

#filename = 'intermediate_xray_skempi_wt_nomask_200.csv'
#filename = 'intermediate_mutation_masked-all-aa.csv'

filename = 'intermediate_xray_skempi_wt_nomask_200.csv'

filename_out = filename.replace('.csv', '_tSNE.csv')

l_embed = 200

df = pd.read_csv(filename, sep='\t')
df[list(range(l_embed))] = df['embeddings'].str.split(',', expand=True)
X_embedded = TSNE(n_components=2,random_state=42, perplexity=500, early_exaggeration=1).fit_transform(df[range(0,l_embed)])
df['PC1'] = X_embedded[:,0]
df['PC2'] = X_embedded[:,1]
df=df.drop(['embeddings'], axis=1)
df.to_csv(filename_out, sep='\t')