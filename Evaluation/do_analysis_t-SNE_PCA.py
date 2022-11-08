# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 02:56:41 2021

@author: awadmin
"""

from sklearn.manifold import TSNE
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

list_aa = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 
           'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 
           'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
list_aa.sort()

list_aa_filter = ['ARG','ASP','GLU', 'GLN', 'ASN','PRO']

hue_order_scr = ['SUP','COR','RIM', 'SUR', 'INT']
hue_color_scr = {'SUP': 'red','COR': '#d4aa00ff','RIM':'blue', 'SUR':'purple', 'INT':'gray'}

hue_order_rl = ['R','L']
hue_color_rl = {'R': 'cyan','L': 'orange'}

                 
hue_color_aa = {}
cp = sns.color_palette('Paired', n_colors=20)
for i, aa in enumerate(list_aa):
    hue_color_aa[aa] = cp[i]
hue_order_aa = list_aa
             
filename_source = 'intermediate_xray_skempi_wt_nomask_200'
#filename_source = 'intermediate_xray_skempi_wt_mask_sidechain_200'
#filename_source = 'intermediate_xray_skempi_wt_mask_sphere5_randomcenter_200'
#filename_source = 'intermediate_backrub_nomask'

#filename_source = 'intermediate_xray_skempi_wt_nomask_200_4channels'
#filename_source = 'intermediate_xray_skempi_wt_mask_sidechain_200_4channels'
#filename_source = 'intermediate_xray_skempi_wt_mask_sphere5_randomcenter_200_4channels'


for algo in ['_tSNE.csv','_PCA.csv']:
    
    filename = filename_source + algo

    df = pd.read_csv(filename, sep='\t')
    
    df['resregion'] = df.resregion.replace({'S': 'SUP', 'C': 'COR', 'R': 'RIM'})
    
    #df = df.loc[df.resname.isin(list_aa_filter)]
    #df = df.loc[df.complex.str.contains('1BRS')]
    
    plt.rcParams.update({'font.size': 24})
    
    sns.jointplot(data=df, x='PC1', y='PC2', hue = 'resname', kind='scatter', height=20, hue_order=hue_order_aa, palette=hue_color_aa)
    plt.tight_layout()
    plt.savefig(filename.replace('.csv', '_aa_scatterplot.png'))
    
    sns.jointplot(data=df, x='PC1', y='PC2', hue = 'resregion', kind='scatter', height=20, hue_order=hue_order_scr, palette=hue_color_scr)
    plt.tight_layout()
    plt.savefig(filename.replace('.csv', '_scr_scatterplot.png'))
    
    sns.jointplot(data=df, x='PC1', y='PC2', hue = 'partner', kind='scatter', height=20, hue_order=hue_order_rl, palette=hue_color_rl)
    plt.tight_layout()
    plt.savefig(filename.replace('.csv', '_rl_scatterplot.png'))
