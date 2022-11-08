# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 23:39:44 2022

@author: awadmin
"""

import pandas as pd
import numpy as np
from sklearn.utils.class_weight import compute_class_weight
import pickle

list_aa = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 
           'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 
           'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
list_aa.sort()

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

a=pd.read_csv('env_mask5_50epoches_completedatabase_nonorm_good/train_log_0_1',sep='\t')
y=list(map(lambda x: list(map(lambda y: float(y), x.split(','))), a.target.to_list()))

y_integers = np.argmax(y, axis=1)
class_weights = compute_class_weight('balanced', np.unique(y_integers), y_integers)
d_class_weights = dict(enumerate(class_weights))


print(list(zip(list_aa, d_class_weights.values())))

save_obj(d_class_weights, 'class_weights')

