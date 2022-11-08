#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 14:42:12 2021

@author: yasser
"""

import glob
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Input, Conv3D, MaxPooling3D, AveragePooling3D, Layer, BatchNormalization, Add, Lambda, Dense, Flatten, Concatenate
from tensorflow.keras.models import Model, Sequential, load_model
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend as K
from tensorflow.keras.optimizers import Adam
from subprocess import CalledProcessError, check_call
from random import shuffle, random, seed, sample
import matplotlib.pyplot as plt
from os import path, remove, system
import pickle, sys
from scipy import stats
import shutil
import gc

from sklearn.preprocessing import OneHotEncoder

from random import shuffle, random, seed, sample
seed(int(np.round(np.random.random()*10)))

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def load_map(filename):
    #system('copy ' + filename + ' file.pkl.lz4')
    #shutil.copyfile(filename, 'file.pkl.lz4')
    """
    check_call(
        [
            'lz4_win64_v1_9_3\lz4.exe', '-d', '-f',
            'file.pkl.lz4'
        ],
        stdout=sys.stdout)
    """
    #X_wt, X_mut, y, scr, _, _ = load_obj('file')
    X_wt, X_mut, y, scr, _, _ = load_obj(filename.replace('.pkl',''))
    
    
    X_wt = X_wt[:,:,:,:,:167]
    X_mut = X_mut[:,:,:,:,:167]
    
    
    #remove('file.pkl')
    #remove('file.pkl.lz4')
    return X_wt, X_mut, y, scr

path_training = '../finetune_mapping_scr_gemme_jet/'

gemme_jet_dict = load_obj(path.join(path_training, '../comp_mut_map_gemme_jet'))

backrub_models = list(map(lambda x: x.strip(), open('../split/test_interfaces.txt', 'r').readlines()))
backrub_models_test = []
for br_model in backrub_models:
    comp = path.basename(br_model)[:4]
    muta = path.basename(br_model).split('--')[1]
    comp_clustid = 'NA'
    backrub_models_test.append((br_model, comp, comp_clustid, muta))

encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(np.asarray([['SUP'], ['COR'], ['RIM'], ['SUR'], ['INT']]))

model_s = load_model(path.join(path_training,'model_s_50'))

test_preds = []
for backrub_model in backrub_models_test:
    br_model = backrub_model[0]
    comp = backrub_model[1]
    muta = backrub_model[3]
    comp_clustid = backrub_model[2]

    try:
        X_wt, X_mut, y, scr = load_map(br_model)
        gemme_value, pc_value, tr_value, freq_value, trace_value = gemme_jet_dict[path.basename(path.dirname(br_model))]
    except:
        continue
    if scr not in ['SUP', 'COR', 'RIM', 'SUR', 'INT']:
        continue
    scr = list(encoder.transform(np.array([scr]).reshape(-1,1)).reshape(-1))
    aux_feat = np.array(scr + [gemme_value, pc_value, tr_value, freq_value, trace_value])
    #aux_feat = np.array(scr)

    print(aux_feat.shape)
    aux_feat = np.array([aux_feat])
    print(aux_feat.shape)
    print('X_wt', X_wt.shape)
    
    test_preds.append((model_s.predict([X_wt, X_mut, aux_feat])[0][0], y, comp, comp_clustid, br_model)) 
    

save_obj(test_preds, 'test_preds')
