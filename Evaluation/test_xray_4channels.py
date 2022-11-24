#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 21:09:56 2020

@author: yasser
"""

import logging
import os
import sys
import gc
from os import path, mkdir, getenv, listdir, remove, system, stat
import pandas as pd
import numpy as np
#from prody import *
import glob
import shutil
#import matplotlib.pyplot as plt
import seaborn as sns
from math import exp
import subprocess
from subprocess import CalledProcessError, check_call
import traceback
from random import shuffle, random, seed, sample
from numpy import newaxis
import matplotlib.pyplot as plt
import time

import collections
#import scr
from numpy import asarray
from sklearn.preprocessing import OneHotEncoder

import tensorflow.keras
from tensorflow.keras import backend as K
from tensorflow.keras import regularizers
from tensorflow.keras.datasets import mnist # subroutines for fetching the MNIST dataset
from tensorflow.keras.models import Model, Sequential,load_model # basic class for specifying and training a neural network
from tensorflow.keras.layers import Input, Conv2D, MaxPooling2D, Dense, Dropout, Activation, Flatten, AveragePooling3D
#from tensorflow.keras.utils import np_utils # utilities for one-hot encoding of ground truth values

from tensorflow.keras.layers import Dot
from tensorflow.keras.backend import ones, ones_like
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, accuracy_score, roc_auc_score, roc_curve, precision_recall_curve, log_loss
from sklearn.preprocessing import MinMaxScaler
import pickle

import scipy

print('Your python version: {}'.format(sys.version_info.major))
USE_TENSORFLOW_AS_BACKEND = True
# IF YOU *DO* HAVE AN Nvidia GPU on your computer, or execute on Google COLAB, then change below to False!
FORCE_CPU = False #False 
if USE_TENSORFLOW_AS_BACKEND:
    os.environ['KERAS_BACKEND'] = 'tensorflow'
else:
    os.environ['KERAS_BACKEND'] = 'theano'
if FORCE_CPU:
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
if USE_TENSORFLOW_AS_BACKEND == True:
    import tensorflow as tf
    print('Your tensorflow version: {}'.format(tf.__version__))
    print("GPU : "+tf.test.gpu_device_name())
    physical_devices = tf.config.experimental.list_physical_devices('GPU')
    tf.config.experimental.set_memory_growth(physical_devices[0], True)
else:
    import theano
    print('Your theano version: {}'.format(theano.__version__))

logging.basicConfig(filename='manager.log', filemode='w', format='%(levelname)s: %(message)s', level=logging.DEBUG)
mainlog = logging.getLogger('main')
logging.Logger

seed(int(np.round(np.random.random()*10)))
#################################################################################################

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS'}

v_dim = 24

map_dir = '../Examples/map_dir'

output_file = 'output_xray_wt_mask_4channels.csv'
intermediate_file = 'intermediate_xray_wt_mask_200_4channels.csv'

model = load_model(path.join('../Models', 'ssDLA_model_4channels'))

def load_map(sample_test):
    check_call(
        [
            'lz4', '-d', '-f',
            sample_test
        ],
        stdout=sys.stdout)
    X, y, reg_type, res_pos, res_name, inter_info = load_obj(sample_test.replace('.pkl.lz4',''))
    
    #Filter features (SCR and RL)
    X = X[:,:,:,:,:4]

    remove(sample_test.replace('.lz4',''))
    return X, y, reg_type, res_pos, res_name, inter_info
        
samples = glob.glob(path.join(map_dir,'*','*'))


output_handler = open(output_file, 'w')
intermediate_handler = open(intermediate_file, 'w')

output_handler.write('complex' + '\t' + 'resname' + '\t' + 'resregion' + '\t' + 'resnumber' + '\t' + 'respos' + '\t' + 'partner' + '\t' + 'prediction' + '\t' + 'target' + '\t' + 'entropy' + '\t' + 'crossentropy' + '\n')
intermediate_handler.write('complex' + '\t' + 'resname' + '\t' + 'resregion' + '\t' + 'resnumber' + '\t' + 'respos' + '\t' + 'partner' + '\t' + 'embeddings' + '\n')

for sample_test in samples:
    try:
        print(sample_test)
        X, y, reg_type, res_pos, res_name, inter_info = load_map(sample_test)
    except Exception as e:
        logging.info("Bad interface!" + '\nError message: ' + str(e) + 
                      "\nMore information:\n" + traceback.format_exc())
        continue
    
    X = np.array(X)
    y = np.array(y)
    
    if X.shape[0] > 400:
        X = X[:400]
        y = y[:400]
        reg_type = reg_type[:400]
        res_pos = res_pos[:400]
        res_name = res_name[:400]
        inter_info = inter_info[:400]
    
    comp_name=path.basename(sample_test).replace('.pkl.lz4', '')            

    start = time.time()
    y_preds = model.predict([X], batch_size=X.shape[0])
    end = time.time()
        
    intermediate_model = Model(inputs=model.input, outputs=model.get_layer('layer1').output)
    intermediate_prediction = intermediate_model.predict([X], batch_size=X.shape[0])
    _ = gc.collect()
    
    
    for i in range(len(X)):
        if y[i].sum() != 1: # No class were assigned to the artificial amino acids like MSE
            continue
        output_handler.write(comp_name + '\t' +
                             res_name[i][0] + '\t' +
                             reg_type[i] + '\t' +
                             str(inter_info[i][2]) + '\t' +
                             ','.join(list(map(lambda x: str(x), res_pos[i]))) + '\t' +
                             inter_info[i][3] + '\t' + 
                             ','.join(list(map(lambda x: str(x), y_preds[i]))) + '\t' + 
                             ','.join(list(map(lambda x: str(x), y[i]))) + '\t' +
                             str(scipy.stats.entropy(y_preds[i])) + '\t' + 
                             str(log_loss(y[i], y_preds[i])) + '\n')
        
        intermediate_handler.write(comp_name + '\t' +
                                   res_name[i][0] + '\t' +
                                   reg_type[i] + '\t' +
                                   str(inter_info[i][2]) + '\t' +
                                   ','.join(list(map(lambda x: str(x), res_pos[i]))) + '\t' +
                                   inter_info[i][3] + '\t' +
                                   ','.join(list(map(lambda x: str(x), intermediate_prediction[i]))) + '\n')
