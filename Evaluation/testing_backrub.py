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
from sklearn.metrics import mean_squared_error, accuracy_score, roc_auc_score, roc_curve, precision_recall_curve
from sklearn.preprocessing import MinMaxScaler
import pickle

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

#map_dir = 'validations'
#map_dir = 'map_dir_backrub_dynamics'
#map_dir = 'skempi_mutations_mask/map_dir_skempi_wt'
#map_dir = 'skempi_mutations_mask/map_dir_mut_sep'
map_dir = '/home/yasser/myThesis/ProtPartDisc/python_scripts/map_dir_mut_sep_nomask'

output_file = 'output_backrub_nomask.csv'
intermediate_file = 'intermediate_backrub_nomask_200.csv'

model = load_model(path.join('models5_nonorm_classweight_Porlineweight08', '0_30_model'))

def load_map(sample_test):
    check_call(
        [
            'lz4', '-d', '-f',
            sample_test
        ],
        stdout=sys.stdout)
    #X, y, y_ddg = load_obj(sample_test.replace('.pkl.lz4',''))
    X, y, y_ddg, region, comp_type, expr_method = load_obj(sample_test.replace('.pkl.lz4',''))    
    
    #Filter features (SCR and RL)
    X = X[:,:,:,:,:167]
    
    
    remove(sample_test.replace('.lz4',''))
    
    #return X, y, y_ddg
    return X, y, y_ddg, region, comp_type, expr_method
        
samples = glob.glob(path.join(map_dir,'*','1','*'))


output_handler = open(output_file, 'w')
intermediate_handler = open(intermediate_file, 'w')

output_handler.write('complex' + '\t' + 'resname' + '\t' + 'resregion' + '\t' + 'resnumber' + '\t' + 'partner' + '\t' + 'type' + '\t' + 'prediction' + '\t' + 'target' + '\t' + 'ddg' + '\n')
intermediate_handler.write('complex' + '\t' + 'resname' + '\t' + 'resregion' + '\t' + 'resnumber' + '\t' + 'partner' + '\t' + 'type' + '\t' + 'embeddings' + '\t' + 'ddg' + '\n')

for sample_test in samples:
    try:
        print(sample_test)
        #X, y, y_ddg = load_map(sample_test)
        X, y, y_ddg, region, comp_type, expr_method = load_map(sample_test)
    except Exception as e:
        logging.info("Bad interface!" + '\nError message: ' + str(e) + 
                      "\nMore information:\n" + traceback.format_exc())
        continue
    
    
    comp_name=path.basename(sample_test).replace('.pkl.lz4', '')            
    if comp_name.split('_')[0] == 'wt':
        res_name = three_letter[comp_name.split('--')[1][0]]
    else:
        res_name = three_letter[comp_name.split('--')[1][-1]]
    inter_info = ('NA','NA',comp_name.split('--')[1][2:-1],comp_name.split('--')[1][1]) #Must be fixed
    #reg_type = 'NA'  #Must be fixed  (the information must have been embedded in the cube generation generate_ddg2.py)
    reg_type = region
    

    start = time.time()
    y_preds = model.predict([X], batch_size=X.shape[0])[0]
    end = time.time()
        
    intermediate_model = Model(inputs=model.input, outputs=model.get_layer('layer1').output)
    intermediate_prediction = intermediate_model.predict([X], batch_size=X.shape[0])[0]
    _ = gc.collect()
    
    y = y[0]
    
    output_handler.write(comp_name + '\t' +
                         res_name + '\t' +
                         reg_type + '\t' +
                         str(inter_info[2]) + '\t' +
                         inter_info[3] + '\t' +
                         comp_type + '\t' + 
                         ','.join(list(map(lambda x: str(x), y_preds))) + '\t' + 
                         ','.join(list(map(lambda x: str(x), y))) + '\t' +
                         str(y_ddg) + '\n')
    
    intermediate_handler.write(comp_name + '\t' +
                               res_name + '\t' +
                               reg_type + '\t' +
                               str(inter_info[2]) + '\t' +
                               inter_info[3] + '\t' +
                               comp_type + '\t' +
                               ','.join(list(map(lambda x: str(x), intermediate_prediction))) + '\t' + 
                               str(y_ddg) + '\n')