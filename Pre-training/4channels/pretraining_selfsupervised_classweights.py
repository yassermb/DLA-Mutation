#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 21:09:56 2020

@author: yasser
"""

import logging
import os
import sys
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
from scipy import interp

import collections
#import scr
from numpy import asarray
from sklearn.preprocessing import OneHotEncoder

import tensorflow.keras
from tensorflow.keras import backend as K
from tensorflow.keras import regularizers
from tensorflow.keras.datasets import mnist # subroutines for fetching the MNIST dataset
from tensorflow.keras.models import Model, Sequential,load_model # basic class for specifying and training a neural network
from tensorflow.keras.layers import Input, Conv2D, MaxPooling2D, Dense, Dropout, Activation, Flatten, AveragePooling3D, Concatenate
from tensorflow.keras.constraints import max_norm
#from tensorflow.keras.utils import np_utils # utilities for one-hot encoding of ground truth values

from tensorflow.keras.layers import Dot
from tensorflow.keras.backend import ones, ones_like
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, accuracy_score, roc_auc_score, roc_curve, precision_recall_curve, log_loss
from sklearn.preprocessing import MinMaxScaler
import pickle
import scipy

from sklearn.preprocessing import OneHotEncoder

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

NB_EPOCH = 30
PATIENCE = 1000
hidden_size1 = 200
hidden_size2 = 20
v_dim = 24

logging.basicConfig(filename='manager.log', filemode='w', format='%(levelname)s: %(message)s', level=logging.DEBUG)
mainlog = logging.getLogger('main')
logging.Logger

seed(int(np.round(np.random.random()*10)))

map_dir_pdb = 'map_dir_pdb'
map_dir_res = 'map_dir_res'

print('Your python version: {}'.format(sys.version_info.major))
# Uncomment lines below only if you need them 
#!{sys.executable} -m pip install -U numpy --user
#!{sys.executable} -m pip install -U matplotlib --user
#!{sys.executable} -m pip install -U keras --user
#!{sys.executable} -m pip install -U tensorflow --user
#!{sys.executable} -m pip install -U theano --user


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
    
    

#print('Your keras version: {}'.format(keras.__version__))
if USE_TENSORFLOW_AS_BACKEND == True:
    import tensorflow
    print('Your tensorflow version: {}'.format(tensorflow.__version__))
    print("GPU : "+tensorflow.test.gpu_device_name())
else:
    import theano
    print('Your theano version: {}'.format(theano.__version__))    
    

from tensorflow.keras.layers import Conv3D, MaxPooling3D, Layer, BatchNormalization, Add, Lambda
import tensorflow as tf

#physical_devices = tf.config.experimental.list_physical_devices('GPU')
#tf.config.experimental.set_memory_growth(physical_devices[0], True)


"""
if path.isdir(map_dir_res):
    shutil.rmtree(map_dir_res)
mkdir(map_dir_res)
list_aa = ['CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 
           'THR', 'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 
           'TRP', 'ALA', 'VAL', 'GLU', 'TYR', 'MET']
for aa in list_aa:
    mkdir(path.join(map_dir_res, aa))
"""

encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(np.asarray([['S'], ['C'], ['R']]))

def Conv_3D_model(input_shape):
    
    X_in = Input(shape=input_shape)
    #aux_input = Input(shape=input_shape_aux)
    
    #w = tf.Variable(name="custom_weight", initial_value=tf.random.normal(shape=(167,15)))
    #X_ret = Lambda(lambda x: tf.matmul(x, w))(X_in)
    #X_ret = Conv3D(15, kernel_size=(1, 1, 1), padding = 'valid', activation='linear', kernel_initializer='he_uniform', input_shape=X_in.shape)(X_in)
    
    H = Conv3D(20, kernel_size=(1, 1, 1), use_bias = True, padding = 'valid', activation='linear', kernel_initializer='he_uniform', input_shape=X_in.shape)(X_in)
    H = BatchNormalization()(H)
    #H = MaxPooling3D(pool_size=(2, 2, 2))(H)
    H = Conv3D(20, kernel_size=(3, 3, 3), use_bias = True, padding = 'valid', activation='elu', kernel_initializer='he_uniform', input_shape=H.shape)(H)
    H = BatchNormalization()(H)
    #H = MaxPooling3D(pool_size=(2, 2, 2))(H)
    H = Conv3D(30, kernel_size=(4, 4, 4), use_bias = True, padding = 'valid', activation='elu', kernel_initializer='he_uniform', input_shape=H.shape)(H)
    H = BatchNormalization()(H)
    #H = MaxPooling3D(pool_size=(2, 2, 2))(H)
    H = Conv3D(20, kernel_size=(4, 4, 4), use_bias = True, padding = 'valid', activation='elu', kernel_initializer='he_uniform', input_shape=H.shape)(H)
    H = BatchNormalization()(H)
    H = AveragePooling3D(pool_size=(4, 4, 4), strides=(4, 4, 4))(H)
    H = Flatten()(H)
    H = Dropout(0.4)(H)
    
    #H = Concatenate()([H, aux_input])
    
    H = Dense(hidden_size1, activation='elu', name='layer1', kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(H)
    H = Dropout(0.2)(H)
    
    H = Dense(hidden_size2, activation='elu', name='layer2', kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(H)
    H = Dropout(0.1)(H)
    
    Y = Dense(20, activation='softmax')(H)
        
    #Y = Lambda(lambda x: 0.5*(x+1))(H)
    _model = Model(inputs=[X_in], outputs=Y)
    _model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=0.0001))
    _model.summary()
    return _model

def load_map(train_path):
    check_call(
        [
            'lz4_win64_v1_9_3\lz4.exe', '-d', '-f',
            train_path
        ],
        stdout=sys.stdout)
    X_train, y_train, reg_type, res_pos, res_name, inter_info = load_obj(train_path.replace('.pkl.lz4',''))
    
    
    
    #Filtering Support regions
    """
    X_train_tmp = []
    reg_type_tmp = []
    y_train_tmp = []
    res_pos_tmp = []
    res_name_tmp = []
    inter_info_tmp = []
    for i in range(min(len(X_train),len(reg_type))):
        if reg_type[i] == 'C' or reg_type[i] == 'R':
            X_train_tmp.append(X_train[i])
            reg_type_tmp.append(reg_type[i])
            y_train_tmp.append(y_train[i])
            res_pos_tmp.append(res_pos[i])
            res_name_tmp.append(res_name[i])
            inter_info_tmp.append(inter_info[i])
    reg_type = reg_type_tmp
    X_train = np.array(X_train_tmp)
    y_train = y_train_tmp
    res_pos = res_pos_tmp
    res_name = res_name_tmp
    inter_info = inter_info_tmp
    """
    
    
    
    #Filter features (SCR and RL)
    X_train = X_train[:,:,:,:,:167]
    
    
    remove(train_path.replace('.lz4',''))
    return X_train, y_train, reg_type, res_name, inter_info


samples = listdir(map_dir_pdb)

fhandler_train = open('train_log', 'w')
fhandler_test = open('test_log', 'w')


#List of pair cluster ids for specific complexes (such as Covid RBD and extreme cases of SKEMPI) to be excluded from the training!
######################################################################
list_clustid = [(178, 440),  #2PCC_A_B
                (2992, 3676),#3SCJ_A_E
                (2889, 264), #1EAW_A_B
                (392, 1596), #1BRS_A_D
                (3746, 8),   #1S1Q_A_B
                (2640, 13008)] #1IAR_A_B
def load_chains_cluster(file_pdb_cluster):
    """Reads all PDB chains in clusters.

    :file_pdb_cluster: File containing all PDB chains categorized in different clusters
    :returns: List of list of PDB chains in cluster

    """
    with open(file_pdb_cluster,'r') as f_handler:
        pdb_clusters = f_handler.readlines()
        pdb_clusters = [cluster.split() for cluster in pdb_clusters]
    return pdb_clusters
file_pdb_cluster = "bc-70.out"
clusters = load_chains_cluster(file_pdb_cluster)
def get_cluster_id(clusters, p1,p2):
    p1_c_id = -1
    p2_c_id = -1
    for c_id, cluster in enumerate(clusters):
        if p1.upper() in cluster:
            p1_c_id = c_id
            break
    for c_id, cluster in enumerate(clusters):
        if p2.upper() in cluster:
            p2_c_id = c_id
            break
    return (p1_c_id, p2_c_id)
#####################################################################

#Class weights
####################################################################
d_class_weights=load_obj('class_weights')
d_class_weights[14] = 0.8
d_class_weights[7] = 0.7
d_class_weights[10] = 0.58
####################################################################

for foldk in range(5):
    seed(int(np.round(np.random.random()*10)))
    
    #input_shape=(v_dim,v_dim,v_dim,167+6)
    input_shape=(v_dim,v_dim,v_dim,167)
    model  = Conv_3D_model(input_shape)
    #model = load_model('0_30_model')
    
    shuffle(samples)
    batch_samples_train_1, batch_samples_test_1 = [], []
    for pair in samples[:int(0.9*len(samples))]:
        map_dir = map_dir_pdb
        positives = listdir(path.join(map_dir,pair,'1'))
        for pos in positives:
            
            #######################################################
            p1,p2 = pos.split('--')[0], pos.split('--')[1]
            p1 = p1.split('_')[0].upper() + '_' + p1.split('_')[3]
            p2 = p2.split('_')[0].upper() + '_' + p2.split('_')[3]
            p1_c_id, p2_c_id = get_cluster_id(clusters, p1,p2)
            if (p1_c_id, p2_c_id) in list_clustid or (p2_c_id, p1_c_id) in list_clustid:
                continue
            ######################################################
            
            
            batch_samples_train_1.append(path.join(map_dir,pair,'1',pos))      
    shuffle(batch_samples_train_1)
    
    for pair in samples[int(0.9*len(samples)):]:
        map_dir = map_dir_pdb
        positives = listdir(path.join(map_dir,pair,'1'))
        for pos in positives:
            batch_samples_test_1.append(path.join(map_dir,pair,'1',pos))      
    shuffle(batch_samples_test_1)
    
    #batch_samples_train = batch_samples_train_1[:2000]
    #batch_samples_test = batch_samples_test_1[:100]
    #batch_samples_train = batch_samples_train_1
    #batch_samples_test = batch_samples_test_1
    
    
    """
    with open(str(foldk) + '_train_interfaces.txt', 'w') as f_handler_trainlist:
        for inter in batch_samples_train:
            f_handler_trainlist.write(inter+'\n')
    with open(str(foldk) + '_test_interfaces.txt', 'w') as f_handler_testlist:
        for inter in batch_samples_test:
            f_handler_testlist.write(inter+'\n')
    """
    with open(str(foldk) + '_train_interfaces.txt', 'r') as f_handler_trainlist:
        batch_samples_train = list(map(lambda x: x.replace('\n',''), f_handler_trainlist.readlines()))
    with open(str(foldk) + '_test_interfaces.txt', 'r') as f_handler_testlist:
        batch_samples_test = list(map(lambda x: x.replace('\n',''), f_handler_testlist.readlines()))
    
    
    shuffle(batch_samples_train)
    shuffle(batch_samples_test)
    
    print(len(batch_samples_train))
    print(len(batch_samples_test_1))
    
    step = 1
    for epoch in range(1, NB_EPOCH+20+1):
        
        fhandler_train_e = open('train_log_'+str(foldk)+'_'+str(epoch), 'w')
        fhandler_test_e = open('test_log_'+str(foldk)+'_'+str(epoch), 'w')
        fhandler_train_e.write('complex' + '\t' + 'resname' + '\t' + 'resregion' + '\t' + 'resnumber' + '\t' + 'partner' + '\t' + 'prediction' + '\t' + 'target' + '\t' + 'entropy' + '\t' + 'crossentropy' + '\n')
        fhandler_test_e.write('complex' + '\t' + 'resname' + '\t' + 'resregion' + '\t' + 'resnumber' + '\t' + 'partner' + '\t' + 'prediction' + '\t' + 'target' + '\t' + 'entropy' + '\t' + 'crossentropy' + '\n')
        
        epoch_log_test_ll = []
        epoch_log_train_ll = []
        epoch_log_val_ll = []
        for batch_i in range(0, len(batch_samples_train), step):
            try:
                t = time.time()
                train_path1  = batch_samples_train[batch_i]
                
                X_train1, y_train1, reg_type1, res_name1, inter_info1 = load_map(train_path1)                
                if len(X_train1)==0:
                    continue
                X = np.array(X_train1)
                y = np.array(y_train1)
                res_name = res_name1
                
                reg_type1_value = encoder.transform(list(map(lambda x: [x], reg_type1)))
                
                X_aux = np.array(reg_type1_value)
                
                X_train = X
                y_train = y
                X_train_aux = X_aux
                
                if len(X_train_aux) != len(X_train) or len(res_name) != len(X_train):
                    raise Exception()
                    
                """    
                if foldk == 0 and epoch == 1:
                    for i in range(len(X)):
                        save_obj((X[i], y[i], res_name[i][0]), path.join(map_dir_res, res_name[i][0], res_name[i][0] + '_' + path.basename(train_path1.replace('.pkl.lz4', '')) + '_' + str(i)))
                """
                
                X_val = np.array([X[0]])
                y_val = np.array([y[0]])
                X_val_aux = np.array([X_aux[0]])
                
                #X_train, X_val, y_train, y_val = train_test_split(X, np.array(y_train1+y_train2), test_size=0.1, random_state=int(np.round(np.random.random()*10)))
                
            except:
                #dddd
                print("pass")
                continue
                
            if X_train.shape[0] > 120:
                X_train = X_train[:120]
                y_train = y_train[:120]
                X_train_aux = X_train_aux[:120]
                res_name = res_name[:120]
                
            model.fit([X_train], y_train, batch_size=X_train.shape[0], epochs=1, verbose = 0, class_weight=d_class_weights)
            
            start = time.time()
            train_preds = model.predict([X_train], batch_size=X_train.shape[0])
            end = time.time()
            
            val_preds = model.predict([X_val], batch_size=X_val.shape[0])
            
            
            train_ll = log_loss(y_train, train_preds)
            val_ll = log_loss(y_val, val_preds)
            epoch_log_train_ll.append(train_ll)
            epoch_log_val_ll.append(val_ll)
            print("Epoch: {:04d}".format(epoch),
                  "train_ll= {:.4f}".format(train_ll),
                  "val_ll= {:.4f}".format(val_ll),
                  "time= {:.4f}".format(time.time() - t))
            
            
            fhandler_train.write(train_path1 + '\t' +
                                 str(foldk)  + '\t' +
                                 str(epoch)  + '\t' +
                                 ','.join(reg_type1) + '\t' +
                                 str(end-start) + '\t' +
                                 str(train_ll) + '\n')
            for i in range(len(train_preds)):
                if y_train[i].sum() != 1: #No class were assigned to the artificial amino acids like MSE
                    continue
                fhandler_train_e.write(path.basename(train_path1).replace('.pkl.lz4', '') + '\t' +
                                       res_name[i][0] + '\t' +
                                       reg_type1[i] + '\t' +
                                       str(inter_info1[i][2]) + '\t' +
                                       inter_info1[i][3] + '\t' +
                                       ','.join(list(map(lambda x: str(x), train_preds[i]))) + '\t' + 
                                       ','.join(list(map(lambda x: str(x), y_train[i]))) + '\t' +
                                       str(scipy.stats.entropy(train_preds[i])) + '\t' + 
                                       str(log_loss(y_train[i], train_preds[i])) + '\n')
            #fhandler_train_e.write('\n\n')
            
            
            
        model.save(str(foldk)+'_'+str(epoch)+'_model')
            
        test_preds = []
        Xts_lbl = []
        for test_interface in batch_samples_test:
            try:
                X_test, y_test, reg_type, res_name, inter_info = load_map(test_interface)
                X_aux = encoder.transform(list(map(lambda x: [x], reg_type)))
                if len(X_test) == 0 or len(X_aux) != len(X_test) or len(res_name) != len(X_test):
                    continue
            except:
                continue
            
            
            start = time.time()
            test_preds = model.predict([X_test], batch_size=X_test.shape[0])
            end = time.time()
            test_ll = log_loss(y_test, test_preds)
            epoch_log_test_ll.append(test_ll)
            
            
            
            fhandler_test.write(test_interface + '\t' +
                                 str(foldk)  + '\t' +
                                 str(epoch)  + '\t' +
                                 ','.join(reg_type) + '\t' +
                                 str(end-start) + '\t' +
                                 str(test_ll) + '\n')
            for i in range(len(test_preds)):
                if y_test[i].sum() != 1: #No class were assigned to the artificial amino acids like MSE
                    continue
                fhandler_test_e.write(path.basename(test_interface).replace('.pkl.lz4', '') + '\t' +
                                       res_name[i][0] + '\t' +
                                       reg_type[i] + '\t' +
                                       str(inter_info[i][2]) + '\t' +
                                       inter_info[i][3] + '\t' + 
                                     ','.join(list(map(lambda x: str(x), test_preds[i]))) + '\t' + 
                                     ','.join(list(map(lambda x: str(x), y_test[i]))) + '\t' +
                                     str(scipy.stats.entropy(test_preds[i])) + '\t' + 
                                     str(log_loss(y_test[i], test_preds[i])) + '\n')
            #fhandler_test_e.write('\n\n')
            
            
            
            
        np.save('epoch_log_train_ll_'+str(foldk)+'_'+str(epoch), np.array(epoch_log_train_ll))
        np.save('epoch_log_val_ll_'+str(foldk)+'_'+str(epoch), np.array(epoch_log_val_ll))
        np.save('epoch_log_test_ll_'+str(foldk)+'_'+str(epoch), np.array(epoch_log_test_ll))
        fhandler_train_e.close()
        fhandler_test_e.close()