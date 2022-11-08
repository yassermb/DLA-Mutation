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

model_pretrained = load_model('../models5_nonorm_classweight_Porlineweight08/0_30_model')


def get_siamese_model(input_shape, input_shape_aux):
    """
    Model architecture
    """
    
    # Define the tensors for the two input images
    left_input = Input(shape=input_shape)
    right_input = Input(shape=input_shape)
    aux_input = Input(shape=input_shape_aux)
    
    # Convolutional Neural Network
    model = Sequential()

    model_intermediate = Model(inputs=model_pretrained.input, outputs=model_pretrained.get_layer('layer1').output)
    model.add(model_intermediate)

    model.add(Dense(100, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3))) #For comparison with no-pretraining I need to remove this layer!
    model.add(Dense(10, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3)))
    
    # Generate the encodings (feature vectors) for the two images
    encoded_l = model(left_input)
    encoded_r = model(right_input)
    
    # Add a customized layer to compute the absolute difference between the encodings
    L1_layer = Lambda(lambda tensors:K.abs(tensors[0] - tensors[1]))
    L1_distance = L1_layer([encoded_l, encoded_r])
    
    feat_vector = Concatenate()([L1_distance, aux_input])
        
    # Add a dense layer with a sigmoid unit to generate the similarity score
    prediction = Dense(1, use_bias=True, activation='linear')(feat_vector)
    
    # Connect the inputs with the outputs
    siamese_net = Model(inputs=[left_input,right_input,aux_input],outputs=prediction)

    siamese_net.compile(loss='mean_squared_error', optimizer=Adam(lr=0.001))
    siamese_net.summary()
    
    # return the model
    return siamese_net


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

gemme_jet_dict = load_obj('../comp_mut_map_gemme_jet')

backrub_models = list(map(lambda x: x.strip(), open('../split/train_interfaces.txt', 'r').readlines()))

skempi_homology = {}
skempi_homology['100'] = load_obj('../skempi_homology/skempi_homology_100')
skempi_homology['95'] = load_obj('../skempi_homology/skempi_homology_95')
skempi_homology['90'] = load_obj('../skempi_homology/skempi_homology_90')
skempi_homology['70'] = load_obj('../skempi_homology/skempi_homology_70')
skempi_homology['50'] = load_obj('../skempi_homology/skempi_homology_50')
skempi_homology['30'] = load_obj('../skempi_homology/skempi_homology_30')

encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(np.asarray([['SUP'], ['COR'], ['RIM'], ['SUR'], ['INT']]))

backrub_models_train = []
for br_model in backrub_models:
    comp = path.basename(br_model)[:4]
    muta = path.basename(br_model).split('--')[1]
    comp_clustid = skempi_homology['100'][comp]
    backrub_models_train.append((br_model, comp, comp_clustid, muta))


v_dim = 24
#input_shape=(v_dim,v_dim,v_dim,167+4+2)
input_shape=(v_dim,v_dim,v_dim,167)
model_s = get_siamese_model(input_shape, 10)

NB_EPOCH = 80
#backrub_models_train = backrub_models_train[:80]

step = 40
history_epoch = {'loss':[], 'val_loss':[]}
for epoch in range(1, NB_EPOCH+1):
    print('epoch: '+str(epoch))
    history_batch = {'loss':[], 'val_loss':[]}
    for batch_i in range(0, len(backrub_models_train), step):
        try:
            X_train, X_train_aux, y_train = [],[],[]
            for step_i in range(step):    
                backrub_model = backrub_models_train[batch_i+step_i]
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
                
                X_train.append([X_wt[0], X_mut[0]])
                X_train_aux.append(aux_feat)
                y_train.append(y)
        except:
            continue
    
    
        val_len = int(np.round(0.2*len(X_train)))
        
        X_valid = np.array(X_train[:val_len])
        X_valid_aux = np.array(X_train_aux[:val_len])
        y_valid = np.array(y_train[:val_len])
        
        X_train = np.array(X_train[val_len:])
        X_train_aux = np.array(X_train_aux[val_len:])
        y_train = np.array(y_train[val_len:])

        #print(X_train.shape)
        #print(X_train[:,0].shape)
        #print(X_train[:,1].shape)
        #print(np.array([X_train[:,0], X_train[:,1]]).shape)
    
        history = model_s.fit([X_train[:,0], X_train[:,1], X_train_aux], y_train, validation_data=([X_valid[:,0], X_valid[:,1], X_valid_aux], y_valid), batch_size=X_train.shape[0], epochs=1, verbose=1)
        
        history_batch['loss'].append(history.history['loss'][0])
        history_batch['val_loss'].append(history.history['val_loss'][0])
        _ = gc.collect()

    history_epoch['loss'].append(np.array(history_batch['loss']).mean())
    history_epoch['val_loss'].append(np.array(history_batch['val_loss']).mean())
    
    model_s.save('model_s_'+str(epoch))
    save_obj(history_epoch, 'history_epoch_'+str(epoch))
