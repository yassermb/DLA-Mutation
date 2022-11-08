#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 14:42:12 2021

@author: yasser
"""

import glob
import pandas as pd
import numpy as np
from numpy import asarray
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Input, Conv3D, MaxPooling3D, AveragePooling3D, Layer, BatchNormalization, Add, Lambda, Dense, Flatten, Concatenate
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend as K
from tensorflow.keras.optimizers import Adam
from subprocess import CalledProcessError, check_call
from random import shuffle, random, seed, sample
import matplotlib.pyplot as plt
from os import path, remove, system, mkdir
import pickle, sys
from scipy import stats
from sklearn.preprocessing import OneHotEncoder

from random import shuffle, random, seed, sample
seed(int(np.round(np.random.random()*10)))

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

labels_scr = [['SUP'], ['COR'], ['RIM'], ['INT'], ['SUR']]
encoder_scr = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot_scr = encoder_scr.fit(asarray(labels_scr))

skempi_homology = {}
skempi_homology['100'] = load_obj('skempi_homology/skempi_homology_100')
skempi_homology['95'] = load_obj('skempi_homology/skempi_homology_95')
skempi_homology['90'] = load_obj('skempi_homology/skempi_homology_90')
skempi_homology['70'] = load_obj('skempi_homology/skempi_homology_70')
skempi_homology['50'] = load_obj('skempi_homology/skempi_homology_50')
skempi_homology['30'] = load_obj('skempi_homology/skempi_homology_30')

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
        
    model.add(Dense(100, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3)))
    model.add(Dense(50, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3)))
    model.add(Dense(10, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3)))
    
    # Generate the encodings (feature vectors) for the two inputs
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

    siamese_net.compile(loss='mean_squared_error', optimizer=Adam(lr=0.01))
    siamese_net.summary()
    
    # return the model
    return siamese_net


df = pd.read_csv('../IntermediateExtraction/getIntermediate/intermediate_ddg', header=None)

#Add GEMME, JET2, and SCR features
gemme_jet_dict = load_obj('comp_mut_map_gemme_jet')
df['gemme_jet'] = df[0].str[:-4].apply(lambda x: list(gemme_jet_dict[x]) if x in gemme_jet_dict else [0,0,0,0,0])
df['scr_features'] = df[1].apply(lambda x: list(encoder_scr.transform([[x]])[0]))
df['auxiliary'] = df.gemme_jet + df.scr_features
df['auxiliary'] = df.auxiliary.apply(lambda x: np.array(x).astype(np.float32))

for seq_sim in ['30', '50', '70', '90', '95', '100']:
    df['ht'+str(seq_sim)] = df[0].str[:4].apply(lambda x: ','.join(map(lambda xx: str(xx), skempi_homology[seq_sim][x])))
clusters_set = df['ht100'].unique()

#df = df.loc[df[3] == 'FL']
df_train = df.loc[df['ht100'].isin(clusters_set[:int(0.3*len(clusters_set))])]
df_test = df.loc[df['ht100'].isin(clusters_set[int(0.3*len(clusters_set)):])]

X_train = df_train[range(5, 400+5)]
X_train.columns = range(0,400)
y_train = df_train[4].to_numpy()

X_test = df_test[range(5, 400+5)]
X_test.columns = range(0,400)
y_test = df_test[4].to_numpy()

model_s = get_siamese_model(200, 10)

train_input = [X_train[range(0, 200)].to_numpy(), X_train[range(200, 400)].to_numpy(), np.array(list(map(lambda x: x.astype(np.float32), df_train.auxiliary.to_numpy())))]
test_input = [X_test[range(0, 200)].to_numpy(), X_test[range(200, 400)].to_numpy(), np.array(list(map(lambda x: x.astype(np.float32), df_test.auxiliary.to_numpy())))]

history = model_s.fit(train_input, y_train, validation_data=(test_input, y_test), batch_size=50, epochs=200, verbose=1)


# --------------------------------------
# Loss functions evolution
# --------------------------------------
plt.figure(figsize=(10,10))
plt.rcParams.update({'font.size': 22})
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss by epoch')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'valid'], loc='right')
plt.show()



train_preds = model_s.predict(train_input)
test_preds = model_s.predict(test_input)

train_preds = np.array([x[0] for x in train_preds])
test_preds = np.array([x[0] for x in test_preds])

print(train_preds)
print(test_preds)
np.save('train_preds', train_preds)
np.save('test_preds', test_preds)
np.save('y_train', y_train)
np.save('y_test', y_test)

corr_pearson_train = stats.pearsonr(y_train, train_preds)
corr_spearmanr_train = stats.spearmanr(y_train, train_preds)

corr_pearson_test = stats.pearsonr(y_test, test_preds)
corr_spearman_test = stats.spearmanr(y_test, test_preds)

print('corr_pearson_train', corr_pearson_train)
print('corr_spearmanr_train',corr_spearmanr_train)
print('corr_pearson_test', corr_pearson_test)
print('corr_spearman_test', corr_spearman_test)

plt.rcParams.update({'font.size': 24})
plt.figure(figsize=(12,12))
plt.scatter(y_train, train_preds, color = 'red', marker = 'o', label = 'Train', alpha=0.5)
plt.scatter(y_test, test_preds, color = 'blue', marker = 's', label = 'Test', alpha=0.5)
plt.xlabel('Experimental ddg', fontsize=24)
plt.ylabel('Predicted ddg', fontsize=24)
plt.legend(fontsize=24)
plt.tight_layout()
plt.savefig('scatter_'+str(corr_pearson_train[0])+'_'+str(corr_pearson_test[0])+'.png')



ddddd






































f_handler_trainlist = open('train_interfaces.txt', 'w')
f_handler_testlist = open('test_interfaces.txt', 'w')

comp_train = []
comp_test = []
#for br_model in backrub_models[:2200]:
for index, row in df.iterrows():
    comp = row[0][:4]
    muta = row[0].split('--')[1]
    comp_clustid = skempi_homology['100'][comp]
        
    if ((random() > 0.2 or (comp_clustid in comp_train)) and (comp_clustid not in comp_test)):
        comp_train.append(comp_clustid)
        X_train_l.append(df[range(5,205)])
        y_train.append(y)
        f_handler_trainlist.write(br_model.replace('\u252c', '__')+'\n')
    else:
        comp_test.append(comp_clustid)
        X_test.append(([X_wt[0], X_mut[0]], comp, muta))
        y_test.append(y)
        f_handler_testlist.write(br_model.replace('\u252c', '__')+'\n')
f_handler_trainlist.close()
f_handler_testlist.close()


v_dim = 24
input_shape=(v_dim,v_dim,v_dim,167+3)
model_s = get_siamese_model(input_shape)

X_train_ss = np.array([x[0] for x in X_train])
X_test_ss = np.array([x[0] for x in X_test])


y_train = np.array(y_train)
y_test = np.array(y_test)
print(X_train_ss.shape)
print(X_train_ss[:,0].shape)
print(X_train_ss[:,1].shape)
print(np.array([X_train_ss[:,0], X_train_ss[:,1]]).shape)

history = model_s.fit([X_train_ss[:,0], X_train_ss[:,1]], y_train, validation_data=([X_test_ss[:,0], X_test_ss[:,1]], y_test), batch_size=10, epochs=50, verbose=1)

# evaluate the model
train_mse = model_s.evaluate([X_train_ss[:,0], X_train_ss[:,1]], y_train, verbose=0)
test_mse = model_s.evaluate([X_test_ss[:,0], X_test_ss[:,1]], y_test, verbose=0)
#save_obj(model_s, 'model_s')
print('Train: %.3f, Test: %.3f' % (train_mse, test_mse))
# plot loss during training
plt.title('Loss / Mean Squared Error')
plt.plot(history.history['loss'], label='train')
plt.plot(history.history['val_loss'], label='test')
plt.legend()
plt.savefig('siameseplot.svg')


model_s.save('model_s')


train_preds = model_s.predict([X_train_ss[:,0], X_train_ss[:,1]])
test_preds = model_s.predict([X_test_ss[:,0], X_test_ss[:,1]])

train_preds = np.array([x[0] for x in train_preds])
test_preds = np.array([x[0] for x in test_preds])

print(train_preds)
print(test_preds)
np.save('train_preds', train_preds)
np.save('test_preds', test_preds)
np.save('y_train', y_train)
np.save('y_test', y_test)

corr_pearson_train = stats.pearsonr(y_train, train_preds)
corr_spearmanr_train = stats.spearmanr(y_train, train_preds)

corr_pearson_test = stats.pearsonr(y_test, test_preds)
corr_spearman_test = stats.spearmanr(y_test, test_preds)

print('corr_pearson_train', corr_pearson_train)
print('corr_spearmanr_train',corr_spearmanr_train)
print('corr_pearson_test', corr_pearson_test)
print('corr_spearman_test', corr_spearman_test)

plt.rcParams.update({'font.size': 24})
plt.figure(figsize=(12,12))
plt.scatter(y_train, train_preds, color = 'red', marker = 'o', label = 'Train', alpha=0.5)
plt.scatter(y_test, test_preds, color = 'blue', marker = 's', label = 'Test', alpha=0.5)
plt.xlabel('Experimental ddg', fontsize=24)
plt.ylabel('Predicted ddg', fontsize=24)
plt.legend(fontsize=24)
plt.tight_layout()
plt.savefig('scatter_'+str(corr_pearson_train[0])+'_'+str(corr_pearson_test[0])+'.svg')




C_train = {}
for seq_sim in ['30', '50', '70', '90', '95', '100']:
    C_train[seq_sim] = [skempi_homology[seq_sim][x[1]] for x in X_train]
del X_train


plt.figure(figsize=(12,12))
#plt.scatter(y_train, train_preds, color = 'red', marker = 'o', label = 'Train R=' + str(np.round(corr_pearson_train[0],3)), alpha=0.3)
plt.scatter(y_train, train_preds, facecolors='none', edgecolors='black', marker = 'o', label = 'Train R=' + str(np.round(corr_pearson_train[0],3)))
color_dict = {'30': 'deeppink', '50': 'purple', '70': 'green', '90': 'cyan', '95': 'magenta', '100': 'blue'}
size_dict = {'30': 60, '50': 50, '70': 40, '90': 30, '95': 20, '100': 10}
already_shown = []
previous_seqsim = ''
for seq_sim in ['30', '50', '70', '90', '95', '100']:
    try:
        chosen_index = [indx for indx in range(len(X_test)) if skempi_homology[seq_sim][X_test[indx][1]] not in C_train[seq_sim]]
        
        X_test_ss = np.array([x[0] for x in np.array(X_test)[chosen_index]])
        y_test_ss = np.array(y_test)[chosen_index]
        test_preds = model_s.predict([X_test_ss[:,0], X_test_ss[:,1]])
        test_preds = np.array([x[0] for x in test_preds])


        np.save('test_preds_seqsim'+seq_sim, test_preds)
        np.save('y_test_seqsim'+seq_sim, y_test_ss)


        corr_pearson_test = stats.pearsonr(y_test_ss, test_preds)
        corr_spearman_test = stats.spearmanr(y_test_ss, test_preds)
        
        chosen_index = list(set(chosen_index) - set(already_shown))
        if not chosen_index == []:
            X_test_ss = np.array([x[0] for x in np.array(X_test)[chosen_index]])
            y_test_ss = np.array(y_test)[chosen_index]
            test_preds = model_s.predict([X_test_ss[:,0], X_test_ss[:,1]])
            test_preds = np.array([x[0] for x in test_preds])
            
            
            plt.scatter(y_test_ss, test_preds, color = color_dict[seq_sim], marker = 's', label = 'Test '+seq_sim+'% R='+str(np.round(corr_pearson_test[0],3)))
            already_shown += chosen_index
            previous_seqsim = seq_sim
        else:
            plt.scatter([], [], color = color_dict[previous_seqsim], marker = 's', label = 'Test '+seq_sim+'% R='+str(np.round(corr_pearson_test[0],3)))
    except:
        fffff
plt.xlabel('Experimental ddg', fontsize=24)
plt.ylabel('Predicted ddg', fontsize=24)
plt.title('Fixed trainset, testset based on sequence similarity')
plt.legend(fontsize=24)
plt.tight_layout()
plt.savefig('scatter_fixedtrainset.svg')