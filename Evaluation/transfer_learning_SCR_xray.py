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
from tensorflow.keras.layers import Input, Conv3D, MaxPooling3D, AveragePooling3D, Layer, BatchNormalization, Add, Lambda, Dense, Concatenate, Flatten
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend as K
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.constraints import max_norm
from subprocess import CalledProcessError, check_call
from random import shuffle, random, seed, sample
import matplotlib.pyplot as plt
from os import path, remove, system, mkdir
import pickle, sys
from scipy import stats
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.utils import class_weight
import seaborn as sns

from random import shuffle, random, seed, sample
seed(int(np.round(np.random.random()*10)))

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)



#Target: SCR
labels = [['SUP'], ['COR'], ['RIM'], ['INT'], ['SUR']]
encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(asarray(labels))


skempi_homology = {}
skempi_homology['100'] = load_obj('skempi_homology/skempi_homology_100')
skempi_homology['95'] = load_obj('skempi_homology/skempi_homology_95')
skempi_homology['90'] = load_obj('skempi_homology/skempi_homology_90')
skempi_homology['70'] = load_obj('skempi_homology/skempi_homology_70')
skempi_homology['50'] = load_obj('skempi_homology/skempi_homology_50')
skempi_homology['30'] = load_obj('skempi_homology/skempi_homology_30')

def get_model(input_shape):
    """
    Model architecture
    """

    X_in = Input(shape=input_shape)
    
    H = Dense(100, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(X_in)
    H = Dense(50, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(H)
    H = Dense(10, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(H)
    Y = Dense(5, activation='softmax')(H)

    _model = Model(inputs=[X_in], outputs=Y)
    _model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=0.00001))
    _model.summary()
    
    return _model


#Read and prepare embeddings
l_embed = 200
embedding_columns = list(map(lambda x: 'e'+str(x),list(range(l_embed))))
df = pd.read_csv('intermediate_xray_skempi_wt_nomask_200.csv', sep='\t')
df[embedding_columns] = df.embeddings.str.split(',', expand=True)
df=df.drop(['embeddings'], axis=1)


#Define the targets
df['target'] = df.resregion.apply(lambda x: encoder.transform([[x]])[0])

class_weights = class_weight.compute_class_weight('balanced',
                                                  np.unique(sorted(df.resregion.to_list())),
                                                  sorted(df.resregion.to_list()))
d_class_weights = dict(enumerate(class_weights))


for seq_sim in ['30', '50', '70', '90', '95', '100']:
    df['ht'+str(seq_sim)] = df.complex.str.replace('SKEMPI2_PDBs_','').apply(lambda x: ','.join(map(lambda xx: str(xx), skempi_homology[seq_sim][x])))

separation_level = 'ht100'
clusters_set = df[separation_level].unique()
df_train = df.loc[df[separation_level].isin(clusters_set[:int(0.3*len(clusters_set))])]
df_test = df.loc[df[separation_level].isin(clusters_set[int(0.3*len(clusters_set)):])]

X_train = df_train[embedding_columns].to_numpy().astype(np.float32)
y_train = np.array(list(map(lambda x: x.astype(np.float32), df_train.target.to_numpy())))

X_test = df_test[embedding_columns].to_numpy().astype(np.float32)
y_test = np.array(list(map(lambda x: x.astype(np.float32), df_test.target.to_numpy())))


print(type(X_train), X_train.shape)
print(type(y_train), y_train.shape)

print(X_train[:2])
print(y_train[:2])



train_input = X_train
test_input = X_test


model = get_model(200)
history = model.fit(train_input, y_train, validation_data=(test_input, y_test), batch_size=50, epochs=200, verbose=1, class_weight=d_class_weights)


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
plt.savefig('tf_trainvalcurve_scr_classificaiton.svg')


train_preds = model.predict(train_input)


test_preds = model.predict(test_input)
"""
#######################################################
#Averaging embeddings over backrub models
test_preds = []
y_test = []
instances = set(map(lambda x: '--'.join(x.split('--')[:-1]), df_test.complex.to_list())) #just removing the backrub model code!
for instance in instances:
    X_test = df_test.loc[df_test.complex.str.contains(instance)][embedding_columns].to_numpy().astype(np.float32)
    X_test_aux = np.array(list(map(lambda x: x.astype(np.float32), df_test.loc[df_test.complex.str.contains(instance)].gemme_jet.to_numpy())))
    test_input = [X_test, X_test_aux]
    test_preds.append(np.mean(model.predict(test_input), axis=0))
    y_test.append(df_test.loc[df_test.complex.str.contains(instance)].iloc[0].target.astype(np.float32))
test_preds = np.array(test_preds)
y_test = np.array(y_test)
########################################################
"""


y_test = np.argmax(y_test, axis=1)
test_preds = np.argmax(test_preds, axis=1)

scr_labels = ['CORE','INTERIOR','RIM','SUPPORT','SURFACE']

#Comprehensive guide to multiclass classification metrics!
#print(classification_report(y_test, test_preds), target_names=scr_labels)

cf_matrix = confusion_matrix(y_test, test_preds)

np.save('tf_confusionmatrix_scr_classificaiton', cf_matrix)

group_counts = ["{0:0.0f}".format(value) for value in
                cf_matrix.flatten()]

"""
group_percentages = ["{0:.2%}".format(value) for value in
                     cf_matrix.flatten()/np.sum(cf_matrix)]
"""
cf_matrix_percentage = cf_matrix/cf_matrix.sum(axis=1)[:,None]
group_percentages = ["{0:.2%}".format(value) for value in
                     cf_matrix_percentage.flatten()]

labels = [f"{v1}\n{v2}\n" for v1, v2 in
          zip(group_counts,group_percentages)]

labels = np.asarray(labels).reshape(5,5)


font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 40}
plt.rc('font', **font)
fig, ax = plt.subplots(figsize=(20,20))
sns.heatmap(cf_matrix_percentage, annot=labels, fmt='', cmap='Blues', ax=ax, vmin=0, vmax=1, annot_kws={'fontsize': 40})

#ax.set_title('');
ax.set_xlabel('Predicted region classes')
ax.set_ylabel('True region classes');
ax.xaxis.set_ticklabels(scr_labels, rotation=45)
ax.yaxis.set_ticklabels(scr_labels, rotation=45)
plt.tight_layout()
plt.savefig('tf_confusionmatrix_scr_classificaiton.svg')