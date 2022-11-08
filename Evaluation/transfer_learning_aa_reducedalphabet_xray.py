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



#Target: AA reduced alphabet
labels = [['ARO'], ['CAST'], ['PHOB'], ['POS'], ['POL-N'], ['GLY'], ['PRO']]
encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(asarray(labels))

labels_scr = [['S'], ['C'], ['R']]
encoder_scr = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot_scr = encoder_scr.fit(asarray(labels_scr))


skempi_homology = {}
skempi_homology['100'] = load_obj('skempi_homology/skempi_homology_100')
skempi_homology['95'] = load_obj('skempi_homology/skempi_homology_95')
skempi_homology['90'] = load_obj('skempi_homology/skempi_homology_90')
skempi_homology['70'] = load_obj('skempi_homology/skempi_homology_70')
skempi_homology['50'] = load_obj('skempi_homology/skempi_homology_50')
skempi_homology['30'] = load_obj('skempi_homology/skempi_homology_30')

def get_model(input_shape, input_shape_aux):
    """
    Model architecture
    """

    X_in = Input(shape=input_shape)
    aux_input = Input(shape=input_shape_aux)
    
    H = Dense(20, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(X_in)
    #H = Dense(30, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(H)
    #H = Dense(10, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(H)
    H = Concatenate()([H, aux_input])
    Y = Dense(7, activation='softmax')(H)

    _model = Model(inputs=[X_in, aux_input], outputs=Y)
    _model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=0.00001))
    _model.summary()
    
    return _model


input_data = 'intermediate_xray_wt_nomask_200.csv'

output_confusionmatrix_svg = 'tf_confusionmatrix_aareducedalphabet_classificaiton_' + input_data.replace('.csv','.svg') 
output_precisionrecall_svg = 'tf_precisionrecall_aareducedalphabet_classificaiton_' + input_data.replace('.csv','.svg')
output_f1_svg = 'tf_f1_aareducedalphabet_classificaiton_' + input_data.replace('.csv','.svg') 
output_trainvalidcurve_svg = 'tf_trainvalcurve_aareducedalphabet_classificaiton_' + input_data.replace('.csv','.svg')
output_matrix = 'tf_confusionmatrix_aareducedalphabet_classificaiton_' + input_data.replace('.csv','')


#Read and prepare embeddings
l_embed = 200
embedding_columns = list(map(lambda x: 'e'+str(x),list(range(l_embed))))
df = pd.read_csv(input_data, sep='\t')
df[embedding_columns] = df.embeddings.str.split(',', expand=True)
df=df.drop(['embeddings'], axis=1)


#Define the targets
reduced_alphabet = {'PHE':'ARO', 'TRP':'ARO', 'TYR':'ARO', 'HIS':'ARO', 
                    'CYS':'CAST', 'ALA':'CAST', 'SER':'CAST', 'THR':'CAST', 
                    'ILE':'PHOB', 'LEU':'PHOB', 'MET':'PHOB', 'VAL':'PHOB',
                    'LYS':'POS', 'ARG':'POS',
                    'ASN':'POL-N', 'GLN':'POL-N', 'ASP':'POL-N', 'GLU':'POL-N',
                    'GLY':'GLY',
                    'PRO':'PRO'}
df['target_class'] = df.resname.apply(lambda x: reduced_alphabet[x])
df['target'] = df.target_class.apply(lambda x: encoder.transform([[x]])[0])


df['scr_features'] = df.resregion.apply(lambda x: list(encoder_scr.transform([[x]])[0]))
df['auxiliary'] = df.scr_features
df['auxiliary'] = df.auxiliary.apply(lambda x: np.array(x).astype(np.float32))


"""
for seq_sim in ['30', '50', '70', '90', '95', '100']:
    df['ht'+str(seq_sim)] = df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].apply(lambda x: ','.join(map(lambda xx: str(xx), skempi_homology[seq_sim][x])))

separation_level = 'ht30'
clusters_set = df[separation_level].unique()
df_train = df.loc[df[separation_level].isin(clusters_set[:int(0.7*len(clusters_set))])]
df_test = df.loc[df[separation_level].isin(clusters_set[int(0.7*len(clusters_set)):])]
"""

train_complexes = list(set(map(lambda x: x.split('_')[1], pd.read_csv('transfer_learning_aa_train_complexes.csv', sep=',').complex.to_list())))
test_complexes = list(set(map(lambda x: x.split('_')[1], pd.read_csv('transfer_learning_aa_test_complexes.csv', sep=',').complex.to_list())))
df_train = df.loc[df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].isin(train_complexes)]
df_test = df.loc[df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].isin(test_complexes)]


class_weights = class_weight.compute_class_weight('balanced',
                                                  np.unique(sorted(df_train.target_class.to_list())),
                                                  sorted(df_train.target_class.to_list()))
d_class_weights = dict(enumerate(class_weights))



X_train = df_train[embedding_columns].to_numpy().astype(np.float32)
X_train_aux = np.array(list(map(lambda x: x.astype(np.float32), df_train.auxiliary.to_numpy())))
y_train = np.array(list(map(lambda x: x.astype(np.float32), df_train.target.to_numpy())))


X_test = df_test[embedding_columns].to_numpy().astype(np.float32)
X_test_aux = np.array(list(map(lambda x: x.astype(np.float32), df_test.auxiliary.to_numpy())))
y_test = np.array(list(map(lambda x: x.astype(np.float32), df_test.target.to_numpy())))


print(type(X_train), X_train.shape)
print(type(X_train_aux), X_train_aux.shape)
print(type(y_train), y_train.shape)

print(X_train[:2])
print(X_train_aux[:2])
print(y_train[:2])

train_input = [X_train, X_train_aux]
test_input = [X_test, X_test_aux]


model = get_model(200, 3)
history = model.fit(train_input, y_train, validation_data=(test_input, y_test), batch_size=50, epochs=1000, verbose=1, class_weight=d_class_weights)


# --------------------------------------
# Loss functions evolution
# --------------------------------------
plt.figure(figsize=(10,10))
plt.rcParams.update({'font.size': 22})
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss by epoch')
plt.ylabel('Log-Loss')
plt.xlabel('Epochs')
plt.legend(['Train', 'Valid'], loc='right')
plt.savefig(output_trainvalidcurve_svg)


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
    X_test_aux = np.array(list(map(lambda x: x.astype(np.float32), df_test.loc[df_test.complex.str.contains(instance)].auxiliary.to_numpy())))
    test_input = [X_test, X_test_aux]
    test_preds.append(np.mean(model.predict(test_input), axis=0))
    y_test.append(df_test.loc[df_test.complex.str.contains(instance)].iloc[0].target.astype(np.float32))
test_preds = np.array(test_preds)
y_test = np.array(y_test)
########################################################
"""



y_test = np.argmax(y_test, axis=1)
test_preds = np.argmax(test_preds, axis=1)

aa_class_labels = ['ARO','CAST','PHOB','POS','POL-N','GLY','PRO']

#Comprehensive guide to multiclass classification metrics!
#print(classification_report(y_test, test_preds), target_names=aa_class_labels)

cf_matrix = confusion_matrix(y_test, test_preds)

np.save(output_matrix, cf_matrix)

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

labels = np.asarray(labels).reshape(7,7)

font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 40}
plt.rc('font', **font)
fig, ax = plt.subplots(figsize=(20,20))
sns.heatmap(cf_matrix_percentage, annot=labels, fmt='', cmap='vlag', ax=ax, vmin=0, vmax=1, annot_kws={'fontsize': 30})

#ax.set_title('');
ax.set_xlabel('Predicted amino acid classes')
ax.set_ylabel('True amino acid classes');
ax.xaxis.set_ticklabels(aa_class_labels)
ax.yaxis.set_ticklabels(aa_class_labels, rotation=0)
plt.tight_layout()
plt.savefig(output_confusionmatrix_svg)



recall = []
a=cf_matrix/cf_matrix.sum(axis=1)[:,None]
for i, v in enumerate(a):
    recall.append(v[i])

precision = []
b=cf_matrix/cf_matrix.sum(axis=0)[:,None]
for i, v in enumerate(b):
    precision.append(v[i])

performance = list(zip(aa_class_labels*2, recall + precision, ['Recall']*7 + ['Precision']*7))
performance_df = pd.DataFrame(performance, columns = ['Classes', 'Performance', 'Metric'])
fig, ax = plt.subplots(figsize=(20,20))
sns.barplot(x='Classes', y='Performance', data=performance_df, hue='Metric', ax = ax)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig(output_precisionrecall_svg)

recall = np.array(recall)
precision = np.array(precision)
performance = list(zip(aa_class_labels, 2*recall*precision / (precision + recall), ['F1']*7))
performance_df = pd.DataFrame(performance, columns = ['Classes', 'Performance', 'Metric'])
fig, ax = plt.subplots(figsize=(20,20))
sns.barplot(x='Classes', y='Performance', data=performance_df, hue='Metric', ax = ax)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig(output_f1_svg)
