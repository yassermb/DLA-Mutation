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
#label_classes = ['Ag-Ab', 'EI', 'GC', 'HR', 'OE', 'RC', 'TCR/pMHC']
label_classes = ['AA', 'EI', 'RX']
#label_classes = ['Ag-Ab', 'EI', 'RC', 'TCR/pMHC']
comp_fun_dict = pd.read_csv('skempi_v2.csv', sep=';').dropna(subset=['Hold_out_type'])
comp_fun_dict.replace({'Pr/PI': 'EI', 'AB/AG': 'AA'}, inplace=True)
comp_fun_dict = dict(zip(comp_fun_dict['#Pdb'].str[:4], comp_fun_dict['Hold_out_type']))


#label_classes_proximate = ['Ag-Ab', 'EI', 'GC', 'HR', 'OE', 'RC']
comp_fun_dict_proximate = pd.read_csv('proximate.csv', sep='\t', header=None).dropna(subset=[5])
comp_fun_dict_proximate.replace({'Ag-Ab': 'AA'}, inplace=True)
comp_fun_dict_proximate = dict(zip(comp_fun_dict_proximate[1], comp_fun_dict_proximate[5]))

comp_fun_dict = {**comp_fun_dict, **comp_fun_dict_proximate}

#comp_fun_dict = {k:v for k,v in comp_fun_dict.items() if v in label_classes}
comp_fun_dict_tmp = {}
for k,v in comp_fun_dict.items():
    if v in label_classes:
        comp_fun_dict_tmp[k]=v
    elif v in ['HR','RC','GC']:
        comp_fun_dict_tmp[k]='RX'
    else:
        pass
comp_fun_dict = comp_fun_dict_tmp

labels = list(map(lambda x: [x], label_classes))
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

def get_model(input_shape):
    """
    Model architecture
    """

    X_in = Input(shape=input_shape)
    
    H = Dense(20, use_bias=True, activation='elu', kernel_initializer='he_uniform', kernel_regularizer=l2(1e-3), kernel_constraint=max_norm(4), bias_constraint=max_norm(4))(X_in)
    Y = Dense(len(label_classes), activation='softmax')(H)

    _model = Model(inputs=[X_in], outputs=Y)
    _model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=0.00001))
    _model.summary()
    
    return _model


input_data = 'intermediate_xray_wt_nomask_200.csv'



output_confusionmatrix_svg = 'tf_confusionmatrix_function_classificaiton_' + input_data.replace('.csv','.svg') 
output_precisionrecall_svg = 'tf_precisionrecall_function_classificaiton_' + input_data.replace('.csv','.svg')
output_f1_svg = 'tf_f1_function_classificaiton_' + input_data.replace('.csv','.svg') 
output_trainvalidcurve_svg = 'tf_trainvalcurve_function_classificaiton_' + input_data.replace('.csv','.svg')
output_matrix = 'tf_confusionmatrix_function_classificaiton_' + input_data.replace('.csv','')


#Read and prepare embeddings
l_embed = 200
embedding_columns = list(map(lambda x: 'e'+str(x),list(range(l_embed))))
df = pd.read_csv(input_data, sep='\t')
df[embedding_columns] = df.embeddings.str.split(',', expand=True)
df=df.drop(['embeddings'], axis=1)


#Define the targets
df['target_class'] = df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].apply(lambda x: comp_fun_dict[x] if x in comp_fun_dict else 'NA')
df = df.loc[df.target_class != 'NA']
df['target'] = df.target_class.apply(lambda x: encoder.transform([[x]])[0])


df['scr_features'] = df.resregion.apply(lambda x: list(encoder_scr.transform([[x]])[0]))
df['auxiliary'] = df.scr_features
df['auxiliary'] = df.auxiliary.apply(lambda x: np.array(x).astype(np.float32))



for seq_sim in ['30', '50', '70', '90', '95', '100']:
    df['ht'+str(seq_sim)] = df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].apply(lambda x: ','.join(map(lambda xx: str(xx), skempi_homology[seq_sim][x])))

separation_level = 'ht30'
clusters_set = df[separation_level].unique()
df_train = df.loc[df[separation_level].isin(clusters_set[:int(0.3*len(clusters_set))])]
df_test = df.loc[df[separation_level].isin(clusters_set[int(0.3*len(clusters_set)):])]


#train_complexes = list(set(map(lambda x: x.split('_')[1], pd.read_csv('transfer_learning_aa_train_complexes.csv', sep=',').complex.to_list())))
#test_complexes = list(set(map(lambda x: x.split('_')[1], pd.read_csv('transfer_learning_aa_test_complexes.csv', sep=',').complex.to_list())))

#train_complexes = list(set(train_complexes) - set(test_complexes))
#test_complexes += train_complexes[30:]
#train_complexes = list(set(train_complexes) - set(test_complexes))

#df_train = df.loc[df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].isin(train_complexes)]
#df_test = df.loc[df.complex.str.replace('SKEMPI2_PDBs_','').str[:4].isin(test_complexes)]

###############################################################################################
X_train = []
y_train = []
instances = list(set(df_train.complex.str.replace('SKEMPI2_PDBs_','').str[:4].to_list()))
classes_fun = []
for instance in instances:
    X_train.append(df_train.loc[df_train.complex.str.contains(instance)][embedding_columns].astype(np.float32).mean().to_numpy())
    y_train.append(df_train.loc[df_train.complex.str.contains(instance)].iloc[0].target.astype(np.float32))
    classes_fun.append(df_train.loc[df_train.complex.str.contains(instance)].iloc[0].target_class)
X_train = np.array(X_train)
y_train = np.array(y_train)

X_test = []
y_test = []
instances = list(set(df_test.complex.str.replace('SKEMPI2_PDBs_','').str[:4].to_list()))
for instance in instances:
    X_test.append(df_test.loc[df_test.complex.str.contains(instance)][embedding_columns].astype(np.float32).mean().to_numpy())
    y_test.append(df_test.loc[df_test.complex.str.contains(instance)].iloc[0].target.astype(np.float32))
X_test = np.array(X_test)
y_test = np.array(y_test)
###############################################################################################

class_weights = class_weight.compute_class_weight('balanced',
                                                  np.unique(sorted(classes_fun)),
                                                  sorted(classes_fun))
d_class_weights = dict(enumerate(class_weights))
d_class_weights[0] = 0.4
d_class_weights[0] = 0.7
d_class_weights[2] = 10

print(type(X_train), X_train.shape)
print(type(y_train), y_train.shape)

print(X_train[:2])
print(y_train[:2])

train_input = [X_train]
test_input = [X_test]


model = get_model(200)
history = model.fit(train_input, y_train, validation_data=(test_input, y_test), batch_size=50, epochs=3500, verbose=1, class_weight=d_class_weights)


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

print(instances, y_test, test_preds)

y_test = np.argmax(y_test, axis=1)
test_preds = np.argmax(test_preds, axis=1)

#Comprehensive guide to multiclass classification metrics!
#print(classification_report(y_test, test_preds), target_names=label_classes)

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

labels = np.asarray(labels).reshape(len(label_classes),len(label_classes))

font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 40}
plt.rc('font', **font)
fig, ax = plt.subplots(figsize=(20,20))
sns.heatmap(cf_matrix_percentage, annot=labels, fmt='', cmap='vlag', ax=ax, vmin=0, vmax=1, annot_kws={'fontsize': 30})

#ax.set_title('');
ax.set_xlabel('Predicted functional classes')
ax.set_ylabel('True functional classes');
ax.xaxis.set_ticklabels(label_classes)
ax.yaxis.set_ticklabels(label_classes, rotation=0)
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

performance = list(zip(label_classes*2, recall + precision, ['Recall']*len(label_classes) + ['Precision']*len(label_classes)))
performance_df = pd.DataFrame(performance, columns = ['Classes', 'Performance', 'Metric'])
fig, ax = plt.subplots(figsize=(20,20))
sns.barplot(x='Classes', y='Performance', data=performance_df, hue='Metric', ax = ax)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig(output_precisionrecall_svg)

recall = np.array(recall)
precision = np.array(precision)
performance = list(zip(label_classes, 2*recall*precision / (precision + recall), ['F1']*len(label_classes)))
performance_df = pd.DataFrame(performance, columns = ['Classes', 'Performance', 'Metric'])
fig, ax = plt.subplots(figsize=(20,20))
sns.barplot(x='Classes', y='Performance', data=performance_df, hue='Metric', ax = ax)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig(output_f1_svg)
