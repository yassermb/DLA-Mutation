#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 10:17:06 2021

@author: mohseni
"""

import logging
import pandas as pd
import numpy as np
from os import path, mkdir, remove, listdir
import glob
import sys
import traceback
import shutil
import itertools

sys.path.insert(1, '../lib/')
import tools as tl

pdbind_path = '/home/yasser/mySoftwares/generaldb/pdbbind/pdbbind_v2019_PP/PP/'
pdbind_preprocessed_path = '/home/yasser/mySoftwares/generaldb/pdbbind/pdbbind_v2019_PP/PP_preprocessed'

def split_pdb(orig_path, dest_path):
    pdb_comp = parsePDB(orig_path).select('protein')
    for ch in set(pdb_comp.getChids()):
        pdb_ch = pdb_comp.select('chain ' + ch)
        writePDB(path.join(dest_path, path.basename(orig_path).replace('.ent.pdb','_') + ch + '.pdb'), pdb_ch)
    
def preprocess_pdbbind():
    if path.exists(pdbind_preprocessed_path):
        shutil.rmtree(pdbind_preprocessed_path)
    mkdir(pdbind_preprocessed_path)
    
    comps = listdir(pdbind_path)
    for comp in comps:
        orig_path = path.join(pdbind_path, comp)
        dest_path = path.join(pdbind_preprocessed_path, comp[:4])
        mkdir(dest_path)
        split_pdb(orig_path, dest_path)
        pdb_chains = glob.glob(path.join(dest_path, '*.pdb'))
        #print(pdb_chains)
        pairs = list(itertools.combinations(pdb_chains, 2))
        for pair in pairs:
            #print(comp, pair)
            tl.call_intbuilder('', dest_path, 
                       path.basename(pair[0]).replace('.pdb', ''), path.basename(pair[1]).replace('.pdb', ''), 
                       5, dest_path, False, True, False)

#preprocess_pdbbind()
#dddddd    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
import logging
import os
import sys
from os import path, mkdir, getenv, listdir, remove, system, stat
import pandas as pd
import numpy as np
from prody import *
import glob
import shutil
#import matplotlib.pyplot as plt
import seaborn as sns
from math import exp
from subprocess import CalledProcessError, check_call, call
import traceback
from random import shuffle, random, seed, sample
from numpy import newaxis
import matplotlib.pyplot as plt
import time
from prody import *
import collections
import scr
from numpy import asarray
from sklearn.preprocessing import OneHotEncoder
import subprocess
from subprocess import CalledProcessError, check_call
import deepScoring.load_data as load
from sklearn.preprocessing import MinMaxScaler

map_dir_pdb = 'map_dir_pdbbind'
inter_dir_pdb = 'inter_dir_pdbbind'
if not path.exists(map_dir_pdb):
    mkdir(map_dir_pdb)
if not path.exists(inter_dir_pdb):
    mkdir(inter_dir_pdb)


bin_path = "./mapsGenerator_mask/build/maps_generator"
v_dim = 24


def save_pdb_interface(comp_inter, name, inter_path):
    #print(comp_inter[0])
    #print(comp_inter[1])
    
    try:
        rec_orig_inter = parsePDB(comp_inter[0]).select('protein')
        rec_orig_inter.setChids('R')
        lig_orig_inter = parsePDB(comp_inter[1]).select('protein')
        lig_orig_inter.setChids('L')
        
        resnum = 'resnum ' + " ".join([str(x[0]) for x in comp_inter[2]])
        rec_orig_inter = rec_orig_inter.select(resnum)
        
        resnum = 'resnum ' + " ".join([str(x[0]) for x in comp_inter[3]])
        lig_orig_inter = lig_orig_inter.select(resnum)
        
        writePDB(path.join(inter_path, name + '.pdb'), rec_orig_inter.toAtomGroup() + lig_orig_inter.toAtomGroup())
    
    except:
        return
    
encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(asarray([['CYS'], ['ASP'], ['SER'], ['GLN'], ['LYS'], ['ILE'], ['PRO'], 
                              ['THR'], ['PHE'], ['ASN'], ['GLY'], ['HIS'], ['LEU'], ['ARG'], 
                              ['TRP'], ['ALA'], ['VAL'], ['GLU'], ['TYR'], ['MET']]))
def mapcomplex(comp_inter, name, map_path):
    print(comp_inter[0])
    print(comp_inter[1])
    
    try:
        rec_orig_inter = parsePDB(comp_inter[0]).select('protein')
        rec_orig_inter.setChids('R')
        lig_orig_inter = parsePDB(comp_inter[1]).select('protein')
        lig_orig_inter.setChids('L')
        
        comp_inter_rec = list(map(lambda x: (x[0], 'R'), comp_inter[2]))
        comp_inter_lig = list(map(lambda x: (x[0], 'L'), comp_inter[3]))
        
        with open('resinfo','w') as fh_res:
            for x in comp_inter_rec:
                fh_res.write(str(x[0])+';'+x[1]+'\n')
            for x in comp_inter_lig:
                fh_res.write(str(x[0])+';'+x[1]+'\n')
        
        
        writePDB(name+'_complex_p1.pdb', rec_orig_inter.toAtomGroup())
        writePDB(name+'_complex_p2.pdb', lig_orig_inter.toAtomGroup())
        writePDB(name+'_complex.pdb', rec_orig_inter.toAtomGroup() + lig_orig_inter.toAtomGroup())
        
        scr.get_scr(name+'_complex_p1.pdb', name+'_complex_p2.pdb', name+'_complex.pdb', name)
        
        
        rimcoresup = pd.read_csv(name+'_rimcoresup.csv', header=None, sep=' ')
        rec_regions = rimcoresup.loc[rimcoresup[4] == 'receptor']
        rec_regions = pd.Series(rec_regions[5].values, index=rec_regions[2]).to_dict()
        lig_regions = rimcoresup.loc[rimcoresup[4] == 'ligand']
        lig_regions = pd.Series(lig_regions[5].values, index=lig_regions[2]).to_dict()
        
        #Unfortunately INTBuilder doesn't distiguish various side chain orientations in the original PDB when it wants
        #to create PDB of docking files! for example we might have two CB atoms in the ARG and in the original PDB they 
        # are name by AARG and BARG! if both of them are name ARG (as INTBuilder causes) then we have a problem with mapGenerator!
        #So here we select those amino acids that has no duplication of side chain atoms!
        #TRIMMING
        
        res_num2name_map_rec = dict(zip(rec_orig_inter.getResnums(),rec_orig_inter.getResnames()))
        res_num2name_map_lig = dict(zip(lig_orig_inter.getResnums(),lig_orig_inter.getResnames()))
        
        
        elements = list(zip(rec_orig_inter.getResnums(), rec_orig_inter.getChids(), 
                            rec_orig_inter.getNames(), rec_orig_inter.getResnames()))
        resnum = set(rec_orig_inter.getResnums()) - set([k[0] for k,v in collections.Counter(elements).items() if (v > 1 or k[0]<0)])
        resnum = 'resnum ' + ' '.join([str(x) for x in resnum])
        rec_orig_inter = rec_orig_inter.select(resnum)
        
        elements = list(zip(lig_orig_inter.getResnums(), lig_orig_inter.getChids(), 
                            lig_orig_inter.getNames(), lig_orig_inter.getResnames()))
        resnum = set(lig_orig_inter.getResnums()) - set([k[0] for k,v in collections.Counter(elements).items() if (v > 1 or k[0]<0)])
        resnum = 'resnum ' + ' '.join([str(x) for x in resnum])
        lig_orig_inter = lig_orig_inter.select(resnum)
        
        
        writePDB(name+'_train.pdb', rec_orig_inter.toAtomGroup() + lig_orig_inter.toAtomGroup())
        
        L1 = list(set(rec_orig_inter.getResnums()))
        res_ind_map_rec = dict([(x,inx) for inx, x in enumerate(sorted(L1))])
        L1 = list(set(lig_orig_inter.getResnums()))
        res_ind_map_lig = dict([(x,inx+len(res_ind_map_rec)) for inx, x in enumerate(sorted(L1))])
        
        
        
        print('comp_inter2', comp_inter_rec)
        print('comp_inter3', comp_inter_lig)
        
        
        res_inter_rec = [(res_ind_map_rec[x[0]], rec_regions[x[0]] if x[0] in rec_regions else 'unknown', x[0], x[1], res_num2name_map_rec[x[0]]) 
                            for x in comp_inter_rec if x[0] in res_ind_map_rec]
        res_inter_lig = [(res_ind_map_lig[x[0]], lig_regions[x[0]] if x[0] in lig_regions else 'unknown', x[0], x[1], res_num2name_map_lig[x[0]])
                            for x in comp_inter_lig if x[0] in res_ind_map_lig]
        res_inter = list(map(lambda x: x[0],res_inter_rec)) + list(map(lambda x: x[0],res_inter_lig))
        reg_type =  list(map(lambda x: x[1],res_inter_rec)) + list(map(lambda x: x[1],res_inter_lig))
        res_name =  list(map(lambda x: [x[4]],res_inter_rec)) + list(map(lambda x: [x[4]],res_inter_lig))
        
        print(res_inter)
        print(reg_type)
        with open('scrinfo','w') as fh_csr:
            for x in res_inter_rec:
                fh_csr.write(str(x[2])+';'+x[3]+';'+x[1]+'\n')
            for x in res_inter_lig:
                fh_csr.write(str(x[2])+';'+x[3]+';'+x[1]+'\n')
    
        if not res_inter:
            return [],[],[]
        
        #tl.coarse_grain_pdb('train.pdb')
        mapcommand = [bin_path, "--mode", "map", "-i", name+'_train.pdb', "--native", "-m", str(v_dim), "-t", "167", "-v", "0.8", "-o", name+'_train.bin']
        call(mapcommand)
        
        dataset_train = load.read_data_set(name+'_train.bin')
        
        
        print(dataset_train.maps.shape)
        
        scaler = MinMaxScaler()
        scaler.fit(dataset_train.maps)
        data_norm = scaler.transform(dataset_train.maps)
        
        X = np.reshape(data_norm, (-1,v_dim,v_dim,v_dim,170))
        #y = [comp_inter[4]]*len(res_inter)
        y = encoder.transform(res_name)
        
        map_name = path.join(map_path, name)
        tl.save_obj((X,y,reg_type), map_name)
        
        check_call(
            [
                'lz4', '-9', '-f',   #, '--rm' because if inconsistency in lz4 versions! 
                map_name + '.pkl'
            ],
            stdout=sys.stdout)
        remove(map_name + '.pkl')
        remove(name+'_train.pdb')
        remove(name+'_train.bin')
        remove(name+'_complex_p1.pdb')
        remove(name+'_complex_p2.pdb')
        remove(name+'_complex.pdb')
        remove(name+'_rimcoresup.csv')
        
        
        print(type(X))
        print(X.shape)
        
        
    except Exception as e:
        #dddd
        print("Bad interface!" + '\ninter_rec: ' +  comp_inter[0] + '\ninter_lig: ' + comp_inter[1] + '\nError message: ' + str(e) + 
                      "\nMore information:\n" + traceback.format_exc())
        print('XXXXXXXXXXXXXXexceptionXXXXXXXXXx')
        return [],[],[]
    
    return X, y, reg_type




already_exist_pdb = listdir(map_dir_pdb)
def process_com(com_path_in, com, report_dict):
    if com in already_exist_pdb:
        return
    
    mkdir(path.join(map_dir_pdb, com))
    pos_path = path.join(map_dir_pdb, com, '1')
    mkdir(pos_path)
    
    mkdir(path.join(inter_dir_pdb, com))
    pos_path_inter = path.join(inter_dir_pdb, com, '1')
    mkdir(pos_path_inter)
    
    int_files = glob.glob(path.join(com_path_in, '*rec_0_dockinter.txt'))
    int_files = [x for x in int_files if stat(x).st_size > 300 and stat(x.replace('rec_0_dockinter.txt', 'lig_0_dockinter.txt')).st_size > 300]

    #already_chosen_pairs = []
    pdb_chain_pairs = []
    pairs = set([(path.basename(x).split('-')[0],
                   path.basename(x).split('-')[1].replace('_rec_0_dockinter.txt', '')) 
                   for x in int_files])
    for pair in pairs:
        #print(pair)
        #ch1 = pair[0].split('_')[0] + ' ' + pair[0].split('_')[3]
        #ch2 = pair[1].split('_')[0] + ' ' + pair[1].split('_')[3]
        #print(ch1, ch2)
        #try:
        #    cls1 = SCOPe.loc[SCOPe[4].str.contains(ch1)][2].iloc[0]
        #    cls2 = SCOPe.loc[SCOPe[4].str.contains(ch2)][2].iloc[0]
        #except:
        #    continue
        #print('classes', cls1, cls2)
        #if (pair[0][-3:], pair[1][-3:]) in already_chosen_pairs or (pair[1][-3:], pair[0][-3:]) in already_chosen_pairs:
        #    continue
        #if (cls1, cls2) in already_chosen_classes or (cls2, cls1) in already_chosen_classes:
        #    continue
        #already_chosen_pairs.append((pair[0][-3:], pair[1][-3:]))
        #already_chosen_classes.append((cls1, cls2))
        pdb_chain_pairs.append((pair[0], pair[1]))
        
    comp_list = []
    for pair in pdb_chain_pairs:
        try:
            interface_file_rec = path.join(com_path_in, pair[0] + '-' + pair[1] + '_rec_0_dockinter.txt')
            interface_file_lig = path.join(com_path_in, pair[0] + '-' + pair[1] + '_lig_0_dockinter.txt')
            interface_comp_rec = tl.read_interface_file(interface_file_rec, "", False, False, "all")[0]
            interface_comp_lig = tl.read_interface_file(interface_file_lig, "", False, False, "all")[0]
            
            rec_pdb = path.join(com_path_in, pair[0] + '.pdb')
            lig_pdb = path.join(com_path_in, pair[1] + '.pdb')
            print(rec_pdb)
            print(lig_pdb)
            
            
            sam_tup = (rec_pdb, lig_pdb, interface_comp_rec, interface_comp_lig, 1)
            mapcomplex(sam_tup, pair[0]+'--'+pair[1], pos_path)
            save_pdb_interface(sam_tup, pair[0]+'--'+pair[1], pos_path_inter)
            comp_list.append(sam_tup)
            
        except Exception as e:
            logging.error("Bad interface!" + '\ninter_rec: ' +  interface_file_rec + '\ninter_lig: ' + interface_file_lig + '\nError message: ' + str(e) + 
                          "\nMore information:\n" + traceback.format_exc())
    #report_dict[path.basename(com_path_in)] = comp_list

def manage_pdb_files(use_multiprocessing):
    complexes = listdir(pdbind_preprocessed_path)
    com_cases = []
    for com in complexes:
        com_path_in = path.join(pdbind_preprocessed_path, com)
        if not path.isdir(com_path_in):
            continue
        com_cases.append((com_path_in, com))
    report_dict = tl.do_processing(com_cases, process_com, use_multiprocessing)
    return report_dict

report_dict = manage_pdb_files(False)