#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:16:06 2021

@author: yasser
"""

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
import subprocess
from subprocess import CalledProcessError, check_call
import traceback
from random import shuffle, random, seed, sample
from numpy import newaxis
import matplotlib.pyplot as plt
import time
import load_data as load

import collections
import scr
from numpy import asarray
from sklearn.preprocessing import OneHotEncoder

sys.path.insert(1, '../lib/')
import tools as tl

logging.basicConfig(filename='manager.log', filemode='w', format='%(levelname)s: %(message)s', level=logging.DEBUG)
mainlog = logging.getLogger('main')
logging.Logger

map_dir_mut = 'map_dir_mut'
map_dir_mut_sep = 'map_dir_mut_sep'
inter_dir_mut = 'inter_dir_mut'

bin_path = "./maps_generator"

v_dim = 24

encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
onehot = encoder.fit(asarray([['CYS'], ['ASP'], ['SER'], ['GLN'], ['LYS'], ['ILE'], ['PRO'], 
                              ['THR'], ['PHE'], ['ASN'], ['GLY'], ['HIS'], ['LEU'], ['ARG'], 
                              ['TRP'], ['ALA'], ['VAL'], ['GLU'], ['TYR'], ['MET']]))

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

# three_letter["S"] will now return "SER"
three_letter = dict([[v,k] for k,v in one_letter.items()])

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS'}


def save_pdb_interface(pdbfile, map_name, interface):
    try:
        inter_dict = {}
        for inter in interface:
            inter_dict.setdefault(inter[0], [])
            inter_dict[inter[0]].append(str(inter[1]))
    
        pdbdata = parsePDB(pdbfile).select('protein')
        
        
        resnum = 'chain ' + list(inter_dict.keys())[0] + ' resnum ' + " ".join(list(inter_dict.values())[0])
        atoms = pdbdata.select(resnum).toAtomGroup()
        for k,v in list(inter_dict.items())[1:]:
            resnum = 'chain ' + k + ' resnum ' + " ".join(v)
            atoms += pdbdata.select(resnum).toAtomGroup()
        
        writePDB(map_name + '.pdb', atoms)
    except:
        return


def mapcomplex(wt_pdb, mut_pdb, protein_complex, mutate_complex, ddg, model_ind, mut_name, interaction_region, complex_type, experimental_method):    
    try:
        name = mut_name + '--' + str(model_ind)
                
        mapcommand = [bin_path, "--mode", "map", "-i", wt_pdb, "--native", "-m", str(v_dim), "-t", "167", "-v", "0.8", "-o", 'wt.bin']
        subprocess.call(mapcommand)

        mapcommand = [bin_path, "--mode", "map", "-i", mut_pdb, "--native", "-m", str(v_dim), "-t", "167", "-v", "0.8", "-o", 'mut.bin']
        subprocess.call(mapcommand)
        
        dataset_trai_wt = load.read_data_set('wt.bin')
        dataset_trai_mut = load.read_data_set('mut.bin')
        
        X_wt = np.reshape(dataset_trai_wt.maps, (-1,v_dim,v_dim,v_dim,167+6))
        X_mut = np.reshape(dataset_trai_mut.maps, (-1,v_dim,v_dim,v_dim,167+6))
        y = ddg
        map_name = path.join(map_dir_mut, mut_name, name)
        tl.save_obj((X_wt, X_mut, y, interaction_region, complex_type, experimental_method), map_name)
        
        check_call(
            [
                'lz4', '-9', '-f',   #, '--rm' because if inconsistency in lz4 versions! 
                map_name + '.pkl'
            ],
            stdout=sys.stdout)
        
        remove(map_name + '.pkl')
        remove('wt.bin')
        remove('mut.bin')
        #remove(wt_pdb.replace('pdb','bin'))
        #remove(mut_pdb.replace('pdb','bin'))
        #remove(name+'_complex_p1.pdb')
        #remove(name+'_complex_p2.pdb')
        #remove(name+'_complex.pdb')
        #remove(name+'_rimcoresup.csv')
        
        
        ##########################################################################
        y_hotencode = encoder.transform([[three_letter[mutate_complex[0]]]])
        map_name = path.join(map_dir_mut_sep, mut_name, '1', 'wt_' + name)
        tl.save_obj((X_wt,y_hotencode,y, interaction_region, complex_type, experimental_method), map_name)
        check_call(
            [
                'lz4', '-9', '-f',   #, '--rm' because if inconsistency in lz4 versions! 
                map_name + '.pkl'
            ],
            stdout=sys.stdout)
        remove(map_name + '.pkl')
        
        y_hotencode = encoder.transform([[three_letter[mutate_complex[-1]]]])
        map_name = path.join(map_dir_mut_sep, mut_name, '1', 'mut_' + name)
        tl.save_obj((X_mut,y_hotencode,y, interaction_region, complex_type, experimental_method), map_name)
        check_call(
            [
                'lz4', '-9', '-f',   #, '--rm' because if inconsistency in lz4 versions! 
                map_name + '.pkl'
            ],
            stdout=sys.stdout)
        remove(map_name + '.pkl')
        ##########################################################################
        
    except Exception as e:
        logging.error("Bad interface!" + '\ninter_rec: ' +  wt_pdb + '\ninter_lig: ' + mut_pdb + '\nError message: ' + str(e) + 
                      "\nMore information:\n" + traceback.format_exc())
        print('XXXXXXXXXXXXXXexceptionXXXXXXXXXx')
        return
    
    return


#already_exist_mutation = listdir(map_dir_mut)
already_exist_mutation = []

def Kd2DeltaG(Kd):
    return (8.314/4184)*(273.15+25)*np.log(Kd)


#Mutations
def process_mut(mut, mut_directory, report_dict):
    try:
        if mut in already_exist_mutation:
            return
        print(mut)
                
        #Read the description of the sample
        with open(path.join(mut_directory, mut, "description")) as infile:
            sample_desc = infile.readline()
        protein_complex, mutate_complex, interaction_region, complex_type, experimental_method, b_affine_wt, b_affine_mt = sample_desc.split('§')
        b_affine_wt = Kd2DeltaG(float(b_affine_wt.replace("<","").replace(">","")))
        b_affine_mt = Kd2DeltaG(float(b_affine_mt.replace("<","").replace(">","")))
        
        mutate_info = mutate_complex.split(',')
        f_handler = open("resinfo", 'w')
        for mi in mutate_info:
            f_handler.write(mi[2:-1]+';'+mi[1]+'\n')
        f_handler.close()
        
        #Skip multiple single point mutations for now!
        if len(mutate_info) > 1:
            return

        mut_name = mut.replace('§','--')
        mkdir(path.join(map_dir_mut, mut_name))
        mkdir(path.join(map_dir_mut_sep, mut_name))
        mkdir(path.join(map_dir_mut_sep, mut_name, '1'))
        mkdir(path.join(inter_dir_mut, mut_name))        
        
        for i in range(1,31):
            wt_pdb = path.join(mut_directory, mut, str(i).zfill(2), "wt_"+str(i).zfill(2)+"_10000.pdb")
            mut_pdb = path.join(mut_directory, mut, str(i).zfill(2), "mut_"+str(i).zfill(2)+"_10000.pdb")
            mapcomplex(wt_pdb, mut_pdb, protein_complex, mutate_complex, b_affine_mt - b_affine_wt, i, mut_name, interaction_region, complex_type, experimental_method)
    except:
        return
        
def manage_mut_files(use_multiprocessing): 
    if not path.exists(map_dir_mut):
        mkdir(map_dir_mut)
    if not path.exists(map_dir_mut_sep):
        mkdir(map_dir_mut_sep)
    if not path.exists(inter_dir_mut):
        mkdir(inter_dir_mut)
    mut_directory = 'backrub_models_directory'
    muts = listdir(mut_directory)
    mut_cases = []
    for mut in muts:
        mut_cases.append((mut, mut_directory))
    report_dict = tl.do_processing(mut_cases, process_mut, use_multiprocessing)
    return report_dict


report_dict = manage_mut_files(False)
