#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 12:39:18 2022

@author: shiva
"""

import os

def list_all_data_files(path, extensions = ['csv']):
    
    '''get all the files with extention in the path where you want to search'''
    
    files = [x for x in os.listdir(path) if not x.startswith('.')]
    files.sort()
    videofile_path = [ fi for fi in files if fi.endswith( tuple(extensions)) ]
    
    return videofile_path


def key_return(X, dictionary):
    for key, value in dictionary.items():
        if X == value:
            return key
        if isinstance(value, list) and X in value:
            return key
    return "Key doesnt exist"


def add_group_to_name(DLC_fresh_filepath_list, mouse_group_dict, path_DLC_fresh):
    for filename in DLC_fresh_filepath_list:
        
        mouse_no = filename.split('_')[0] [1:] # first charecter is M, the rest is number
        corr_group = key_return( int(mouse_no) , 
                                mouse_group_dict)
        new_filename = 'M' + mouse_no + '_' + corr_group + '.csv' # add the corresponding group to the mouse number to generate new name
        src_path = os.path.join(path_DLC_fresh, filename)
        des_path = os.path.join(path_DLC_fresh, new_filename)
        os.rename(src_path, des_path)
    
path_DLC_fresh = '/Users/jasminebutler/Desktop/DLC_Analysis_org'

DLC_fresh_filepath_list = list_all_data_files(path_DLC_fresh, extensions = ['csv'])

mouse_group_dict = {'1' :[5,7,9,11,17,20,28,32,59,67], '2': [8,10,12,47,54,56,60,62,64,68,70], '3': [16,18,19,22,27,30,37,39,46,53,55], '4': [13,14,15,26,29,34,42,48,63,65], '5': [21,25,35,41,45,49,52,58,61,66,69],'6': [23,24,36,38,40,44,50,51,57]} # example, mouse numbers corresponding to each group

add_group_to_name(DLC_fresh_filepath_list, mouse_group_dict, path_DLC_fresh)

