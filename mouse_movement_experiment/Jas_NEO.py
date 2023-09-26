#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 15:28:51 2022

@author: shiva
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import csv
import subprocess

# implement pip as a subprocess:
try:
    import neo
except ImportError or ModuleNotFoundError:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'neo'])
try:
    import sonpy
except ImportError or ModuleNotFoundError:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'sonpy'])
try:
    import scipy
except ImportError or ModuleNotFoundError:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'scipy'])


def read_neo_file_return_analogsignals(filename, neo_obj):
    
    block = neo_obj.read()[0] # read the file 
    analogsignals = block.segments[0].analogsignals
    report_info_on_file(filename, block, analogsignals)
    
    return analogsignals

def report_info_on_file(filename, block, analogsignals):
    
    print("file = ", filename)
    print('number of segments = ', len(block.segments))
    print('number of analog signals = ', len(analogsignals))
    
    for i in range(len(analogsignals)):
        print('signal {} contains series of shape {} '.format(i+1, analogsignals[i].shape))
     
def enforce_neo_version_satisfaction():
    v1, v2, v3 = neo_version.split('.')
    if v1 == 0 and v2 < 10:
        path_to_package = input("neo version >= 0.10.0 is required in order to read .smrx files. \n \
                                Please download the latest release from : \
                                https://github.com/NeuralEnsemble/python-neo/releases/tag/0.10.0, \n \
                                Then, extract the zip and input the absolute path \
                                to the extracted folder here(e.g. \n /home/User-name/Downloads/python-neo-0.10.0): \n")
        subprocess.check_call([sys.executable, '-m', 'pip', 'install',
                               path_to_package])
        
def read_laser_file(filepath, neo_version_check):
    
    filename = os.path.basename(filepath)
    file_extension = filename.split('.') [-1]
    
    ## Depending on the file extension, the laser information is stored differently
    if file_extension == 'smrx':
        if not neo_version_check:
            enforce_neo_version_satisfaction()
            neo_version_check = True
        neo_obj = neo.CedIO(filepath)
        analogsignals = read_neo_file_return_analogsignals(filename, neo_obj)
        
        # keep the signal as a 16 bit float
        laser_series = np.float16(analogsignals[1]) # the laser information is stored as the second analog signal
    else:
        neo_obj = neo.Spike2IO(filepath)
        analogsignals = read_neo_file_return_analogsignals(filename, neo_obj)
        
        # keep the signal as a 16 bit float
        laser_series = np.float16(analogsignals[1][:,1]) # the laser information is stored as the second column in the 
                                             # the second analog signal
        
    return laser_series  



if __name__ == '__main__':
    
    global neo_version_check
    neo_version_check = False
    neo_version = neo.__version__
    print("Neo package version is : {}".format(neo_version) )
    filepath = input("Enter the full path under which the .smr/.smrx file hierarchy resides:  \n")     

    read_laser_file(filepath, neo_version_check)
