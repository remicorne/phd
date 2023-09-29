#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 14:56:14 2022

@author: shiva
"""
#%%

import os
from typing import Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from DLC_module import (create_empty_df_summaries, 
                        read_DLC_file,
                        check_video_length_long_enough,
                        get_cage_coords,
                        normalize_cage_coords,
                        pre_process_df,
                        cal_time_spent_percentage,
                        calculate_distance_traveled,
                        fill_df_summaries,
                        drop_extra_rows_df,
                        list_all_data_files,
                        create_dict_of_file_names,
                        grubbs_test_DLC,
                        save_box_plot_quadrent_data,
                        save_table_of_quadrent_time,
                        save_dist_traveled_box_plot,
                        save_table_of_distance_traveled,
                        produce_heat_map_all_groups,
                        produce_trajectory_all_groups
                        ) ## must be in the same directory


# %% Loop File Function


### this function should be monadic
def loop_over_files(column_names : list, column_names_sns : list, n_files : int, 
                    filename_dict : dict, folder_path:str, n_min_to_keep:int, 
                    side_length_cm : int, fps : int = 30):
    
    df_summary, df_summary_sns  = create_empty_df_summaries(column_names, column_names_sns, n_files)
    
    count : int = 0
    
    for group, filename_list in filename_dict.items():

        n_files_in_group = len(filename_list)
        
        df_one_group : pd.DataFrame = pd.DataFrame( np.zeros(( int(n_files_in_group * n_pts_to_keep), 3)), columns = ['x', 'y', 'mouse_no'])
        last_time_pts : int = 0
        
        for n_f, filename in enumerate(filename_list):

            mouse_no : str =  filename[1:3]
            
            filepath : str = os.path.join(folder_path, filename)
            
            df : pd.DataFrame = read_DLC_file(filepath)
            
            video_not_long_enough : bool = check_video_length_long_enough(df, threshold = n_min_to_keep, fps = fps)
            
            if video_not_long_enough: # skip if video not long enough
                print(' insufficient length')
                continue

            cage_coord : dict = get_cage_coords(df, coord_list = [ 'x', 'y'], 
                                label_list = ['Center', 'ULCorner',  'URCorner',  
                                              'LLCorner',  'LRCorner'],
                                p_cutoff = 0.95)



            cage_coord_normalized, scale_pix_to_cm = normalize_cage_coords(cage_coord, side_length_cm)


            
            df_filtered,  n_pts_lost, longest_gap, cage_coord_normalized = pre_process_df(df, n_pts_to_keep, cage_coord_normalized, 
                                                                                          scale_pix_to_cm, p_cutoff = 0.8, 
                                                                                          plot_missing_dist = False)
            

            time_spent_points, time_spent_percentage = cal_time_spent_percentage(df_filtered)

            distance_traveled = calculate_distance_traveled (df_filtered, cage_coord_normalized)
            
            df_summary, df_summary_sns = fill_df_summaries(df_summary, df_summary_sns, mouse_no, group, time_spent_percentage, 
                                                            n_pts_lost, longest_gap, distance_traveled, count)
            
            
            df_one_group[ last_time_pts: 
                          last_time_pts +  n_pts_to_keep - n_pts_lost] = pd.Series({'x': df_filtered['CM', 'x'] ,
                                                                                      'y': df_filtered['CM', 'y'], 
                                                                                     'mouse_no': [mouse_no] * 
                                                                                     (n_pts_to_keep - n_pts_lost)})
                                                                                    
            
            last_time_pts += n_pts_to_keep
            
            count += 1
            
            
         
            
        df_one_group : pd.DataFrame = df_one_group[ df_one_group['x'] != 0 ]   
        df_one_group.to_csv(os.path.join(folder_path, 'CM_analysis', group + '_CM_x_y.csv')) # saving CM data from heatmapping later 
            
            
    
    print('number of videos analyzed = ', count)
    
    df_summary : pd.DataFrame = drop_extra_rows_df(df_summary, col_name_to_check_if_filled = 'mouse_no')
    df_summary_sns = drop_extra_rows_df(df_summary_sns, col_name_to_check_if_filled = 'mouse_no')
    
    return df_summary, df_summary_sns, cage_coord_normalized, df_filtered




if __name__ == "__main__":

    #%% VAIRABLES 
    fps : int = 30
    n_min_to_keep : int  = 5 # minutes to keep from video
    n_pts_to_keep : int = n_min_to_keep * fps * 60
    side_length_cm : int = 33
    folder_path : str = './input/'
    


    all_data_files_list : list = list_all_data_files(folder_path, extensions = ['csv'])
    n_files : int = len(all_data_files_list)
    filename_dict : dict = create_dict_of_file_names(all_data_files_list)   



    column_names : list = ["mouse_no", 'group', "quad_0", "quad_1", "quad_2", "quad_3", 'n_pts_lost', 'longest_gap', 'distance_traveled']
    column_names_sns : list = ["mouse_no", 'group', "quad", "time"]


    #%% RUN

    '''
    (needs normalsation between treatment groups for lost data due to DLC poor detection)

    '''


    df_summary, df_summary_sns, cage_coord_normalized, df_filtered = loop_over_files(column_names, column_names_sns, n_files, 
                                                filename_dict, folder_path,
                                                n_min_to_keep, side_length_cm, 
                                                fps = fps)   



    treatment_dict : dict = { 1: 'vehicles', 2 : '10mg/kgTCB', 3 : '3mg/kgTCB', 4 : '0.3mg/kgTCB',5 : 'TCB+MDL', 6 : '0.2mg/kgMDL' }     #TCB2



    palette_labeled : dict ={'vehicles': "white",      #TCB2
                    '10mg/kgTCB': "firebrick", 
                    '3mg/kgTCB': "red", 
                    '0.3mg/kgTCB': "salmon",
                    'TCB+MDL': "grey",
                    '0.2mg/kgMDL': "black"}


    list_of_groups : np.ndarray = np.unique( df_summary['group'])

    # outlier detection

    df_summary_raw : pd.DataFrame = df_summary.copy()
    df_summary : pd.DataFrame = grubbs_test_DLC(df_summary_raw, treatment_dict, list_of_groups, p_value = 0.05, remove_all = False)


    # Producing 

    '''
    Quadrent Analysis: 
        save_box_plot_quadrent_data - generates a box and wisker plot for number of data points in each quadrent, 
                                        for each treatent group followed by a dunns test 
                                        
        save_table_of_quadrent_time - generates table of time spent in each quadrent for each treatent group including SEM
    '''
    save_box_plot_quadrent_data(df_summary_sns)  #OUTLIERS ARE NOT REMOVED FROM df_summary_sns
    save_table_of_quadrent_time(df_summary)

    '''
    Distance Traveled: 
            save_dist_traveled_box_plot - creates a box plot of the distance traveled for each treatent group 
                                        (plot 1 : all 6 groups, plot 2: dose responce, plot 3: agonist antagonist)
                                        
            save_table_of_distance_traveled - enerates table of distance traveled for each treatent group
    '''

    save_dist_traveled_box_plot(df_summary) # also saves histogram 

    save_table_of_distance_traveled(df_summary, plot_hist = True)



    '''
    Heat Mapping: 
            using previously save CM (x,y) values to create and save heat maps for  each group and trajectory plots 
    '''

    X_Y_df_file_list : list = list_all_data_files(os.path.join(folder_path, 'CM_analysis')) #creates list of files to be accessed for CM data 

    plt.close('all')


    produce_heat_map_all_groups(X_Y_df_file_list, treatment_dict, folder_path, name = '5_min_')

    produce_trajectory_all_groups(X_Y_df_file_list, cage_coord_normalized, folder_path, name = '_5_min_')


#%%## run all funcs for one file

# folder_path = '/Users/jasminebutler/Desktop/DLC_Analysis_org'
# all_data_files_list = list_all_data_files(folder_path, extensions = ['csv'])
# n_files = len(all_data_files_list)
# filename_dict = create_dict_of_file_names(all_data_files_list)   

# column_names = ["mouse_no", 'group', "quad_0", "quad_1", "quad_2", "quad_3", 'n_pts_lost', 'longest_gap']
# df_summary = pd.DataFrame(np.zeros((n_files, len(column_names))), columns = column_names)
# filename = 'M68_2.csv'
# mouse_no =  filename[1:3]
# group = filename[4:5]
# # i use filename later to keep track of each mouse refering filename[0:5] therefore need it to be defined in the loop over files
# filepath = os.path.join(folder_path, filename)

# df = read_DLC_file(filepath)
# cage_coord = get_cage_coords(df, coord_list = [ 'x', 'y'], 
#                     label_list = ['Center', 'ULCorner',  'URCorner',  
#                                   'LLCorner',  'LRCorner'],
#                     p_cutoff = 0.8)
# fps = 30
# n_min_to_keep = 3 # minutes to keep from video
# n_pts_to_keep = n_min_to_keep * fps * 60
# df_trimmed = trim_df(df, n_pts_to_keep)
# df_filtered, n_pts_lost, longest_gap = filter_body_parts_based_on_p_cutoff(df_trimmed, p_cutoff = 0.7, 
#                                                                                 plot_missing_dist =True)

# ##### normalize here using cage coords 

# df_filtered = calculate_and_add_mouse_CM(df_filtered, plot_missing_time_pts = True, 
#                                           plot_trajectory = False)

# df_filtered = find_which_quad(df_filtered, cage_coord)

# time_spent = time_spent_in_each_quad (df_filtered)
# df_summary = fill_df_of_each_group_analysis(mouse_no, group, time_spent, 
#                                             n_pts_lost, longest_gap, i = 0)



# %% imobility ??????

''' 
calculate speed < threshold for a period of time?
'''

