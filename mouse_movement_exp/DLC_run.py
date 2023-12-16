#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 14:56:14 2022

@author: shiva
"""
#%%

import os
from typing import Tuple
import pickle

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from DLC_module import (
                        read_DLC_file,
                        check_video_length_long_enough,
                        get_cage_coord_full,
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
                        produce_trajectory_all_groups,
                        trim_df,
                        filter_df,
                        frames_missing_consec,
                        plot_missing_hist,
                        scale_df_vals_to_cm,
                        plot_time_points,
                        calculate_and_add_mouse_CM,
                        plot_trajectory,
                        change_coord_rel_to_cage,
                        constrain_to_cage_coord,
                        find_which_quad
                        ) ## must be in the same directory


# %% Loop File Function


### this function should be monadic
def loop_over_files(column_names : list, column_names_sns : list, n_files : int, 
                    filename_dict : dict, folder_path:str, n_min_to_keep:int, 
                    side_length_cm : int, fps : int = 30):
    
    df_summary : pd.DataFrame = pd.DataFrame(np.zeros((n_files, len(column_names))), columns = column_names)
    df_summary_sns : pd.DataFrame  = pd.DataFrame(np.zeros(( int( n_files * 4) , 4)), columns = column_names_sns)


    count : int = 0
    
    for group, filename_list in filename_dict.items():

        n_files_in_group : int = len(filename_list)
        last_time_pts : int = 0
        df_one_group : pd.DataFrame = pd.DataFrame( np.zeros(( int(n_files_in_group * n_pts_to_keep), 3)), columns = ['x', 'y', 'mouse_no'])
        
        
        for filename in filename_list:

            mouse_no : str =  filename[1:3]
            
            filepath : str = os.path.join(folder_path, filename)
            
            df : pd.DataFrame = read_DLC_file(filepath)
            
            video_not_long_enough : bool = check_video_length_long_enough(df, n_min_to_keep, fps)
            
            if video_not_long_enough: 
                print(' insufficient length')
                continue


            cage_coord_full : dict = get_cage_coord_full(df)


            df_trimmed : pd.DataFrame = trim_df(df, n_pts_to_keep)

            
            def get_high_p_ind(df : pd.DataFrame) -> tuple:
                p_cutoff : float = 0.8
                return  (( df['Nose']['likelihood'] > p_cutoff) & 
                        ( df['Tail']['likelihood'] > p_cutoff) &
                        ( df['LEar']['likelihood'] > p_cutoff) &
                        ( df['REar']['likelihood'] > p_cutoff) )


            high_p_ind : tuple = get_high_p_ind(df_trimmed)

            df_filtered : pd.DataFrame = filter_df(df_trimmed,high_p_ind)

            df_filtered.longest_gap : float = frames_missing_consec( np.where( high_p_ind ))
            df_filtered.n_pts_lost : int = len(df_trimmed) - len(df_filtered)

            # the last header is 'which' we don't want that -> your missing columns is bodyparts 

            df_filtered  = scale_df_vals_to_cm(df_filtered, cage_coord_full['scale_pix_to_cm'], coord_list = ['x', 'y'])
            
            df_filtered = calculate_and_add_mouse_CM(df_filtered)

            df_filtered = change_coord_rel_to_cage(df_filtered, cage_coord_full['zero'])
    
            df_filtered = constrain_to_cage_coord(df_filtered, cage_coord_full['cage_coord_normalized']) ############# BEWARE yo won't see the miss detections anymorte with this

            df_filtered = find_which_quad(df_filtered, cage_coord_full['cage_coord_normalized']) ### some reoccuring warning try df_filtered[][] = [quad ] * len(df_filtered[][]) to fix SHIVA




            time_spent_points, time_spent_percentage = cal_time_spent_percentage(df_filtered)

            distance_traveled = calculate_distance_traveled (df_filtered, cage_coord_full['cage_coord_normalized'])

            plot_missing_time_pts : bool = False
            if plot_missing_time_pts:
                plot_time_points(df_filtered)

            is_plot_trajectory: bool = False
            if is_plot_trajectory:
                plot_trajectory(df_filtered)

            plot_missing_dist: bool = False
            if plot_missing_dist:
                total_length_of_video : int = len(df.index)
                bin_length: int = 30
                plot_missing_hist(high_p_ind, int( total_length_of_video / bin_length))

            
            df_summary, df_summary_sns = fill_df_summaries(df_summary, df_summary_sns, mouse_no, group, time_spent_percentage, 
                                                            df_filtered.n_pts_lost, df_filtered.longest_gap, distance_traveled, count)
            
            
            df_one_group[ last_time_pts: 
                          last_time_pts +  n_pts_to_keep - df_filtered.n_pts_lost] = pd.Series({'x': df_filtered['CM', 'x'] ,
                                                                                      'y': df_filtered['CM', 'y'], 
                                                                                     'mouse_no': [mouse_no] * 
                                                                                     (n_pts_to_keep - df_filtered.n_pts_lost)})
                                                                                    
            
            last_time_pts += n_pts_to_keep
            
            count += 1
            
            
         
            
        df_one_group : pd.DataFrame = df_one_group[ df_one_group['x'] != 0 ]   
        df_one_group.to_csv(os.path.join(folder_path, 'CM_analysis', group + '_CM_x_y.csv')) # saving CM data from heatmapping later 
            
            
    
    print('number of videos analyzed = ', count)
    
    df_summary : pd.DataFrame = drop_extra_rows_df(df_summary, col_name_to_check_if_filled = 'mouse_no')
    df_summary_sns = drop_extra_rows_df(df_summary_sns, col_name_to_check_if_filled = 'mouse_no')
    
    return df_summary, df_summary_sns, cage_coord_full['cage_coord_normalized'], df_filtered




if __name__ == "__main__":

    #%% VAIRABLES 
    fps : int = 30
    n_min_to_keep : int  = 5 # minutes to keep from video
    n_pts_to_keep : int = n_min_to_keep * fps * 60
    side_length_cm : float = 33
    folder_path : str = './mouse_movement_exp/input/'
    


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


    ### quality test
    ### pickle dump X_Y_df_file_list
    with open(os.path.join(".\\mouse_movement_exp\\test","X_Y_df_file_list"),"wb+") as fid:
        fid.write(pickle.dumps(X_Y_df_file_list))
    ### df_summary
    with open(os.path.join(".\\mouse_movement_exp\\test","df_summary"),"wb+") as fid:
        fid.write(pickle.dumps(df_summary))
    ### df_summary_sns
    with open(os.path.join(".\\mouse_movement_exp\\test","df_summary_sns"),"wb+") as fid:
        fid.write(pickle.dumps(df_summary_sns))

    ### data are the same between two run
    # with open(os.path.join(".\\mouse_movement_exp\\test","X_Y_df_file_list"),"rb") as fid:
    #     result2assert = fid.read()
    #     print(pickle.dumps(X_Y_df_file_list) == result2assert)

    # ### data are not the same between two run
    # with open(os.path.join(".\\mouse_movement_exp\\test","df_summary"),"rb") as fid:
    #     result2assert = fid.read()
    #     print(pickle.dumps(df_summary) == result2assert)

    # ### data are not the same between two run
    # with open(os.path.join(".\\mouse_movement_exp\\test","df_summary_sns"),"rb") as fid:
    #     result2assert = fid.read()
    #     print(pickle.dumps(df_summary_sns) == result2assert)