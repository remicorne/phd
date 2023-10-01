#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 13:48:56 2022

@author: jasminebutler
"""

'''
DLC_module

Functions for the analysis of mouse video data 
'''
#%% INSTALATIONS


import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
import statistics
import seaborn as sns
from scipy import stats
import scikit_posthocs as sp
from outliers import smirnov_grubbs as grubbs



#%% FUNCTIONS


def list_all_data_files(path : str, extensions : list = ['csv']) -> list:
    
    '''get all the files with extention in the path where you want to search'''
    
    files : list = [x for x in os.listdir(path) if not x.startswith('.') and x.endswith( tuple(extensions))]
    files.sort()

    return files

def create_dict_of_file_names(files : list) -> dict:
    
    '''
    create a dictionary of filenames with keys corresponding to the group number
    '''

    def get_group_from_filename(filename):
        return filename[4:5]

    groups : set = set( [ get_group_from_filename(val) for val in files])
    
    filenames_groupby_group = {group: [filename for filename in files if get_group_from_filename(filename) == group] 
                     for group in groups}
    
    return filenames_groupby_group

def avg_high_likelihoods(df : pd.DataFrame, label_coord_pair : tuple, p_cutoff : float = 0.8) -> pd.DataFrame:
    
    ''' Return the average of the <coord> coordinate of the <label> trackings 
        having detection likelihood above a certain p_cutoff'''
    label : str = label_coord_pair[0] 
    coord : str = label_coord_pair[1] 

    ind : pd.Series = df[label]['likelihood'] > p_cutoff 

    return df[label][coord][ind].mean()



def get_cage_coords(df : pd.DataFrame, coord_list : list = [ 'x', 'y'], 
                    label_list : list = ['Center', 'ULCorner',  'URCorner',  
                                  'LLCorner',  'LRCorner'],
                    p_cutoff : float = 0.8) -> dict:
    ''' Return a dictionary of the high likelihood averaged coordinates of cage corners and center'''
    

    coord_label_pair : list = [(label,coord) for label in label_list for coord in coord_list]
    
    return {pair : avg_high_likelihoods(df, pair, p_cutoff) for pair in coord_label_pair}

def read_DLC_file(filepath : str) -> pd.DataFrame:
    
    ''' Read the DLC file, discard the first Row and consider two lines as headers'''
    
    return pd.read_csv(filepath, header = [0,1], skiprows= [0])



def filter_df (df : pd.DataFrame, high_p_ind:tuple) -> pd.DataFrame:

    df_filtered : pd.DataFrame = df[ high_p_ind ]
    return df_filtered.drop( columns = df_filtered.columns[1:10], axis = 1) # drop info on cage



def frames_missing_consec(chosen_inds : np.ndarray , fps : int = 30) -> float:
    
    '''
        rerturn the longest gap (in s) produced by filter in likliohood
    '''
    
    frame_diffs = np.diff( chosen_inds)
    return np.max(frame_diffs) / fps
    
def plot_missing_hist(high_p_ind, bins, fps : int = 30):
    '''
    plot histogram in seconds

    '''
    fig, ax = plt.subplots()
    # bins here  is the number of frames for each bin 1sec = 30
    ax.hist( np.where(high_p_ind)[0] / fps, bins = bins)

def plot_trajectory(df_filtered : pd.DataFrame):
        
    fig, ax = plt.subplots()
    ax.plot(df_filtered['CM','x'], df_filtered['CM','y'], '-o', markersize= 7)
    
def plot_time_points(df_filtered : pd.DataFrame):
    
    ''' to see discrepencies - location of removed data '''
    
    coords = df_filtered['bodyparts', 'coords']
    fig, ax = plt.subplots()
    ax.plot(coords, np.full_like(coords, 1), 'o', markersize= 1)

def trim_df(df : pd.DataFrame, n_pts_to_keep : int) -> pd.DataFrame:
    
    return df.drop( np.arange(n_pts_to_keep, len(df.index)))

##### Find CM
def calculate_and_add_mouse_CM(df_filtered : pd.DataFrame) -> pd.DataFrame:
        
    body_parts : list = ['Nose', 'Tail', 'REar', 'LEar']    
    x_columns : list = [[value, 'x']  for value in body_parts]
    y_columns : list = [[value, 'y']  for value in body_parts]
    
    df_filtered['CM','x'] = df_filtered[x_columns].mean(axis = 1) # filling info on CM
    df_filtered['CM','y'] = df_filtered[y_columns].mean(axis = 1)
    
    # plot CM x,y points - note: floor  90 degrees rotation anticlockwise 
    return df_filtered





def find_mid_lines(cage_coord):
    
    cage_coord['UM', 'x'] = statistics.mean([cage_coord['ULCorner', 'x'], 
                                       cage_coord['URCorner', 'x']])
    cage_coord['UM', 'y'] = statistics.mean([cage_coord['ULCorner', 'y'], 
                                       cage_coord['URCorner', 'y']])
    cage_coord['LeM', 'x'] = statistics.mean([cage_coord['ULCorner', 'x'], 
                                       cage_coord['LLCorner', 'x']])
    cage_coord['LeM', 'y'] = statistics.mean([cage_coord['ULCorner', 'y'], 
                                       cage_coord['LLCorner', 'y']])
    cage_coord['LM', 'x'] = statistics.mean([cage_coord['LRCorner', 'x'], 
                                       cage_coord['LLCorner', 'x']])
    cage_coord['LM', 'y'] = statistics.mean([cage_coord['LRCorner', 'y'], 
                                       cage_coord['LLCorner', 'y']])
    cage_coord['RM', 'x'] = statistics.mean([cage_coord['LRCorner', 'x'], 
                                       cage_coord['URCorner', 'x']])
    cage_coord['RM', 'y'] = statistics.mean([cage_coord['LRCorner', 'y'], 
                                       cage_coord['URCorner', 'y']])
    return cage_coord



def find_which_quad(df_filtered : pd.DataFrame, cage_coord : dict) -> pd.DataFrame:
    
    '''
        finds the quadrent of the data point and appends a column [which][quad]
        0 = upper left = white
        1 = upper right = brown
        2 = lower left = orange
        3 = lower right = red 
        
    '''
    # print( 'n time points = ', len(df_filtered.index))
    cage_coord = find_mid_lines(cage_coord)
    df_filtered['which','quad'] = np.zeros(df_filtered.index.shape)
    
    ind_q_0 = ((df_filtered['CM','x'] < statistics.mean([cage_coord['UM', 'x'], 
                                                  cage_coord['Center', 'x']])) &
               (df_filtered['CM','y'] < statistics.mean([cage_coord['LeM', 'y'], 
                                                  cage_coord['Center', 'y']])))
    df_filtered['which','quad'][ind_q_0] = 0
    
    
    ind_q_1 = ((df_filtered['CM','x'] > statistics.mean([cage_coord['UM', 'x'], 
                                              cage_coord['Center', 'x']])) &
               (df_filtered['CM','y'] < statistics.mean([cage_coord['LeM', 'y'], 
                                              cage_coord['Center', 'y']])))
    df_filtered['which','quad'][ind_q_1] = 1
    
    ind_q_2 = ((df_filtered['CM','x'] < statistics.mean([cage_coord['LM', 'x'], 
                                              cage_coord['Center', 'x']])) &
               (df_filtered['CM','y'] > statistics.mean([cage_coord['LeM', 'y'], 
                                              cage_coord['Center', 'y']])))
    df_filtered['which','quad'][ind_q_2] = 2
    
    ind_q_3 = ((df_filtered['CM','x'] > statistics.mean([cage_coord['LM', 'x'], 
                                              cage_coord['Center', 'x']])) &
               (df_filtered['CM','y'] > statistics.mean([cage_coord['RM', 'y'], 
                                              cage_coord['Center', 'y']])))
    df_filtered['which','quad'][ind_q_3] = 3
    
    
    return df_filtered

# df_filtered = find_which_quad(df_filtered, cage_coord)

#### time spent in each quadrent

def time_spent_in_each_quad (df_filtered : pd.DataFrame):

    '''calculates the tpercentage of time spent in each quadrent from the filtered data'''
    quads = np.arange(4)
    time_spent = [ sum(df_filtered['which','quad'] == val) for val in quads]
    # print(time_spent)
    #adding mouse and grpup name to time_spent 
    # list_to_append = [filename[0:5]]  + time_spent      
    #plotting bar graph for each quadrent
    #fig, ax = plt.subplots()
    # Use automatic FuncFormatter creation
    #ax.bar(quads, time_spent)

    return time_spent#, list_to_append

def fill_df_of_each_group_analysis(df_summary, mouse_no, group, time_spent, n_pts_lost, longest_gap, distance_traveled, i = 0):
        
    df_summary.loc[i] = [mouse_no] + [group] + time_spent + [n_pts_lost] + [longest_gap] +[distance_traveled]

    
    return df_summary

# df_quadrents = create_df_of_time_spent_in_quad(df_filtered)

def fill_df_of_each_group_quad_concat(df_summary_sns, mouse_no, group, time_spent,i = 0):
        
    df_summary_sns.loc[i] = [mouse_no] + [group] + ['0']+ [time_spent[0]]
    df_summary_sns.loc[i+1] = [mouse_no] + [group] + ['1']+ [time_spent[1]] 
    df_summary_sns.loc[i+2] = [mouse_no] + [group] + ['2']+ [time_spent[2]] 
    df_summary_sns.loc[i+3] = [mouse_no] + [group] + ['3']+ [time_spent[3]] 

    
    return df_summary_sns

#### box plot for each treatment group

def save_box_plot_quadrent_data(df_summary_sns):

    fig, ax = plt.subplots()
    
    #convert noumber of points to percentage
    group_list = ["1", "2","3", "4","5", "6"]
    pdf = PdfPages(os.path.join("./output",'box_plots_by_group.pdf'))
    colours = ['c', 'burlywood', 'darkorange','crimson']
    
    for group in group_list:
        # fig, ax = plt.subplots(figsize=(20, 10))
        
        statistic, p_value = stats.kruskal(df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='0']['time'], 
                                           df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='1']['time'], 
                                           df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='2']['time'],
                                           df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='3']['time'])
    
        if p_value > 0.05:
            print ('no significance    p= ', p_value)
        else: 
            print('significant, continue with post hoc %%%%% add meDunn testing, p = ', p_value)
            data = [df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='0']['time'], 
                                           df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='1']['time'], 
                                           df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='2']['time'],
                                           df_summary_sns[df_summary_sns['group']==group][df_summary_sns['quad']=='3']['time']]
            
            p_quad_matrix = sp.posthoc_dunn(data, p_adjust = 'bonferroni')
            p_quad_matrix.insert(0, 'quadrents', np.arange(4))  # insert column at beginning 0,1,2,3
            #change column names from 1,2,3,4 to 0,1,2,3,
            dict = {1: '0',
                    2: '1',
                    3: '2', 
                    4:'3'}
             
            p_quad_matrix.rename(columns=dict,
                      inplace=True)
            
            
            # save  table to next page of pdf
            fig_p = plt.figure(figsize=(9,2))
            ax = plt.subplot(111)
            ax.axis('off')
            ax.table(cellText=p_quad_matrix.values, colLabels=p_quad_matrix.columns, bbox=[0,0,1,1])
            ax.set_title('Dunn Test for Group '+group)
            
            pdf.savefig(fig_p)
            # print(p_quad_matrix) ## change to save plotting starts and p values on graph ### SHIVA
        
        fig, ax = plt.subplots(figsize=(20, 10))   #new spot     
        sns.boxplot(x="quad", y="time", data=df_summary_sns [df_summary_sns['group'] == group], whis=np.inf, order = ['0','1','2','3'], palette=colours, saturation=0.6 ).set_title('Treatment_Group_'+group+'          p = '+ str(round(p_value, 4)))
        sns.swarmplot(x="quad", y="time", data= df_summary_sns [df_summary_sns['group'] == group], order = ['0','1','2','3'], palette=colours, marker = 'o', size = 6)
        plt.ylim(0,1)
        
        ### something wrong with the way this is done but i think it is just comparing between any two groups not like I have done above 
        # add_stat_annotation(ax, data=df_summary_sns,x="quad", y="time", order = ['0','1','2','3'],
        #             box_pairs=[("0", "1"), ("0", "2"), ("0", "3"), ("1", "2"), ("1", "3"), ("2", "3")],
        #             test='Kruskal', text_format='star', loc='outside', verbose=2)
        
        pdf.savefig(fig)
        
    
    pdf.close()
    plt.close('all') 
    
    return 


def remove_groups_from_df (df : pd.DataFrame,  groups_to_drop : list = []) -> pd.DataFrame:
    # to use groups to drop must be a list of strings refering to group i.e.  for dose responce  = ['5'.'6']
    df_croped : pd.DataFrame = df.copy()   
    for i in groups_to_drop:
        
        df_croped.drop(df_croped.index[df_croped['group'] == i], inplace=True)
    return df_croped

def save_dist_traveled_box_plot (df_summary : pd.DataFrame) -> None:
    
    pdf = PdfPages(os.path.join("./output",'distance_traveled_by_group.pdf'))
    
    fig, ax = plt.subplots(figsize=(20, 10))
    
    #alll 6 groups
    sns.boxplot(x="group", y="distance_traveled", data=df_summary, order = ['1', '2','3','4','5','6'], whis=np.inf, palette='pastel', saturation=0.6 ).set_title('Distance_traveled')
    sns.swarmplot(x="group", y="distance_traveled", data= df_summary, order = ['1', '2','3','4','5','6'], palette='pastel', marker = 'o', size = 6, edgecolor='k')
    pdf.savefig(fig)
    
    #for dose responce groups 1,2,3,4 
    df_dose_responce_dist_traveled = remove_groups_from_df(df_summary,  groups_to_drop = ['5','6'])

    ####SHIVA : there must be a more generic way to do this
    
    # stats non-parametric Kruskal-Wallis Test         
    data_dose_responce = [df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='1']['distance_traveled'], 
                  df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='2']['distance_traveled'], 
                  df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='3']['distance_traveled'],
                  df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='4']['distance_traveled']]
    
    statistic, p_value = stats.kruskal(df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='1']['distance_traveled'], 
                  df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='2']['distance_traveled'], 
                  df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='3']['distance_traveled'],
                  df_dose_responce_dist_traveled[df_dose_responce_dist_traveled['group']=='4']['distance_traveled'])
    
    if p_value > 0.05:
        print ('no significance    p= ', p_value)
    else: 
        print('significant, continue with post hoc Dunn testing, p = ', p_value)

        p_dose_responce_matrix = sp.posthoc_dunn(data_dose_responce, p_adjust = 'bonferroni')
        
        
        #save table of Dunn test
        fig_p = plt.figure(figsize=(9,2))
        ax = plt.subplot(111)
        ax.axis('off')
        ax.table(cellText=p_dose_responce_matrix.values, colLabels=p_dose_responce_matrix.columns, bbox=[0,0,1,1])
        ax.set_title('Dunn Test for Dose Responce')
        pdf.savefig(fig_p)
        
        ### SHIVA I would then like to add a star where it is significant on the boxplot below also reporting the p_value
        
    #boxplot 
    fig, ax = plt.subplots(figsize=(20, 10))
    sns.boxplot(x="group", y="distance_traveled", data=df_dose_responce_dist_traveled,   order = ['1','2','3','4'], whis=np.inf, palette='pastel', saturation=0.6 ).set_title('Distance_traveled_dose_responce p='+str(round(p_value,4)))
    sns.swarmplot(x="group", y="distance_traveled", data= df_dose_responce_dist_traveled, order = ['1','2','3','4'], palette='pastel', marker = 'o', size = 6, edgecolor='k')
    pdf.savefig(fig)
    
    
    #for agonist antagonist groups 1,3,5,6
    df_agg_antag_dist_traveled = remove_groups_from_df (df_summary,  groups_to_drop = ['2','4'])
    
    #non-parametric tetsing Kruskal-Wallis
    data_agg_antag = [df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='1']['distance_traveled'], 
                  df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='3']['distance_traveled'], 
                  df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='5']['distance_traveled'],
                  df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='6']['distance_traveled']]
    
    statistic, p_value = stats.kruskal(df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='1']['distance_traveled'], 
                  df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='3']['distance_traveled'], 
                  df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='5']['distance_traveled'],
                  df_agg_antag_dist_traveled[df_agg_antag_dist_traveled['group']=='6']['distance_traveled'])
    
    if p_value > 0.05:
        print ('no significance    p= ', p_value)
    else: 
        print('significant, continue with post hoc testing, p = ', p_value)
        #data is a list of lists
        
        p_agg_antag_matrix = sp.posthoc_dunn(data_agg_antag, p_adjust = 'bonferroni')  
        fig_p = plt.figure(figsize=(9,2))
        ax = plt.subplot(111)
        ax.axis('off')
        ax.table(cellText=p_agg_antag_matrix.values, colLabels=p_agg_antag_matrix.columns, bbox=[0,0,1,1])
        ax.set_title('Dunn Test for Agonist Antagonist')
        pdf.savefig(fig_p)
        
        
        ### SHIVA I would then like to add a star where it is significant on the boxplot below also reporting the p_value
        
    #boxplot     
    fig, ax = plt.subplots(figsize=(20, 10))
    sns.boxplot(x="group", y="distance_traveled", data=df_agg_antag_dist_traveled, order = ["1", "3", "5","6"],  whis=np.inf, palette='pastel', saturation=0.6 ).set_title('Distance_traveled_agg_antag p='+str(round(p_value,4)))
    sns.swarmplot(x="group", y="distance_traveled", data= df_agg_antag_dist_traveled, order = ["1", "3", "5","6"], palette='pastel', marker = 'o', size = 6, edgecolor='k')
    pdf.savefig(fig)
    
    
    pdf.close()
    plt.close('all') 
    return

#### Save table of mean percentage time and SEM in each quadrent for each treatment group

def save_table_of_quadrent_time (df_summary : pd.DataFrame):
    
    #get average of groupwise proportionality for each quadrent +- SEM
    mean_quadrent_percentage = df_summary.groupby('group', as_index=False)[['quad_0', 'quad_1','quad_2', 'quad_3']].mean()
    #calculate SEM and append to each group
    SEM_treatment = df_summary.groupby('group', as_index=False)[['quad_0', 'quad_1','quad_2', 'quad_3']].sem()
    SEM_treatment = SEM_treatment.drop(['group'], axis = 1)
    SEM_treatment = SEM_treatment.rename(columns={'quad_0': 'quad_0_SEM', 'quad_1': 'quad_1_SEM', 'quad_2': 'quad_2_SEM', 'quad_3': 'quad_3_SEM'})
    ####need to append sem data and check its really sem
    average_SEM_percentage_quadrent_data = pd.concat([mean_quadrent_percentage, SEM_treatment], axis=1, join='inner')
    average_SEM_percentage_quadrent_data.to_csv('./output/quadrent_percentage_+SEM_by_group.csv')
    
    
    return average_SEM_percentage_quadrent_data

def save_table_of_distance_traveled (df_summary : pd.DataFrame, plot_hist : bool = True):
    
    #get average of groupwise proportionality for each quadrent +- SEM
    mean_distance_traveled = df_summary.groupby('group', as_index=False)[['distance_traveled']].mean()
    #calculate SEM and append to each group
    n_each_treatment = df_summary.groupby('group', as_index=False)[['group']].size()
    n_each_treatment = n_each_treatment.drop(['group'], axis = 1)
    n_each_treatment = n_each_treatment.rename(columns={'size': 'n'})
    
    SEM_treatment = df_summary.groupby('group', as_index=False)[['distance_traveled']].sem()
    SEM_treatment = SEM_treatment.drop(['group'], axis = 1)
    SEM_treatment = SEM_treatment.rename(columns={'distance_traveled': 'distance_traveled_SEM'})
    
    average_SEM_distance_traveled_data = pd.concat([mean_distance_traveled, SEM_treatment, n_each_treatment], axis=1, join='inner')
    average_SEM_distance_traveled_data.to_csv('./output/distance_traveled_+SEM_by_group.csv')

    if plot_hist:
            
        plot_hist_distance_traveled(average_SEM_distance_traveled_data)

    return average_SEM_distance_traveled_data

#average_SEM_distance_traveled_data = save_table_of_distance_traveled (df_summary)


def plot_hist_distance_traveled(average_SEM_distance_traveled_data):

    pdf = PdfPages(os.path.join("./output",'distance_traveled_by_group_bar.pdf'))
    
    fig, ax = plt.subplots(figsize=(20, 10))

    '''plot histogram with SEM'''
    groups = list(average_SEM_distance_traveled_data['group'])
    CTEs = list(average_SEM_distance_traveled_data['distance_traveled'])
    error = list(average_SEM_distance_traveled_data['distance_traveled_SEM'])
    
    
    ax.bar(groups, CTEs, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
    ax.set_ylabel('Distance Traveled in cm')
    ax.set_xlabel('Treatment Group')
    ax.set_xticks(groups)
    ax.set_xticklabels(groups)
    ax.set_title('Mean Distance Traveled +/- SEM')
    ax.yaxis.grid(True)
    pdf.savefig(fig)
    
    pdf.close()
    plt.close('all') 
    return

### distance traveled 

def calculate_distance_traveled (df_filtered, cage_cords):
    
    y = df_filtered['CM','y'] 
    x = df_filtered['CM','x'] 
    y_2 = np.power( np.diff (y) , 2)
    x_2 = np.power( np.diff (x) , 2)
    distance_travelled = np.sqrt( x_2 + y_2 )
    sum_distance_travelled = np.sum (distance_travelled)

    sum_distance_travelled_cm = sum_distance_travelled # /pix_per_cm
    
    return sum_distance_travelled_cm



#### Discard short videos

def check_video_length_long_enough(df, threshold = 8, fps = 30):
    
    ''' Return True if length of video > threshold, False othervise'''
    
    length_in_min = round( len(df.index) / fps / 60, 1)
    video_not_long_enough = length_in_min < threshold
    
    # print('duration =', length_in_min , ' min', 
    #       ' Discarded = ',  video_not_long_enough)
    
    return video_not_long_enough

def create_empty_df_summaries_sns(column_names_sns, n_files):
    
    ''' Create empty df to be filled over the loop of videos'''
    

    return pd.DataFrame(np.zeros(( int( n_files * 4) , 4)), columns = column_names_sns)

def create_empty_df_summaries(column_names, n_files):
    
    ''' Create empty df to be filled over the loop of videos'''

    return pd.DataFrame(np.zeros((n_files, len(column_names))), columns = column_names)

### normalize to cm

def cal_side_length_in_pix( cage_coord ):
    ''' calculate the length of the side as the average of four sides and then return scale of pix to cm'''
    
    corner_pairs = [ 
        
                ['LRCorner', 'URCorner'],
                ['LLCorner', 'LRCorner'],
                ['LLCorner', 'ULCorner'],
                ['ULCorner', 'URCorner']
                
              ]
    ####### SHIVA previously this said LLCorner was paired with URCorner second pair
    
    side_lengths_pix = np.zeros(4)
    
    for i, corners in enumerate(corner_pairs):
        
        if np.sum(np.isnan(np.array ([
                                    cage_coord[ corners[0],'x'],
                                    cage_coord[ corners[1],'x'],
                                    cage_coord[ corners[0],'y'],
                                    cage_coord[ corners[1],'y']]
                                    )
                            )
                ) > 0:
            
            continue ### that side doesn't have high likelihood detections, skip
        
        side_lengths_pix[i] = (  ( cage_coord[ corners[0],'x'] - cage_coord[ corners[1], 'x'] ) ** 2 +                 
                                 ( cage_coord[ corners[0],'y'] - cage_coord[ corners[1], 'y'] ) ** 2
                              ) ** 0.5
    
    acceptable_length_inds = np.nonzero(side_lengths_pix)[0]
    n_acceptable_sides_detected = len( acceptable_length_inds )
    
    if  n_acceptable_sides_detected < 1:
        
        raise "no high likelihood for any two corners to get side length"
        
    # print('side lengths in pix = ', side_lengths_pix)    
    
    return np.average( side_lengths_pix[acceptable_length_inds] ) 

def cal_scale_pix_to_cm(side_length_cm, side_length_pix):
    
    scale =  side_length_cm / side_length_pix 
    # print( 'scale = ', round(scale, 5))
    
    return scale 

def scale_df_vals_to_cm(df : pd.DataFrame, scale_pix_to_cm : int, label_list : list , coord_list : list = ['x', 'y']) -> pd.DataFrame:
    
    ''' multiply each pixel value of each coordinate and label pair by the scale'''
    coord_label_pair : list = [(label,coord) for label in label_list for coord in coord_list]
    df[coord_label_pair] = df[coord_label_pair].apply(lambda x: x * scale_pix_to_cm)
    return df

def normalize_cage_coords(cage_coord, side_length_cm):
    
    side_length_pix = cal_side_length_in_pix( cage_coord )
    scale_pix_to_cm = cal_scale_pix_to_cm(side_length_cm, side_length_pix) #gives a scaler for each video of cm to pixels
    cage_coord_normlized = {key : v * scale_pix_to_cm for key, v in cage_coord.items()} # multiplies all valuesof the cage cords by scaler

    return cage_coord_normlized, scale_pix_to_cm

#### 
def drop_extra_rows_df(df, col_name_to_check_if_filled = 'mouse_no'):
    ''' remove extra rows in df'''
    
    ind_extra = np.where(df[col_name_to_check_if_filled] == 0) [0]
    return df.drop(ind_extra)


def cal_time_spent_percentage(df_filtered):
    
    time_spent_points = time_spent_in_each_quad (df_filtered)
    
    time_spent_percentage = [i / len(df_filtered) for i in time_spent_points]       #make time spent proportional not no. of points
            
    return time_spent_points, time_spent_percentage




def pre_process_df(df : pd.DataFrame, n_pts_to_keep, cage_coord_normalized : dict, scale_pix_to_cm, p_cutoff = 0.7, 
                   plot_missing_dist = False):
    

    ''' trim data according to time,  
                                filter with P, 
                                            scale pixels to cm, 
                                                    calculate  CM x,y and append, 
                                                                    calculate which quad and append'''

    def get_high_p_ind(df : pd.DataFrame) -> tuple:

        return  (( df['Nose']['likelihood'] > p_cutoff) & 
                ( df['Tail']['likelihood'] > p_cutoff) &
                ( df['LEar']['likelihood'] > p_cutoff) &
                ( df['REar']['likelihood'] > p_cutoff) )

    df_trimmed : pd.DataFrame = trim_df(df, n_pts_to_keep)


    high_p_ind : tuple = get_high_p_ind(df_trimmed)
    
    df_filtered : pd.DataFrame = filter_df(df_trimmed,high_p_ind)


    longest_gap : float = frames_missing_consec( np.where( high_p_ind ))

    n_pts_lost : int = len(df_trimmed) - len(df_filtered)

    if plot_missing_dist:
        total_length_of_video : int = len(df.index)
        bin_length: int = 30
        plot_missing_hist(high_p_ind, int( total_length_of_video / bin_length))


    
    label_list : list = list(set([key[0] for key in df_filtered.columns if key[0] != 'bodyparts'])) # the last header is 'which' we don't want that -> your missing columns is bodyparts 

    df_filtered : pd.DataFrame = scale_df_vals_to_cm(df_filtered, scale_pix_to_cm, label_list, coord_list = ['x', 'y'])
    
    plot_missing_time_pts : bool = False
    if plot_missing_time_pts:
        plot_time_points(df_filtered)

    df_filtered = calculate_and_add_mouse_CM(df_filtered)

    plot_trajectory: bool = False
    if plot_trajectory:
        plot_trajectory(df_filtered)

    zero: dict = get_zero_coord(cage_coord_normalized)

    cage_coord_normalized = {k : v - zero[k[1]] for k,v in cage_coord_normalized.items()}

    df_filtered = change_coord_rel_to_cage(df_filtered, zero)
    
    df_filtered = constrain_to_cage_coord(df_filtered, cage_coord_normalized) ############# BEWARE yo won't see the miss detections anymorte with this

    df_filtered = find_which_quad(df_filtered, cage_coord_normalized) ### some reoccuring warning try df_filtered[][] = [quad ] * len(df_filtered[][]) to fix SHIVA
    
    

    return df_filtered,  n_pts_lost, longest_gap, cage_coord_normalized


def get_zero_coord(cage_coord_normalized):
    return {'x': min(cage_coord_normalized['ULCorner', 'x'], cage_coord_normalized['LLCorner', 'x']),
            'y': min(cage_coord_normalized['ULCorner', 'y'], cage_coord_normalized['URCorner', 'y'])}

def change_coord_rel_to_cage(df : pd.DataFrame, zero : dict) -> pd.DataFrame:
    
    #subtract mean ULCorner x and y from df_filtered to noemalise for heat map
    label_list : list = list(set([key[0] for key in df.columns if key[0] != 'bodyparts'])) # the last headers are 'bodyparts' and 'which' zhich do not have x and y values to be normalised
    coord_list : list = ['x', 'y']
    coord_label_pairs : list = [[(label,coord) for label in label_list] for coord in coord_list]

    for coord,coord_label_pair in zip(coord_list,coord_label_pairs):
        df[coord_label_pair] = df[coord_label_pair].apply(lambda x : x - zero[coord])
            
    return df
    
def constrain_to_cage_coord(df : pd.DataFrame, cage_coord_normalized : dict) -> pd.DataFrame:
    
    right_bound : float = max(cage_coord_normalized['URCorner','x'], 
                        cage_coord_normalized['LRCorner','x'])
    left_bound : float = min(cage_coord_normalized['ULCorner','x'], 
                        cage_coord_normalized['LLCorner','x'])
    bottom_bound : float = max(cage_coord_normalized['LRCorner','y'], 
                        cage_coord_normalized['LLCorner','y'])
    upper_bound : float = min(cage_coord_normalized['URCorner','y'], 
                         cage_coord_normalized['ULCorner','y'])
    

    
    ind_x_over : int = df['CM','x'] > right_bound - left_bound
    
    ind_y_over : int = df['CM','y'] > bottom_bound - upper_bound
    
    ind_x_less : int = df['CM','x'] < 0
    
    ind_y_less : int = df['CM','y'] < 0
    
    print(sum(ind_x_over), sum(ind_y_over),
          sum(ind_x_less), sum(ind_y_less))
    
    df['CM','x'][ind_x_over] = right_bound
    df['CM','y'][ind_y_over] = bottom_bound
    df['CM','x'][ind_x_less] = left_bound
    df['CM','y'][ind_y_less] = upper_bound

    return df

def fill_df_summaries(df_summary, df_summary_sns, mouse_no, group, time_spent_percentage, 
                      n_pts_lost, longest_gap, distance_traveled, i):
    
    df_summary = fill_df_of_each_group_analysis(df_summary,mouse_no, group, time_spent_percentage, 
                                                        n_pts_lost, longest_gap, distance_traveled, i = i)
            
            
    df_summary_sns = fill_df_of_each_group_quad_concat(df_summary_sns, mouse_no, group, time_spent_percentage,
                                                                i = int( i * 4))
    return df_summary, df_summary_sns

def grubbs_test_DLC(df_summary, treatment_dict, list_of_groups, p_value = 0.05, remove_all = False):
    ''' 
    For every treatment Grubbs test is run, excluded data reported and removed from df_inc_ratios    
    per treatment group and a single outlier per group 
    '''
    # count = 0
    for treatment in list_of_groups:
        
        print('working on treatment ', treatment)
        
        df2=df_summary[df_summary['group']== treatment ]['distance_traveled']
        x = np.array(df2)
        
        
        
                    
        if remove_all == True :
            x_ = grubbs.test(x, alpha= p_value ) #removing all outliers not just one #https://pypi.org/project/outlier_utils/
        
        else:
            G_calculated, outlier = grubbs_calculted_value (x)
            G_critical = grubbs_critical_value(len(x), p_value)
            
            if G_critical < G_calculated :
                ''' we can reject the Ho: there are no outliers in the set '''
                x_ = np.delete(x, np.where(x == outlier))
            else:
                x_ = x
        
        if len(x) != len(x_):
            print ('Grubbs Test exclusion for ', treatment_dict[int(treatment)])
            # excluded = float(np.setdiff1d(x,x_)
            print('mean of new set = ', x_.mean(), 'excluded point = ', np.setdiff1d(x,x_))
            
            mouse_outliers = []
            for outlier in np.setdiff1d(x,x_):
                ind_of_outllier = df_summary.loc[df_summary['distance_traveled'] == float(outlier)].index[0]
                mouse_outliers.append( df_summary.loc[ind_of_outllier]['mouse_no'] )
                
            # print(ind_of_outllier)
            print('mouse number ', mouse_outliers , ' excluded')
            
            df_summary['distance_traveled'].replace( np.setdiff1d(x,x_) , np.NaN, inplace=True)
                
    return df_summary


def grubbs_calculted_value (y):# https://www.youtube.com/watch?v=Hn_lMUaMcak&ab_channel=BhaveshBhatt
    avg_y = np.mean(y)
    abs_val_minus_avg = abs(y - avg_y)
    max_of_deviations = max(abs_val_minus_avg)
    outlier = max_of_deviations+avg_y
    # outlier_ind = np.where(abs_val_minus_avg == max_of_deviations)
    s = np.std(y)
    G_calculated = max_of_deviations/ s
    
    return G_calculated, outlier

def grubbs_critical_value(size, alpha):
    """Calculate the critical value with the formula given for example in
    https://en.wikipedia.org/wiki/Grubbs%27_test_for_outliers#Definition
    Args:
        ts (list or np.array): The timeseries to compute the critical value.
        alpha (float): The significance level.
    Returns:
        float: The critical value for this test.
    """
    t_dist = stats.t.ppf(1 - alpha / (2 * size), size - 2)
    numerator = (size - 1) * np.sqrt(np.square(t_dist))
    denominator = np.sqrt(size) * np.sqrt(size - 2 + np.square(t_dist))
    G_critical_value = numerator / denominator
    # print("Grubbs Critical Value: {}".format(critical_value))
    return G_critical_value


def draw_cage(cage_coord_normalized, ax):
    
    sides = [['ULCorner', 'URCorner'],
             ['ULCorner', 'LLCorner'],
             ['URCorner', 'LRCorner'],
             ['LRCorner', 'LLCorner']]
    # side_coords = [[[corn[0], coord] for coord in ['x','y']] for corn in sides]
    for corner in sides:
        
        ax.plot([ cage_coord_normalized[corner[0],'x'], cage_coord_normalized[corner[1],'x']],
                [ cage_coord_normalized[corner[0],'y'], cage_coord_normalized[corner[1],'y']], '-o', c = 'k')
    ##
    
def produce_heat_map_all_groups(X_Y_df_file_list : list, teratment_dict,  folder_path, name = '_') -> None:
    
    pdf = PdfPages( os.path.join("./output",name +'Heat_Maps_per_treatment.pdf'))
    
    for filename in X_Y_df_file_list:
        
        filepath = os.path.join(folder_path, 'CM_analysis', filename)
        df = pd.read_csv(filepath, header = [0])
        
        fig, ax = plt.subplots()
        
        ###### need to scale for the number of points lost # use max occurnce with 1 square for max of colour map np.hist2D
        print ('mouse = ', filename, ' n = ', len(df))
   
        # fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(20, 8))

        
        nbins = 50
        # hist1 = axes[0].hist2d(df['x'], df['y'], bins=nbins, density = True, cmap='magma')
        
        hist1 = ax.hist2d(df['x'], df['y'], bins=nbins, density = True, cmap= plt.cm.BuPu)
        # fig.colorbar(hist1[3], ax = ax)
        
        
        H, xedges, yedges = np.histogram2d(df['x'], df['y'], bins=(nbins, nbins),density=True )
        # H_normalized = H/float(az1.shape[0]) # the integral over the histogrm is 1
        H_normalized = H/H.max((0,1)) # the max value of the histogrm is 1
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
        # im = axes[1].imshow(H_normalized, extent=extent, cmap='magma', interpolation='none',origin ='lower')
        # fig.colorbar(im, ax=axes[1])
        
        im = ax.imshow(H_normalized, extent=extent, cmap=plt.cm.BuPu, interpolation='none',origin ='lower')
        fig.colorbar(im, ax=ax)
        
        plt.show()

        treatment = int(filename[0:1])
        ax.set_title(teratment_dict[treatment])
            
        ax.invert_yaxis() # making same as video
        
        pdf.savefig(fig)
        
        plt.close(fig)

    pdf.close()
    
    

    return 
        


def produce_trajectory_all_groups(X_Y_df_file_list : list, cage_coord_normalized, folder_path: str, name : str = '_') -> None:
    
    pdf = PdfPages(os.path.join("./output", name +'Trajectory_per_treatment.pdf'))

    for filename in X_Y_df_file_list:
        
        filepath = os.path.join(folder_path, 'CM_analysis', filename)
        df = pd.read_csv(filepath, header = [0])
        mice_no_set = df['mouse_no'].unique()
        
        fig, ax = plt.subplots()
        
        color = iter(cm.rainbow(np.linspace(0, 1, len(mice_no_set))))
        
        for mouse_no in mice_no_set:
        # for mouse_no in mice_no_set[:1]:    # plot for a single mouse 
            c = next(color)
            mouse_ind = df['mouse_no'] == mouse_no
            ax.plot(df['x'][mouse_ind], df['y'][mouse_ind], lw = 0.5, c = c, label = mouse_no)
            draw_cage(cage_coord_normalized, ax)

        ax.set_title('Tretment Group ' + filename[0:1])
        ax.legend(fontsize= 5 )
        ax.invert_yaxis() # making same rotation as video
        
        pdf.savefig(fig)
        
        plt.close(fig)

    pdf.close()
    
    return

# produce_trajectory_all_groups(X_Y_df_file_list)

