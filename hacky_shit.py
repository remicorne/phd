#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 16:55:10 2023

@author: jasminebutler
"""

# database language if needed -  SQL

# Jupiternotebook * just displaying (matplotlib) * not generating PDF's saves time

#%% how to change working directory in console 

# pwd #pwd is the working directory of the terminal/consile
# Out[7]: '/Users/jasminebutler/Desktop'

# cd phd/ #cd = change working directory
# /Users/jasminebutler/Desktop/phd


# ls / # get a list of everything bellow 
# Applications/ Volumes/      etc@          sbin/
# Library/      bin/          home@         tmp@
# System/       cores/        opt/          usr/
# Users/        dev/          private/      var@

# cd .. #to go up 
# /Users/jasminebutler/Desktop


#%% df PLAY 

df.query("'drug' == 'PSIL' and 'data_type' == 'FP'") #returns T or F

#single slicing
drd_df = drd_df.loc[drd_df['data_type'] == 'FP']

#saving a df 
feature_df_ex_tau.to_csv('middle_TCB_data.csv')

#find max in col 
drd_df['rheobased_threshold'].idxmax() #e.g.268
# call col value by index
drd_df['rheobased_threshold'][268]

#%%

'''
Output from the code as of August 2023:
    
    'ouliers_by_treatment_compounds.xlsx' : mice removed from each treatment group using gubbs test for compounds and ratrios
    'ouliers_by_treatment_ratios.xlsx'
    
    'df_compound_mean', 'df_compound_SD', 'df_compound_SEM', 'df_comp_ratios_mean', 'df_comp_ratios_SD', 'df_comp_ratios_SEM'': mean SEM and SD tables
    
    'shapiro_data_all_treatments.xlsx' : shapiro wilk normality test for each treatment group and comp_BR
    
    'TUKEY_dose_resp_compounds.xlsx' : tukey test of psignificant one way anovas
    'TUKEY_dose_resp_ratios.xlsx'
    
    'one_way_ANOVA_ifsig_Tukey_compounds.xlsx' : each tukey test following sig ANOVA and last sheet is all ANOVA
    'one_way_ANOVA_ifsig_Tukey_ratios.xlsx'
    
    
'''


#%% creating pdf's

# with PdfPages('max_firing_by_mouseline.pdf') as pdf:
#     or
# pdf = PdfPages(name+'_histograms_.pdf')


#%% # #MAXs soloution ALOTHER DF EXPANSION SOLOUTION



# #add emptyy columns to append extracted values to FP data 
# feature_df = feature_df.assign(max_firing=np.NaN, 
#                                rheobased_threshold = np.NaN, 
#                                voltag_threshold = np.NaN, 
#                                AP_height = np.NaN, 
#                                AP_width = np.NaN, 
#                                AP_slope = np.NaN, 
#                                FI_slope = np.NaN, 
#                                tau_rc = np.NaN, 
#                                sag = np.NaN)



# # combinations = feature_df.folder_file.unique()

# combinations = product(feature_df.data_type.unique(),
#                 feature_df.drug.unique())


# for data_type, drug in combinations:
    
#     #take slice of data 
#     sliced_df= feature_df.loc[(feature_df.data_type== 'FP') &
#                               (feature_df.drug == 'PRE')]
    
#     # sliced_df= feature_df.loc[(feature_df.folder_file== file ) ] # single row
    
#     print(sliced_df.data_type)

#     if sliced_df.data_type[0] == 'FP' :
#         print('hi')
    
    
#     or 'FP_AP':
        
#         path = 
#         func
#         folder_file = ' '  
#         out = []
#         ...
        
#     elif data_type == 'AP':
#         continue


#     # call on a specific set of rows and set to whatever 
#     feature_df = feature_df.loc[(feature_df.data_type==a) &
#                    (feature_df.drug == b), 'max_firing'] = max_firing

              