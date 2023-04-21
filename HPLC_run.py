#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 14 15:15:02 2022

@author: jasminebutler
"""

# %% LE PLAN

'''
    1. create rations for each brain region: DOP/DAP, HIA/5HT, 3MT/DAP
    2. calculate and save to table average, SD and SEM
    3. ShapiroWilk test for normailty (note if not noirmal) and statistical analysis
    4. plot histograms with SEM and ** for each treatment, compound and BR
    5. test normaility for correlation if normal : 
                                                    --> pearson product coirrelation plots 
                            non-parametric data :
                                                    --> spearman coirrelation plots      
                                                    
            produced for each BR comparing compunds and compound comparing BR
    6. PCA for  complete compound_BR/ ratio_BR sets between groups 


'''


# %% OPEN DATA

from sklearn.preprocessing import StandardScaler # mean = 0 vairance =1
from sklearn.decomposition import PCA
from HPLC_module import *  # must be in the same directory

path = os.getcwd() + '/TCB2_data_APRIL23.csv'  # TCB2 #using current working directory plus file name 


df_raw = pd.read_csv(path, header=0)

df = df_raw.copy() #cworking df 
df = df.replace(np.nan, 0) # to set all 0 to Nan

#%%  DEVELOPING POS2.0 new df structure

# ### restructure df = mouse_id, treatment, BR, compound, ng_mg
col_names = ['mouse_id' , 'treatment' , 'BR' , 'compound' , 'ng_mg']
raw_col_names = df_raw.columns.copy()
result_rows = []

for ind, row in df_raw.iterrows(): #loop for row get mouse and group
    
    mouse_id = row[0]
    treatment = row[1]

    for val, col_name in zip(row[2:], raw_col_names[2:]): #loop within row

        if val > 0:
            compound, BR = col_name.split('_')
            result_rows.append([mouse_id, treatment, BR, compound, val])

restructured_df = pd.DataFrame(result_rows, columns=col_names)


# %% CREATE WORKING DF
print(ind, row)
index_brain_region = list(df.head(0))

# first two elements are mouse and group info
ideal_headers = [i.split('_') for i in index_brain_region]
compound_list, brain_region_list = zip(*ideal_headers[2:])


# separate column headers into multi index
df.columns = pd.MultiIndex.from_tuples(ideal_headers)

list_of_brain_regions = np.unique(brain_region_list)  # TCB2 == 32
list_of_compounds = np.unique(compound_list)  # TCB2 == 26
list_of_groups = np.unique(df['group', 'no'])
list_of_mice = np.unique(df['mouse', 'no'])


# create dict of existing indicies to handel missing values in data set
exist_dict = get_summary_of_exsisting_data(
    df, list_of_brain_regions, list_of_groups, list_of_mice)



# %%% 1. create rations for each brain region: DOPAC/DA, HIA/5HT, 3MT/DA

'''
Define compound_ratio_dict as the pairs of compounds to be created for analysis 

ratios are then ceated and appended to df_inc_ratios and a dictionary is created with the same format as for quantative data 
'''

# ratios to be included in analysis

compound_ratio_dict = {'DOPAC': ['DA'],     # TCB2
                       '5HIAA': ['5HT'],
                       '3MT': ['DA'],
                       'HVA': ['DA', '3MT', 'DOPAC'],
                       'GLN': ['GLU']}


# compound_ratio_dict = {'DO' : [ 'DA'],     #SNORD115116         #ALZ_mice
#                   '5HI' : ['5HT'],
#                   'HVA' : ['DA', 'DO']}

df_inc_ratios_raw, list_of_ratios, ratio_match_ind_dict = create_ratios_for_each_BR(df,
                                                                                    list_of_brain_regions,
                                                                                    compound_ratio_dict,
                                                                                    exist_dict,
                                                                                    len(df.index),
                                                                                    list_of_groups)



'''
##### TREATMENT LABLES   and   COLOUR PALLET #####
 #to be used for all plots - treatment dict name e.g. WT and colour
'''
# #ALZ_mice
# treatment_dict = { 1: 'young_WT', 2 : 'young_AD', 3: 'old_WT', 4:'old_AD'}
# palette_labeled = { 'young_WT': 'white',  'young_AD': 'orange',                                       #ALZ_mice
#                     'old_WT': 'grey', 'old_AD': 'red'}

treatment_dict = {1: 'vehicles', 2: '10mg/kgTCB', 3: '3mg/kgTCB',
                  4: '0.3mg/kgTCB', 5: 'TCB+MDL', 6: '0.2mg/kgMDL'}  # TCB2

# treatment_dict = { 1: 'WT', 2 : 'SNORD_115116_KO'}       #SNORD115116

palette_labeled = {'vehicles': "white",  # TCB2
                   '10mg/kgTCB': "firebrick",
                   '3mg/kgTCB': "red",
                   '0.3mg/kgTCB': "salmon",
                   'TCB+MDL': "grey",
                   '0.2mg/kgMDL': "black"}

# palette_labeled = { 'WT': 'white',        #SNORD115116
#                     'SNORD_115116_KO': 'red'}


# %% TEST FOR OUTLIERS, Grubbs test (same as graphpad outliers)

df_inc_ratios = df_inc_ratios_raw.copy()

# GRUBBS TEST

'''
input: 
        raw data : df_inc_ratios_raw
        dictionary of data : exist_dict or ratio_match_ind_dict
        treatment groups : list_of_groups
        p_value / alpha : 0.05
        remove_all: if TRUE Grubbs test will be applied as in outlier_utils (remove 1-2 outliers max and min?) otherwise if FALSE will only remove most extreme
'''

# reporting and removing data from all  BR and compound combinations
df_inc_ratios, outlier_dict = grubbs_test(
    df_inc_ratios, exist_dict, list_of_groups, treatment_dict, p_value=0.05, remove_all=True)

# reporting and removing data from all BR ratio combinations
df_inc_ratios, outlier_dict_ratio = grubbs_test(
    df_inc_ratios, ratio_match_ind_dict, list_of_groups, treatment_dict, p_value=0.05, remove_all=False)


save_outlier_info(outlier_dict, name='compounds')
save_outlier_info(outlier_dict_ratio, name='ratios')


# %% 2. calculate and save to table average, SD and SEM
'''
Just seems to make global vairables not saved, do we use them later or why generate?  ## FIX ME

So, this is not currently in use but was use to graph with SEM rathere than with CI (as our data is largly not parametric) POS2.0 will resolve this issue for now only use it to get number if needed 
'''

df_compound_mean, df_compound_SD, df_compound_SEM = calculate_mean_SD_SEM(
    df_inc_ratios, exist_dict)
df_comp_ratios_mean, df_comp_ratios_SD, df_comp_ratios_SEM = calculate_mean_SD_SEM(
    df_inc_ratios, ratio_match_ind_dict)

# WARNING
# See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
# summary['SEM'][comp, BR] [int(group) - 1] =  df_inc_ratios[comp, BR][group_ind].sem()


# current bug re dictionary recognising missing data and excluding from plotting
# agonist_antagonist: compounds and ratios

# plot_hist_comparing_treatment_SEM (treatment_dict, exist_dict,  palette_labeled, df_compound_mean, df_compound_SEM,      # missing the significance stars
#                                     treatments_to_compare = [1,3,5,6], name = 'tester',
#                                     order = ['vehicles', '3mg/kgTCB','0.2mg/kgMDL','TCB+MDL' ],
#                                     test_path = '/Users/jasminebutler/Desktop/TUKEY_ag_ant_compounds.xlsx',
#                                     ratio = False)

#saving as csv
list_df_to_save = [df_compound_mean, df_compound_SD, df_compound_SEM, df_comp_ratios_mean, df_comp_ratios_SD, df_comp_ratios_SEM ]
list_names = ['df_compound_mean', 'df_compound_SD', 'df_compound_SEM', 'df_comp_ratios_mean', 'df_comp_ratios_SD', 'df_comp_ratios_SEM']

for  df, name  in zip(list_df_to_save, list_names):
    
    df.to_csv(name + 'TCB2.xls')

# %% 3. ShapiroWilk test for normailty (note if not noirmal) and statistical analysis
'''' 
##############   NORMAILITY TESTING FORALL     ##############
#calculate normaility for each treatment group, save excel doc with F and p value if p<0.05 the data set will be printed in console
# p<0.05 we reject null hypothesis (that data is normaly distributed)

Outputs : 
    shapiro_data_all_treatments.xlsx      
'''

save_shapiro_data_for_all_treatments(
    exist_dict, list_of_brain_regions, list_of_compounds, df_inc_ratios, list_of_groups)


# %% 4. quantative statistical analysis between groups

'''
##############       MULTIPLE GROUPS       ##############
#DOSE/RESPONCE (e.g. groups 1,2,3,4) == one way anova followed by post hoc tukey test (a function for ratios and compounds sepraly)
#AGONIST/ANTAGONIST (e.g groups 1,3,5,6) == two way anova followed by one way anova and post hoc tukey (combined function for MA and ratios)

Outputs:
    dose_responce_anova_ratios/MA/AA.csv           single table with all F and p values from one way anova
    dose_responce_tukey_table_ratios/MA/AA.xlsx    significant one way anova tukey tests
    
    two_way_ANOVA_agonist_antagonist_significant.xlsx   
    TUKEY_agonist_antagonist_significant.xlsx
'''
#saving ALL one way ANOVAs for PDD
onewayANOVA_Tukeyposthoc(df_inc_ratios, exist_dict, groups_to_drop=[5, 6], name='compounds') 
onewayANOVA_Tukeyposthoc(df_inc_ratios, ratio_match_ind_dict, groups_to_drop=[5, 6], name='ratios')

# #old way with spit out seperate MAYBE NEED FOR PLOTS
one_way_ANOVA_post_hoc_between_groups(
    df_inc_ratios, exist_dict, groups_to_drop=[5, 6], name='compounds')

one_way_ANOVA_post_hoc_between_groups(
    df_inc_ratios, ratio_match_ind_dict, name='ratios', groups_to_drop=[5, 6])

##############
# #ALZ_mice
# two_way_ANOVA_ALZ (exist_dict, df_inc_ratios, treatments_to_compare = [1,2,3,4], name = 'compounds')  #FIX ME?
# two_way_ANOVA_ALZ (ratio_match_ind_dict, df_inc_ratios, treatments_to_compare = [1,2,3,4], name = 'ratios')

#saving ALLLLLLLL two way and one way anova and tukey only for sig one way  RUN LAST 
twoway_ANOVA_oneway_ANOVA_ifsig_Tukey(exist_dict, df_inc_ratios, treatments_to_compare=[1, 3, 5, 6], name='compounds')
twoway_ANOVA_oneway_ANOVA_ifsig_Tukey(exist_dict, df_inc_ratios, treatments_to_compare=[1, 3, 5, 6], name='compounds')

#OLD way maybe NEED FOR HISTOGRAM PLOTS 
two_way_ANOVA(exist_dict, df_inc_ratios, treatments_to_compare=[
              1, 3, 5, 6], name='compounds')

two_way_ANOVA(ratio_match_ind_dict, df_inc_ratios,
              treatments_to_compare=[1, 3, 5, 6], name='ratios')


'''
##############       TWO GROUPS       ##############
#unpaired student t-test e.g. comparing two genetic lines 

Outputs:
    compounds_t_test_values.csv          t test t and p values between two groups
    ratios_t_test_values.csv          t test t and p values between two group
'''

# df_t_test_results = student_t_test(exist_dict, df_inc_ratios, name = 'compounds')

# df_t_test_results_ratios = student_t_test(ratio_match_ind_dict, df_inc_ratios, name = 'ratios')


# %% 4. plot histograms with SEM and ** for each treatment, compound and BR


''''
##############        MULTIPLE TREATMENTS/GROUPS       ##############
# plot histograms of quantative data witrh statistical significance (calculated prior) indicatde with stars +/- 98% CI --> SEM
Outputs : 

    
    agonist-antagonist 
    
    dose-responce
    
'''

start_time = time.perf_counter()


# ######ALZ_mice######
# plot_hist_comparing_treatment_CI (treatment_dict, exist_dict, df_inc_ratios, palette_labeled,
#                                 treatments_to_compare = [1,2,3,4], name = 'ALZ_mice_compounds',
#                                 order = ['young_WT',  'young_AD', 'old_WT','old_AD'],
#                                 test_path = '/Users/jasminebutler/Desktop/TUKEY_ag_ant_compounds.xlsx', ratio = False, mouse_id = False )
# plot_hist_comparing_treatment_CI (treatment_dict, ratio_match_ind_dict, df_inc_ratios, palette_labeled,
#                                 treatments_to_compare = [1,2,3,4], name = 'ALZ_mice_ratios',
#                                 order = ['young_WT',  'young_AD', 'old_WT','old_AD'],
#                                 test_path = '/Users/jasminebutler/Desktop/TUKEY_ag_ant_ratios.xlsx', ratio = False, mouse_id = False )
# ALZ_mice###### ^

plot_hist_comparing_treatment_CI(treatment_dict, exist_dict, df_inc_ratios, palette_labeled,
                                 treatments_to_compare=[1, 3, 5, 6], name='Agonist_Antagonist_compounds',
                                 order=['vehicles', '3mg/kgTCB',
                                        '0.2mg/kgMDL', 'TCB+MDL'],
                                 test_path=os.getcwd() + '/TUKEY_ag_ant_compounds.xlsx', ratio=False, mouse_id=False)

plot_hist_comparing_treatment_CI(treatment_dict, ratio_match_ind_dict, df_inc_ratios, palette_labeled,
                                 treatments_to_compare=[1, 3, 5, 6], name='Agonist_Antagonist_ratios',
                                 order=['vehicles', '3mg/kgTCB',
                                        '0.2mg/kgMDL', 'TCB+MDL'],
                                 test_path=os.getcwd() + '/TUKEY_ag_ant_ratios.xlsx', ratio=True, mouse_id=False)

# dose_responce: compounds and ratios
plot_hist_comparing_treatment_CI(treatment_dict, exist_dict, df_inc_ratios, palette_labeled,
                                 treatments_to_compare=[1, 2, 3, 4], name='Dose_Responcxe_compounds',
                                 order=['vehicles', '0.3mg/kgTCB',
                                        '3mg/kgTCB', '10mg/kgTCB'],
                                 test_path=os.getcwd() + '/TUKEY_dose_resp_compounds.xlsx', ratio=False, mouse_id=False)

plot_hist_comparing_treatment_CI(treatment_dict, ratio_match_ind_dict, df_inc_ratios, palette_labeled,
                                 treatments_to_compare=[1, 2, 3, 4], name='Dose_Responce_ratios',
                                 order=['vehicles', '0.3mg/kgTCB',
                                        '3mg/kgTCB', '10mg/kgTCB'],
                                 test_path=os.getcwd() + '/TUKEY_dose_resp_ratios.xlsx', ratio=True, mouse_id=False)


end_time = time.perf_counter()

print(f"Execution Time :  {end_time - start_time:0.6f}")
'''
##############    TWO GROUPS     ##############
#histogram of tissue quantity showing significance using unpaired student t-test

Outputs:
    componds_histograms_.pdf            
    ratios_histograms_.pdf

'''

# plot_two_groups (treatment_dict, exist_dict, df_inc_ratios, palette_labeled, df_t_test_results,
#                      name = 'componds', treatments_to_compare = [1,2], p = 0.05)

# plot_two_groups (treatment_dict, ratio_match_ind_dict, df_inc_ratios, palette_labeled, df_t_test_results_ratios,
#                      name = 'ratios', treatments_to_compare = [1,2], p = 0.05)


# %% 5.  --> pearson product coirrelation plots

######### JJB testing/ checking

# >>> l = ['a', 'b', 'c']
# >>> all(x in l for x in ['a', 'b'])
# True
# >>> all(x in l for x in ['a', 'x'])
# False


# all(x in list_of_compounds for x in column_order )

# Correlations within a BR for each treatment group

# PN: HSER has been removed manuly and will not run in a general way! using if function?

'''
Output: 
         compound_correlograms_within_BR.pdf
         ratio_correlograms_within_BR.pdf
'''

# Pearson product correlograms with only compunds in column_order  ###### JAS23 something wrong with ordering
pearson_correlations_within_BR(list_of_groups, exist_dict, df_inc_ratios, treatment_dict, name='compound', p_value=0.05,
                               column_order=['A', 'NA', 'VMA', 'DA', 'DOPAC', 'LDOPA', 'HVA', '3MT', '5HT', '5HIAA', '5HTP'], method='pearson')  # ordered data



# pearson_correlations_within_BR (list_of_groups, exist_dict, df_inc_ratios,treatment_dict, name = 'compound', p_value = 0.05,            ### all data  unordered in correlogram
#                                 column_order = False)

# Pearson product correlograms with ordered ratios #something wrong
pearson_correlations_within_BR(list_of_groups, ratio_match_ind_dict, df_inc_ratios, treatment_dict, name='ratio', p_value=0.05,
                               column_order=False, method='pearson')
#TypeError: ufunc 'invert' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''

'''
Output: 
            compound_correlograms_between_BRs.pdf
            ratio_correlograms_between_BRs.pdf
'''
#########
# Correlations between BR's for each treatment group

# Pearson product correlograms with all compunds

# col_order_old = ['OF', 'PL', 'IC', 'SJ', 'SL1', 'SR1', 'SL6', 'SR6', 'AC', 'V', 'M', 'CC', 'NAc', 'VM', 'DM',
#              'VL', 'DL', 'Am', 'dH', 'vH', 'Y',   'SN', 'VTA', 'DR', 'MR', 'VPL', 'VPR', 'MD', 'SC', 'DG', 'CE']

col_order = ['OF', 'PL', 'CC', 'IC', 'M', 'SJ', 'SL1', 'SR1', 'SL6', 'SR6', 'AC', 'V' , 
       'Am', 'dH', 'vH', 'NAc', 'VM', 'DM', 'VL', 'DL', 'MD', 'VPL', 'VPR', 'DG', 
       'Y', 'SC', 'SN', 'VTA', 'DR', 'MR', 'CE' ]  #pdd august 2023

pearson_correlations_within_compound(list_of_groups, exist_dict, df_inc_ratios, list_of_compounds, treatment_dict,
                                     name='compound', p_value=0.05, column_order=col_order, method='pearson')


# for SNORD? I dont know old list
# col_order = ['OF', 'PL','Cg', 'IC', 'SJ', 'SL1', 'SR1', 'SL6', 'SR6', 'NAc', 'VM','DM', 'VL', 'DL', 'SP','Am', 'HD', 'HV','Th','DY',  'VY',  'SN', 'VTA','DR',  'MR']
# ? Cg > CC , SP?, DY > Y (dorsal and ventral hypothalamus?)


# pearson_correlations_within_compound (list_of_groups, exist_dict, df_inc_ratios, list_of_compounds, treatment_dict,
#                                       name = 'compound', p_value = 0.05,n_minimum=8, column_order = False)

# ALZ mice
# col order for ALZ: ['MO', 'LO', 'M2', 'PL', 'IL', 'aCg', 'pCg', 'S1', 'Ent', 'A', 'dH', 'vH', 'sH', 'Co', 'VM', 'DM', 'VL', 'DL', 'GP', 'vIT', 'MD', 'HYP', 'CB', 'SN', 'VTA', 'DR', 'MR']
# pearson_correlations_within_compound (list_of_groups, exist_dict, df_inc_ratios, list_of_compounds, treatment_dict,
#                                       name = 'compound', p_value = 0.05, column_order = ['MO', 'LO', 'M2', 'PL', 'IL', 'aCg', 'pCg', 'S1', 'ENT', 'A', 'DH', 'VH', 'SH', 'CO', 'VM', 'DM', 'VL', 'DL', 'GP', 'vlT', 'MD', 'HYP', 'CB', 'SN', 'VTA', 'DR', 'MR'], method = 'pearson')
# pearson_correlations_within_compound (list_of_groups, ratio_match_ind_dict, df_inc_ratios, list_of_ratios, treatment_dict,
#                                       name = 'ratios', p_value = 0.05, column_order = ['MO', 'LO', 'M2', 'PL', 'IL', 'aCg', 'pCg', 'S1', 'ENT', 'A', 'DH', 'VH', 'SH', 'CO', 'VM', 'DM', 'VL', 'DL', 'GP', 'vlT', 'MD', 'HYP', 'CB', 'SN', 'VTA', 'DR', 'MR'], method = 'pearson')


# ALZ ^mice

# Pearson product correlograms with all ratios

pearson_correlations_within_compound(list_of_groups, ratio_match_ind_dict, df_inc_ratios, list_of_ratios, treatment_dict,
                                     name='ratios', p_value=0.05, column_order=col_order, method='pearson')

# pearson_correlations_within_compound(list_of_groups, ratio_match_ind_dict, df_inc_ratios, list_of_ratios, treatment_dict,
#                                      name='ratios', p_value=0.05, column_order=False, method='pearson')


#########
# Correlations between TWO BR's for each treatment group SQUARE CORRELOGRAM

# pearson_correlations_between_two_BR (treatment_dict, list_of_groups, exist_dict,
#                                            df_inc_ratios, name = 'neurotransmitters_between_VL_VM', p_value = 0.05,
#                                            compounds_to_analise = ['NA','5HT','DA', 'GLU', 'GABA', 'GLY'],
#                                            BR_to_analise = ['VL', 'VM'])

# %% 5.  --> Principal Component Analysis

# colapse multi index headding
# df_inc_ratios_PCA = df_inc_ratios.copy()


# MA_list = ['3MT', '5HIAA', '5HT', '5HTP', 'DA', 'DOPAC', 'HVA', 'NA']
# AA_list = ['ALA', 'ARG', 'ASN', 'ASP', 'GABA', 'GLN', 'GLU',
#            'GLY', 'HIS', 'HSER', 'LSER', 'TAU', 'THR', 'TYR']
# # list_of_ratios

# df_inc_ratios_PCA.columns = ['_'.join(col).strip()
#                              for col in df_inc_ratios_PCA.columns.values]


# ratio_features = [
#     col for ratio in list_of_ratios for col in df_inc_ratios_PCA.columns if ratio in col]

# # MA ratios included
# MA_cols = [col for MA in MA_list for col in df_inc_ratios_PCA.columns if MA in col]
# MA_features = [
#     features for features in MA_cols if features not in ratio_features]

# # AA ratios
# AA_cols = [col for AA in AA_list for col in df_inc_ratios_PCA.columns if AA in col]
# AA_features = [
#     features for features in AA_cols if features not in ratio_features]


# # create df cntaining only complete data sets and report which mice exist for each group
# # will also include a a filter of MA , ratio or AA that would like to be kept for all BR i.e. 3MT is missing for a lot of groups and should not be used for PCA

# # def normalisation_of_features (MA_features, df_inc_ratios_PCA):
# # normalise/ standardise data

# features = MA_features

# # Separating out the features
# x = df_inc_ratios_PCA.loc[:, features].values
# # Separating out the target
# y = df_inc_ratios_PCA.loc[:, ['group_no']].values
# # Standardizing the features
# x = StandardScaler().fit_transform(x)


# # PCA

# pca = PCA(n_components=2)
# principalComponents = pca.fit_transform(x)
# principalDf = pd.DataFrame(data=principalComponents, columns=[
#                            'principal component 1', 'principal component 2'])
# # mean centering for each comp and BR = subtract mean for all features  i.e.comp, BR
# # NOTE: MA, AA and ratios should all be handeled seperatly


# # PCA

# # %% OLD SHIT


# # '''
# # Output: custom correlograms or comparing two features either:
# #                                                                 two BR all comp or
# #                                                                             two comp all BR

# # '''
# # ### PRODUCTS TO PRESENT CUSTOM STUFF ###

# # ####CUSTOM HBRID CORRELOGRAM within BR (cuaston plot style for figure also )
# # neurotransmitters = [ 'NA', 'DA', '5HT', 'GLU', 'GABA', 'GLY']   # 'ASP', 'GLY'
# # pearson_correlations_within_BR_custom (treatment_dict, list_of_groups, exist_dict, df_inc_ratios, name = 'neurotransmitters', p_value = 0.05, compounds_to_analise = neurotransmitters )

# # #comparison between two BR of neurotranmitters
# # pearson_correlations_between_two_BR (treatment_dict, list_of_groups, exist_dict,
# #                                            df_inc_ratios, name = 'neurotransmitters_between_VL_VM', p_value = 0.05,
# #                                            compounds_to_analise = ['NA','5HT','DA', 'GLU', 'GABA', 'GLY'],
# #                                            BR_to_analise = ['VL', 'VM'])

# # ### NEED TO DEBUG

# # #   File "/Users/jasminebutler/Desktop/HPLC_module.py", line 679, in pearson_correlations_between_two_BR
# # #     coulmns  = list(map("_".join, itertools.product(compounds_to_analise, BR_to_analise)))

# # # TypeError: sequence item 0: expected str instance, list found

# # pearson_correlations_between_two_BR (treatment_dict, list_of_groups, ratio_match_ind_dict,
# #                                            df_inc_ratios, name = 'ratios_between_VL_VM', p_value = 0.05,
# #                                            compounds_to_analise = ['DOPAC_DA', '5HIAA_5HT'],
# #                                            BR_to_analise = ['VL', 'VM'])
