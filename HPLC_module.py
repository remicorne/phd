#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 12:46:03 2022

@author: jasminebutler   
"""

# %% INSTALL
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import scipy

# conda install -c conda-forge pingouin
import pingouin as pg

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
from scipy.stats import ttest_ind
# conda install -c conda-forge weasyprint
import weasyprint
import pprint
from scipy.stats import pearsonr
from statannotations.Annotator import Annotator
from statannot import add_stat_annotation
import itertools
import time as time

import numpy as np
from outliers import smirnov_grubbs as grubbs
import scipy.stats as stats

from matplotlib.pyplot import cm
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

# %% FUNCTIONS


def check_data_input(list_of_compounds, list_of_brain_regions):
    '''Print summary of data in console'''

    print('number of brain regions: ', len(list_of_brain_regions))
    print('regions include:')

    for i in list_of_brain_regions:

        print(i)

    print('number of compounds: ', len(list_of_compounds))
    print('compounds include:')

    for i in list_of_compounds:

        print(i)

    return


def create_ratios_for_each_BR(df, list_of_brain_regions, compound_ratio_dict, exist_dict, n_mice, groups):
    '''Add ratios from ratio dict to working df, create list of ratios and dict of ratios'''

    ratio_match_ind_dict = {}
    list_of_ratios = []

    for BR in list_of_brain_regions:

        ratio_match_ind_dict[BR] = {}

        for MA_num, ratio_list_denum in compound_ratio_dict.items():

            for MA_denum in ratio_list_denum:

                # need to check that the ratio can be made and if not then not included in dict
                # if MA_num and MA_denum in exist_dict[BR].keys():
                if MA_num in exist_dict[BR].keys() and MA_denum in exist_dict[BR].keys():
                    # print (MA_num, 'and', MA_denum, 'exist in', BR)

                    ratio_name = MA_num + '_' + MA_denum
                    num_ind = exist_dict[BR][MA_num]['viable_indices']
                    denum_ind = exist_dict[BR][MA_denum]['viable_indices']

                    overlap = np.intersect1d(num_ind, denum_ind)

                    if len(overlap) > 0:

                        ratio_match_ind_dict[BR][ratio_name] = {}

                        df[ratio_name, BR] = np.zeros(n_mice)
                        df[ratio_name, BR][overlap] = df[MA_num][BR][overlap] / \
                            df[MA_denum][BR][overlap]

                        non_overlapping = np.delete(np.arange(n_mice), overlap)
                        df[ratio_name, BR][non_overlapping] = np.full(
                            len(non_overlapping), np.nan)

                        ratio_match_ind_dict[BR][ratio_name]['viable_indices'] = overlap

                        group_spec_ind_dict, missing_groups = find_group_specific_indices(
                            df, overlap, groups)
                        ratio_match_ind_dict[BR][ratio_name]['group_specific_ind'] = group_spec_ind_dict
                        ratio_match_ind_dict[BR][ratio_name]['missing_groups'] = missing_groups

                        list_of_ratios.append(ratio_name)

                    else:

                        print(MA_num, 'and', MA_denum,
                              "don't have the same data")

                else:
                    print('unable to create ratio', MA_num, '_',  MA_denum)

    return df, np.unique(list_of_ratios), ratio_match_ind_dict


def remove_groups_from_df(df_inc_ratios,  groups_to_drop=[]):

    df_dose_responce = df_inc_ratios.copy()
    for i in groups_to_drop:
        # need to check here that when i drop values i am not messing up the indexes in the dictionary
        df_dose_responce.drop(
            df_dose_responce.index[df_dose_responce['group', 'no'] == i], inplace=True)

    return df_dose_responce

# def get_summary_of_exsisting_data( df, list_of_brain_regions, list_of_groups, list_of_mice):

#     headers  = df.columns[2:]
#     exist_dict = {}

#     for BR in list_of_brain_regions:
#         # print('BR = ', BR)
#         corr_compunds = [ h[0] for h in headers if h[1] == BR]

#         exist_dict [ BR ] = {}

#         for comp in corr_compunds:
#             # print('compound = ', comp)
#             exist_ind = np.where( df[ comp, BR] != 0 ) [0]     # worked fine for TCB2 couldnt detect NaN in SNORD115116
#             # exist_ind = np.where(df[comp, BR].isnull() == False)[0]      # worked somewhat for SNORD115116 but didnt allow for correct dictionary of valid indicies to be created

#             if len(exist_ind) > 0 :


#                 group_spec_dict, missing_groups = find_group_specific_indices(df, exist_ind, list_of_groups)

#                 exist_dict [ BR ] [comp]  = {
#                                             'viable_indices' : exist_ind,
#                                             'missing_groups' : missing_groups,
#                                             'group_specific_ind' : group_spec_dict
#                                             }
#             # print(list_of_mice, exist_ind)
#             if len( list_of_mice) - len( exist_ind ) > 0:
#                 print(' BR : ', BR,
#                       ' Comnpound : ', comp,
#                       ', nb miss mice:', len( list_of_mice) - len( exist_ind ),
#                       ', nb miss groups', len(missing_groups) )

#     return exist_dict


def grubbs_test(df_inc_ratios, exist_dict, list_of_groups, treatment_dict, p_value=0.05, remove_all=False):
    ''' 
    For every BR and compound Grubbs test is run, excluded data reported and removed from df_inc_ratios

    per treatment group and a single outlier per group 
    '''
    # count = 0

    # df_outliers_removed = pd.DataFrame({'Treatment': , } )
    outlier_dict = {}

    for treatment in list_of_groups:
        outlier_dict[treatment_dict[treatment]] = {}

        for BR, BR_dict in exist_dict.items():
            # count = count + 1
            # if count == 2:
            #     break
            outlier_dict[treatment_dict[treatment]][BR] = {}

            for comp, comp_dict in BR_dict.items():
                # print(comp_dict)
                print('working on ...   ', BR, comp)
                if treatment in comp_dict['missing_groups']:
                    break
                x = np.array(
                    df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][treatment]])

                if remove_all == True:
                    # removing all outliers not just one #https://pypi.org/project/outlier_utils/
                    x_ = grubbs.test(x, alpha=p_value)

                else:
                    G_calculated, outlier = grubbs_calculted_value(x)
                    G_critical = grubbs_critical_value(len(x), p_value)
                    if G_critical < G_calculated:
                        ''' we can reject the Ho: there are no outliers in the set '''
                        x_ = np.delete(x, np.where(x == outlier))
                    else:
                        x_ = x

                if len(x) != len(x_):  # if outlier detected
                    print('Grubbs Test exclusion for ', BR, ' ', comp,
                          'treatment group: ', treatment_dict[treatment])
                    print('mean of new set = ', x_.mean(),
                          'excluded point = ', np.setdiff1d(x, x_))

                    mouse_outliers = []
                    for outlier in np.setdiff1d(x, x_):
                        ind_of_outllier = df_inc_ratios.loc[df_inc_ratios[comp, BR] == float(
                            outlier)].index[0]
                        mouse_outliers.append(
                            df_inc_ratios.loc[ind_of_outllier]['mouse', 'no'])

                    outlier_dict[treatment_dict[treatment]][BR][comp] = {'mice_removed': mouse_outliers,
                                                                         'excluded_values': np.setdiff1d(x, x_),
                                                                         'mean_new_set': x_.mean()}

                    # print(ind_of_outllier)
                    print('mouse number ', mouse_outliers, ' excluded')

                    df_inc_ratios[comp, BR].replace(
                        np.setdiff1d(x, x_), np.NaN, inplace=True)

    return df_inc_ratios, outlier_dict


# https://www.youtube.com/watch?v=Hn_lMUaMcak&ab_channel=BhaveshBhatt
def grubbs_calculted_value(y):
    avg_y = np.mean(y)
    abs_val_minus_avg = abs(y - avg_y)
    max_of_deviations = max(abs_val_minus_avg)
    outlier = max_of_deviations+avg_y
    # outlier_ind = np.where(abs_val_minus_avg == max_of_deviations)
    s = np.std(y)
    G_calculated = max_of_deviations / s

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


def save_outlier_info(outlier_dict, name='compounds'):
    writer = pd.ExcelWriter('ouliers_by_treatment_' + name + '.xlsx')
    for treatment, treatment_BRs_dict in outlier_dict.items():

        df_outlier_dict = pd.DataFrame.from_dict(treatment_BRs_dict)
        treatment_replaced = treatment.replace('/', '_')
        df_outlier_dict.to_excel(writer, sheet_name=treatment_replaced)
        
    writer.save()
    # writer._save()  # TODO make compatible with jasmine code that uses 'save()' and not '_save()' DEAR REMI: i have commented this change everywhere
    writer.close()
    return

# old way!
def one_way_ANOVA_post_hoc_between_groups(df_inc_ratios, exist_dict, groups_to_drop=[], name='MA'):
    '''calculates a one way anova betweentreatment 1 ,2,3,4 outputs .csv of anova data
                    if p<0.05 Tukeys test is performed, printed and saved to an excel spreadsheet
                                                                        compound_type is for the labeling of outputed files'''
    df_dose_responce = remove_groups_from_df(
        df_inc_ratios,  groups_to_drop=groups_to_drop)
    df_one_way_anova = pd.DataFrame(data=['F', 'p'])  # create empy df
    df_one_way_anova.columns = pd.MultiIndex.from_tuples([['compound', 'BR']])

    writer = pd.ExcelWriter('TUKEY_dose_resp_' + name + '.xlsx')

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():

            # if exist_dict[BR][comp]['missing_groups']) #dealing with missing groups
            # need to debug and to add correct missing groups! SHIVA
            if len(exist_dict[BR][comp]['missing_groups'].intersection(set([1, 2, 3, 4]))) > 0:
                continue
            # if treatment in exist_dict[BR][comp]['group_specific_ind'].keys():

            F_value, p_value = scipy.stats.f_oneway(df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][1]],
                                                    df_inc_ratios[comp,
                                                                  BR][comp_dict['group_specific_ind'][2]],
                                                    df_inc_ratios[comp,
                                                                  BR][comp_dict['group_specific_ind'][3]],
                                                    df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][4]])

            df_F_p_temp = pd.DataFrame(data=[F_value, p_value])
            df_F_p_temp.columns = pd.MultiIndex.from_tuples([[comp, BR]])
            df_one_way_anova = df_one_way_anova.join(df_F_p_temp[comp, BR])

            if p_value < 0.05:  # followed by post hoc

                # could this be done through the dictionary
                mc = MultiComparison(
                    df_dose_responce[comp, BR], df_dose_responce['group', 'no'])
                mc_results = mc.tukeyhsd()
                print(mc_results)
                df_Tukey = pd.DataFrame(
                    data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
                df_Tukey.to_excel(writer, sheet_name=comp+'_'+BR)

    df_one_way_anova = df_one_way_anova.round(6)
    df_one_way_anova.to_csv(os.getcwd() + '/' + name + '.csv')
    writer.save()  # TODO make compatible with jasmine code that uses 'save()' and not '_save()'
    writer.close()

    return df_one_way_anova

def onewayANOVA_Tukeyposthoc(df_inc_ratios, exist_dict, groups_to_drop=[], name='MA'): #POS2.0 still HARD CODED GROUPS
    '''calculates a one way anova betweentreatment 1 ,2,3,4 outputs .csv of anova data
                    if p<0.05 Tukeys test is performed, printed and saved to an excel spreadsheet
                                                                        compound_type is for the labeling of outputed files'''
    df_dose_responce = remove_groups_from_df(
        df_inc_ratios,  groups_to_drop=groups_to_drop)
    df_one_way_anova = pd.DataFrame(data=['F', 'p'])  # create empy df
    df_one_way_anova.columns = pd.MultiIndex.from_tuples([['compound', 'BR']])

    writer = pd.ExcelWriter('one_way_ANOVA_ifsig_Tukey_' + name + '.xlsx')

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():


            if len(exist_dict[BR][comp]['missing_groups'].intersection(set([1, 2, 3, 4]))) > 0: #if group is missking skip
                continue

            F_value, p_value = scipy.stats.f_oneway(df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][1]],
                                                    df_inc_ratios[comp,
                                                                  BR][comp_dict['group_specific_ind'][2]],
                                                    df_inc_ratios[comp,
                                                                  BR][comp_dict['group_specific_ind'][3]],
                                                    df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][4]])

            df_F_p_temp = pd.DataFrame(data=[F_value, p_value])
            df_F_p_temp.columns = pd.MultiIndex.from_tuples([[comp, BR]])
            df_one_way_anova = df_one_way_anova.join(df_F_p_temp[comp, BR])

            if p_value < 0.05:  # followed by post hoc

                # could this be done through the dictionary
                mc = MultiComparison(
                    df_dose_responce[comp, BR], df_dose_responce['group', 'no'])
                mc_results = mc.tukeyhsd()
                print(mc_results)
                df_Tukey = pd.DataFrame(
                    data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
                df_Tukey.to_excel(writer, sheet_name=comp+'_'+BR)

    df_one_way_anova = df_one_way_anova.round(6)
    df_one_way_anova.to_excel(writer, sheet_name= name+'_one_way_ANOVA')
    # df_one_way_anova.to_csv(os.getcwd() + '/' + name + '.csv') #saving to csv randomly not nice will put in excel  *^
    
    writer.save()  # TODO make compatible with jasmine code that uses 'save()' and not '_save()'
    writer.close()

    return df_one_way_anova


def student_t_test(exist_dict, df_inc_ratios, name='compounds'):
    '''     perform student t-test between groups, 
       if significant saved.       '''

    # create df for p and t values
    df_t_test_results = pd.DataFrame(
        'NaN', index=['t_value', 'p_value'], columns=df_inc_ratios.columns[2:])

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():

            a = df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][1]]
            b = df_inc_ratios[comp, BR][comp_dict['group_specific_ind'][2]]

            t_stat, p = ttest_ind(a.dropna(), b.dropna())

            df_t_test_results[comp, BR] = t_stat, p

    df_t_test_results.to_csv(name + '_t_test_values.csv')

    return df_t_test_results


def calculate_shapiro_for_treatment(exist_dict, list_of_brain_regions, list_of_compounds, df_inc_ratios, treatment=2):

    df_shapiro = pd.DataFrame(data=['F', 'p'])  # create empty df tobe filled
    df_shapiro.columns = pd.MultiIndex.from_tuples([['compound', 'BR']])

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():     # to deal with groups that have n < 3

            if treatment in exist_dict[BR][comp]['group_specific_ind'].keys():

                x = df_inc_ratios[comp,
                                  BR][comp_dict['group_specific_ind'][treatment]]

                # if len(x) > 3:     #only do stats on groups  >3   #### if exist_dict is fixed then wont need SHIVA

                shapiro_stats = pd.DataFrame(scipy.stats.shapiro(x))
                F, p = scipy.stats.shapiro(x)
                shapiro_stats.columns = pd.MultiIndex.from_tuples([[comp, BR]])

                if p < 0.05:
                    print(comp, 'in', BR,
                          'is not normaly distributed for treatment ', treatment)
                    print('p = ', p)

                df_shapiro = df_shapiro.join(shapiro_stats[comp, BR])

    return df_shapiro


def save_shapiro_data_for_all_treatments(exist_dict, list_of_brain_regions, list_of_compounds, df_inc_ratios, list_of_groups):

    writer = pd.ExcelWriter('shapiro_data_all_treatments.xlsx')

    for treatment in list_of_groups:

        df_shapiro = calculate_shapiro_for_treatment(
            exist_dict, list_of_brain_regions, list_of_compounds, df_inc_ratios, treatment=treatment)
        treatment_str = str(treatment)
        df_shapiro.to_excel(
            writer, sheet_name='treatmenmt_group_'+treatment_str)
    # writer._save()  # TODO make compatible with jasmine code that uses 'save()' and not '_save()'
    writer.save()
    writer.close()
    return


def two_way_ANOVA(exist_dict, df_inc_ratios, treatments_to_compare=[1, 3, 5, 6], name='compounds'):
    '''create columns for the factors i.e. Agonist: True/False Antagomist: True/False
                        colapse multiindex hedders for anova'''
    df_factors = df_inc_ratios.copy()

    # structure df for two way anova- create columns for each factor ag = TCB2 and ant = MDL
    MDL_positive_index = df_factors.index[(
        df_factors['group', 'no'] == 5) | (df_factors['group', 'no'] == 6)]
    TCB2_positive_index_list = df_factors.index[(
        df_factors['group', 'no'] == 3) | (df_factors['group', 'no'] == 5)]

    bool_arr_TCB2 = np.full(len(df_factors), False)
    bool_arr_TCB2[TCB2_positive_index_list] = True
    # add boolian column with TCB2 positive = true
    df_factors['AG', 'TCB2'] = bool_arr_TCB2

    bool_arr_MDL = np.full(len(df_factors), False)
    bool_arr_MDL[MDL_positive_index] = True
    # add boolian column with MDL positive = true
    df_factors['ANT', 'MDL'] = bool_arr_MDL

    '''     perform two way anova between groups, 
        if significant saved to excel (seperate sheets) and followed by a one way anova, 
    if significant post hoc tukes is performed and saved.       '''

    writer_tukey = pd.ExcelWriter('TUKEY_ag_ant_' + name + '.xlsx')
    writer_two_way = pd.ExcelWriter(
        name + 'two_way_ANOVA_agonist_antagonist.xlsx')

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():

            # print(exist_dict[BR][comp]['missing_groups'])

            # need to debug and to add correct missing groups! SHIVA
            if len(exist_dict[BR][comp]['missing_groups'].intersection(set(treatments_to_compare))) > 0:
                print('missing group data to analise')
                continue

            # create df for analysis
            indexes_to_keep = []
            for treatment in treatments_to_compare:
                indexes_to_keep.extend(
                    list(comp_dict['group_specific_ind'][treatment]))

            d = {'dv': df_factors[comp, BR][indexes_to_keep], 'ANT_MDL': df_factors['ANT', 'MDL'][indexes_to_keep],
                  'AG_TCB2': df_factors['AG', 'TCB2'][indexes_to_keep], 'group': df_factors['group', 'no'][indexes_to_keep]}
            df_anova_working = pd.DataFrame(data=d)

            two_way_anova = pg.anova(data=df_anova_working, dv='dv', between=[
                                      ('ANT_MDL'), ('AG_TCB2')], detailed=True).round(3)

            p_list = two_way_anova['p-unc'][0:3]

            if min(p_list) < 0.05:
                print(comp, ' ', BR, ' is significant for TWO-way anova',
                      '   p=', min(p_list))
                two_way_anova.to_excel(writer_two_way, sheet_name=comp+'_'+BR)

                F_value, p_value_one_way = scipy.stats.f_oneway(df_factors[comp, BR][comp_dict['group_specific_ind'][1]],
                                                                df_factors[comp,
                                                                            BR][comp_dict['group_specific_ind'][3]],
                                                                df_factors[comp,
                                                                            BR][comp_dict['group_specific_ind'][5]],
                                                                df_factors[comp, BR][comp_dict['group_specific_ind'][6]])

                if p_value_one_way < 0.05:
                    print(comp, ' ', BR, ' is significant for ONE-way anova',
                          '   p=', p_value_one_way)

                    # Nan_free_data = df_agonist_antagonist_colapsed_header.filter(items=['GABA_VM', 'group_no']).dropna()

                    mc = MultiComparison(
                        df_anova_working['dv'], df_anova_working['group'])
                    mc_results = mc.tukeyhsd()
                    print(mc_results)
                    df_Tukey = pd.DataFrame(
                        data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
                    df_Tukey.to_excel(writer_tukey, sheet_name=comp+'_'+BR)
                else:
                    print(comp, ' ', BR, ' is NOT significant for ONE-way-anova',
                          '    p=', p_value_one_way)
            else:
                print(comp, ' ', BR, ' is NOT significant',
                      '    p=', min(p_list))

    # TODO: replace '_save()' by 'save()' for compat with jasmine's pandas version?
    writer_tukey.save()
    writer_tukey.close()
    # TODO: replace '_save()' by 'save()' for compat with jasmine's pandas version?
    writer_two_way.save()
    writer_two_way.close()
    return


def twoway_ANOVA_oneway_ANOVA_ifsig_Tukey(exist_dict, df_inc_ratios, treatments_to_compare=[1, 3, 5, 6], name='compounds'):
    '''create columns for the factors i.e. Agonist: True/False Antagomist: True/False
                        colapse multiindex hedders for anova'''
    df_factors = df_inc_ratios.copy()

    # structure df for two way anova- create columns for each factor ag = TCB2 and ant = MDL
    MDL_positive_index = df_factors.index[(
        df_factors['group', 'no'] == 5) | (df_factors['group', 'no'] == 6)]
    TCB2_positive_index_list = df_factors.index[(
        df_factors['group', 'no'] == 3) | (df_factors['group', 'no'] == 5)]

    bool_arr_TCB2 = np.full(len(df_factors), False)
    bool_arr_TCB2[TCB2_positive_index_list] = True
    # add boolian column with TCB2 positive = true
    df_factors['AG', 'TCB2'] = bool_arr_TCB2

    bool_arr_MDL = np.full(len(df_factors), False)
    bool_arr_MDL[MDL_positive_index] = True
    # add boolian column with MDL positive = true
    df_factors['ANT', 'MDL'] = bool_arr_MDL

    '''     perform two way anova between groups, 
       if significant saved to excel (seperate sheets) and followed by a one way anova, 
    if significant post hoc tukes is performed and saved.       '''

    writer_tukey = pd.ExcelWriter('oneway_ANOVA_Tukey_' + name + '.xlsx')
    df_one_way_anova = pd.DataFrame(data=['F', 'p'])  # create empy df
    df_one_way_anova.columns = pd.MultiIndex.from_tuples([['compound', 'BR']])
    
    writer_two_way = pd.ExcelWriter(name + '_two_way_ANOVA.xlsx')

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():

            # print(exist_dict[BR][comp]['missing_groups'])

            # need to debug and to add correct missing groups! SHIVA
            if len(exist_dict[BR][comp]['missing_groups'].intersection(set(treatments_to_compare))) > 0:
                print('missing group data to analise')
                continue

            # create df for analysis
            indexes_to_keep = []
            for treatment in treatments_to_compare:
                indexes_to_keep.extend(
                    list(comp_dict['group_specific_ind'][treatment]))

            d = {'dv': df_factors[comp, BR][indexes_to_keep], 'ANT_MDL': df_factors['ANT', 'MDL'][indexes_to_keep],
                 'AG_TCB2': df_factors['AG', 'TCB2'][indexes_to_keep], 'group': df_factors['group', 'no'][indexes_to_keep]}
            df_anova_working = pd.DataFrame(data=d)
            # do two way ANOVA
            two_way_anova = pg.anova(data=df_anova_working, dv='dv', between=[
                                     ('ANT_MDL'), ('AG_TCB2')], detailed=True).round(3)

            p_list = two_way_anova['p-unc'][0:3]
            
            two_way_anova.to_excel(writer_two_way, sheet_name=comp+'_'+BR) #save all two way ANOVA
            
            #do one way ANOVA for all not just for sig two way ANOVA
            F_value, p_value_one_way = scipy.stats.f_oneway(df_factors[comp, BR][comp_dict['group_specific_ind'][1]],
                                                                df_factors[comp,
                                                                           BR][comp_dict['group_specific_ind'][3]],
                                                                df_factors[comp,
                                                                           BR][comp_dict['group_specific_ind'][5]],
                                                                df_factors[comp, BR][comp_dict['group_specific_ind'][6]])
            
            df_F_p_temp = pd.DataFrame(data=[F_value, p_value_one_way])
            df_F_p_temp.columns = pd.MultiIndex.from_tuples([[comp, BR]])
            df_one_way_anova = df_one_way_anova.join(df_F_p_temp[comp, BR])
            
            if p_value_one_way < 0.05: #if oneway ANOVA significant 
                print(comp, ' ', BR, ' is significant for ONE-way anova',
                      '   p=', p_value_one_way)

                # Nan_free_data = df_agonist_antagonist_colapsed_header.filter(items=['GABA_VM', 'group_no']).dropna()

                mc = MultiComparison(
                    df_anova_working['dv'], df_anova_working['group'])
                mc_results = mc.tukeyhsd()
                print(mc_results)
                df_Tukey = pd.DataFrame(
                    data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
                df_Tukey.to_excel(writer_tukey, sheet_name=comp+'_'+BR)
            else:
                print(comp, ' ', BR, ' is NOT significant for ONE-way-anova',
                      '    p=', p_value_one_way)
            
            if min(p_list) < 0.05:
                print(comp, ' ', BR, ' is significant for TWO-way anova',
                      '   p=', min(p_list))
                # two_way_anova.to_excel(writer_two_way, sheet_name=comp+'_'+BR) #only save significant two way ANOVA

                # F_value, p_value_one_way = scipy.stats.f_oneway(df_factors[comp, BR][comp_dict['group_specific_ind'][1]],
                #                                                 df_factors[comp,
                #                                                            BR][comp_dict['group_specific_ind'][3]],
                #                                                 df_factors[comp,
                #                                                            BR][comp_dict['group_specific_ind'][5]],
                #                                                 df_factors[comp, BR][comp_dict['group_specific_ind'][6]])

                # if p_value_one_way < 0.05:
                #     print(comp, ' ', BR, ' is significant for ONE-way anova',
                #           '   p=', p_value_one_way)

                #     # Nan_free_data = df_agonist_antagonist_colapsed_header.filter(items=['GABA_VM', 'group_no']).dropna()

                #     mc = MultiComparison(
                #         df_anova_working['dv'], df_anova_working['group'])
                #     mc_results = mc.tukeyhsd()
                #     print(mc_results)
                #     df_Tukey = pd.DataFrame(
                #         data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
                #     df_Tukey.to_excel(writer_tukey, sheet_name=comp+'_'+BR)
                # else:
                #     print(comp, ' ', BR, ' is NOT significant for ONE-way-anova',
                #           '    p=', p_value_one_way)
            else:
                print(comp, ' ', BR, ' is NOT significant',
                      '    p=', min(p_list))
                
    df_one_way_anova = df_one_way_anova.round(6)
    df_one_way_anova.to_excel(writer_tukey, sheet_name= name+'_one_way_ANOVA')
    # TODO: replace '_save()' by 'save()' for compat with jasmine's pandas version?
    writer_tukey.save()
    writer_tukey.close()
    # TODO: replace '_save()' by 'save()' for compat with jasmine's pandas version?
    writer_two_way.save()
    writer_two_way.close()
    return


def two_way_ANOVA_ALZ(exist_dict, df_inc_ratios, treatments_to_compare=[1, 2, 3, 4], name='compounds'):
    '''create columns for the factors i.e. Agonist: True/False Antagomist: True/False
                        colapse multiindex hedders for anova'''
    df_factors = df_inc_ratios.copy()

    # structure df for two way anova- create columns for each factor ag = TCB2 and ant = MDL
    # structure df for two way anova- create columns for each factor old = TCB2 and AD = MDL
    MDL_positive_index = df_factors.index[(df_factors['group', 'no'] == 2) | (
        df_factors['group', 'no'] == 4)]  # AD positive (false for WT)
    TCB2_positive_index_list = df_factors.index[(df_factors['group', 'no'] == 3) | (
        df_factors['group', 'no'] == 4)]  # old positive (false for young)

    bool_arr_TCB2 = np.full(len(df_factors), False)
    bool_arr_TCB2[TCB2_positive_index_list] = True
    # add boolian column with TCB2 positive = true
    df_factors['AG', 'TCB2'] = bool_arr_TCB2

    bool_arr_MDL = np.full(len(df_factors), False)
    bool_arr_MDL[MDL_positive_index] = True
    # add boolian column with MDL positive = true
    df_factors['ANT', 'MDL'] = bool_arr_MDL

    '''     perform two way anova between groups, 
       if significant saved to excel (seperate sheets) and followed by a one way anova, 
    if significant post hoc tukes is performed and saved.       '''

    writer_tukey = pd.ExcelWriter('TUKEY_ag_ant_' + name + '.xlsx')
    writer_two_way = pd.ExcelWriter(
        name + 'two_way_ANOVA_agonist_antagonist.xlsx')

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():

            # print(exist_dict[BR][comp]['missing_groups'])

            # need to debug and to add correct missing groups! SHIVA
            if len(exist_dict[BR][comp]['missing_groups'].intersection(set(treatments_to_compare))) > 0:
                print('missing group data to analise')
                continue

            # create df for analysis
            indexes_to_keep = []
            for treatment in treatments_to_compare:
                indexes_to_keep.extend(
                    list(comp_dict['group_specific_ind'][treatment]))

            d = {'dv': df_factors[comp, BR][indexes_to_keep], 'ANT_MDL': df_factors['ANT', 'MDL'][indexes_to_keep],
                 'AG_TCB2': df_factors['AG', 'TCB2'][indexes_to_keep], 'group': df_factors['group', 'no'][indexes_to_keep]}
            df_anova_working = pd.DataFrame(data=d)

            two_way_anova = pg.anova(data=df_anova_working, dv='dv', between=[
                                     ('ANT_MDL'), ('AG_TCB2')], detailed=True).round(3)

            p_list = two_way_anova['p-unc'][0:3]

            if min(p_list) < 0.05:
                print(comp, ' ', BR, ' is significant for TWO-way anova',
                      '   p=', min(p_list))
                two_way_anova.to_excel(writer_two_way, sheet_name=comp+'_'+BR)

                F_value, p_value_one_way = scipy.stats.f_oneway(df_factors[comp, BR][comp_dict['group_specific_ind'][1]],
                                                                df_factors[comp,
                                                                           BR][comp_dict['group_specific_ind'][2]],
                                                                df_factors[comp,
                                                                           BR][comp_dict['group_specific_ind'][3]],
                                                                df_factors[comp, BR][comp_dict['group_specific_ind'][4]])

                if p_value_one_way < 0.05:
                    print(comp, ' ', BR, ' is significant for ONE-way anova',
                          '   p=', p_value_one_way)

                    # Nan_free_data = df_agonist_antagonist_colapsed_header.filter(items=['GABA_VM', 'group_no']).dropna()

                    mc = MultiComparison(
                        df_anova_working['dv'], df_anova_working['group'])
                    mc_results = mc.tukeyhsd()
                    print(mc_results)
                    df_Tukey = pd.DataFrame(
                        data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
                    df_Tukey.to_excel(writer_tukey, sheet_name=comp+'_'+BR)
                else:
                    print(comp, ' ', BR, ' is NOT significant for ONE-way-anova',
                          '    p=', p_value_one_way)
            else:
                print(comp, ' ', BR, ' is NOT significant',
                      '    p=', min(p_list))

    # TODO: replace '_save()' by 'save()' for compat with jasmine's pandas version?
    writer_tukey._save()
    writer_tukey.close()
    # TODO: replace '_save()' by 'save()' for compat with jasmine's pandas version?
    writer_two_way._save()
    writer_two_way.close()
    return


def plot_hist_comparing_treatment_SEM(treatment_dict, exist_dict,  palette_labeled, df_compound_mean, df_compound_SEM,      # note: if there is a missing group in the set there will be no plot and no stats! SHIVA
                                      treatments_to_compare=[1, 3, 5, 6], name='name me',
                                      order=['vehicles', '3mg/kgTCB',
                                             '0.2mg/kgMDL', 'TCB+MDL'],
                                      test_path=os.getcwd() + '/TUKEY_ag_ant_compounds.xlsx',
                                      ratio=False):
    '''


order of teatments to compare is the order of plotting 

    '''

    # colors = ['white', 'red', 'grey', 'black']
    # sns.set_palette(sns.color_palette(colors))
    # palette_labeled ={'vehicles': "white", '10mg/kgTCB': "firebrick", '3mg/kgTCB': "red", '0.3mg/kgTCB': "salmon",'TCB+MDL': "grey",'0.2mg/kgMDL': "black"}

    pdf = PdfPages(name+'_histograms_.pdf')

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():
            print(comp + '_' + BR)

            if len(exist_dict[BR][comp]['missing_groups'].intersection(set(treatments_to_compare))) > 0:
                print('missing treatment group ',
                      exist_dict[BR][comp]['missing_groups'])
                continue

            mean_treatment = []  # CTE's

            SEM_treatment = []  # error
            labels = []
            color = []

            for treatment in treatments_to_compare:

                mean_treatment.append(np.mean(
                    df_compound_mean[df_compound_mean['group', 'no'] == treatment][comp, BR]))
                SEM_treatment.append(
                    np.mean(df_compound_SEM[df_compound_SEM['group', 'no'] == treatment][comp, BR]))
                labels.append(treatment_dict[treatment])
                color.append(palette_labeled[treatment_dict[treatment]])

            x_pos = np.arange(len(labels))

            fig, ax = plt.subplots(figsize=(20, 10))

            ax.bar(x_pos, np.array(mean_treatment),
                   yerr=SEM_treatment, color=color, edgecolor='black',
                   align='center',
                   alpha=0.5,
                   ecolor='black',
                   capsize=10)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylabel("ng/mg of tissue", fontsize=24)
            ax.tick_params(axis='y', which='major', labelsize=20)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, fontsize=20)
            ax.set_title(comp.replace('_', '/') +
                         ' in ' + BR, y=1.04, fontsize=34)
            # ax.yaxis.grid(True)
            if ratio == True:
                ax.set_ylabel(" ", fontsize=24)

            # ax = put_significnce_stars(ax, treatment_dict, test_path, data = df_working_labeled,
            #                            x = "group", y = 'dv',
            #                            order = order, sheet = comp + '_' + BR)#, p_values)

            # ax.tick_params(labelsize=24)

            ax.set_xlabel(" ", fontsize=20)  # treatments
            ax.set_title(comp.replace('_', '/') + ' in ' + BR,
                         y=1.04, fontsize=34)  # '+/- 68%CI'
            pdf.savefig(fig)
            plt.close(fig)

    pdf.close()
    plt.close('all')
    return


def plot_hist_comparing_treatment_CI(treatment_dict, exist_dict, df_inc_ratios, palette_labeled,
                                     treatments_to_compare=[1, 2, 3, 4], name='Dose_Responcxe_compounds',
                                     order=['vehicles', '0.3mg/kgTCB',
                                            '3mg/kgTCB', '10mg/kgTCB'],
                                     test_path=os.getcwd() + '/TUKEY_dose_resp_compounds.xlsx', ratio=False, mouse_id=False):

    # def plot_hist_comparing_treatment_CI (treatment_dict, exist_dict, df_inc_ratios, palette_labeled,       # note: if there is a missing group in the set there will be no plot and no stats! SHIVA
    #                                    treatments_to_compare = [1,3,5,6], name = 'name me',
    #                                    order = ['vehicles', '3mg/kgTCB','0.2mg/kgMDL','TCB+MDL' ],
    #                                    test_path = '/Users/jasminebutler/Desktop/dose_responce_tukey_tables_ratios.xlsx',
    #                                    ratio = False, mouse_id = False):

    # colors = ['white', 'red', 'grey', 'black']
    # sns.set_palette(sns.color_palette(colors))
    # palette_labeled ={'vehicles': "white", '10mg/kgTCB': "firebrick", '3mg/kgTCB': "red", '0.3mg/kgTCB': "salmon",'TCB+MDL': "grey",'0.2mg/kgMDL': "black"}

    pdf = PdfPages(name+'_histograms_.pdf')

    # count = 0

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():
            # count = count + 1
            print(comp + '_' + BR)
            # if exist_dict[BR][comp]['missing_groups']
            # if len(  exist_dict[BR][comp]['missing_groups'].intersection(set([1,3,5,6]))  ) > 0:###### need to debug and to add correct missing groups! SHIVA
            # if one of the treatment groups is missing skip entire plot for that compound and brain region
            if len(exist_dict[BR][comp]['missing_groups'].intersection(set(treatments_to_compare))) > 0:

                continue

            indexes_to_keep = []

            for treatment in treatments_to_compare:
                indexes_to_keep.extend(
                    list(comp_dict['group_specific_ind'][treatment]))

            # create relevant df
            d = {'dv': df_inc_ratios[comp, BR][indexes_to_keep], 'group': df_inc_ratios['group',
                                                                                        'no'][indexes_to_keep], 'mouse': df_inc_ratios['mouse', 'no']}
            df_working = pd.DataFrame(data=d)
            df_working_labeled = df_working.replace(
                {'group': treatment_dict})  # label groups

            # col_name = str(comp)+'_'+str(BR)

            fig, ax = plt.subplots(figsize=(20, 10))
            ax = sns.barplot(x="group", y='dv', data=df_working_labeled, palette=palette_labeled,
                             ci=68, order=order, capsize=.1,
                             alpha=0.8, errcolor=".2", edgecolor=".2")

            if mouse_id == True:
                color = iter(cm.rainbow(np.linspace(0, 1, len(df_inc_ratios))))
                ax = sns.swarmplot(x="group", y="dv",  order=order, palette=color,
                                   data=df_working_labeled, hue='mouse', edgecolor='k', linewidth=1, linestyle='-')
                # https://seaborn.pydata.org/tutorial/color_palettes.html
                ax.legend_.remove()
                plt.legend(bbox_to_anchor=(0.99, 1),
                           loc='upper left', borderaxespad=0, ncol=2)

            else:

                ax = sns.swarmplot(x="group", y="dv", palette=palette_labeled, order=order,
                                   data=df_working_labeled,  edgecolor='k', linewidth=1, linestyle='-')

            ax = put_significnce_stars(ax, treatment_dict, test_path, data=df_working_labeled,
                                       x="group", y='dv',
                                       order=order, sheet=comp + '_' + BR)  # , p_values)
            # ci=95 is ~1.96 standard errors, and ci=68 is ~1 standard error for normaly distriibuted data
            ax.tick_params(labelsize=24)
            # ax.set(xlabel='Treatment', ylabel='ng/mg of tissue')
            ax.set_ylabel("ng/mg of tissue", fontsize=24)

            if ratio == True:
                ax.set_ylabel(" ", fontsize=24)

            ax.set_xlabel(" ", fontsize=20)  # treatments
            ax.set_title(comp.replace('_', '/') + ' in ' + BR,
                         y=1.04, fontsize=34)  # '+/- 68%CI'
            sns.despine(left=False)

            # if count >  1:
            #     break

            pdf.savefig(fig)
            # plt.cla()
            plt.close(fig)

    pdf.close()
    plt.close('all')
    return


def put_significnce_stars(ax, treatment_dict, test_path, data=None,
                          x=None, y=None, order=None,
                          sheet='5HT_DL'):  # , p_values):

    try:
        df = pd.read_excel(test_path, sheet_name=sheet, header=0)
        df_significant = df[df['reject']]  # only take the significant rows

        if len(df_significant.index) > 0:
            print(df_significant)
            p_values = df_significant['p-adj'].values
            pairs = [(treatment_dict[i[1]["group1"]], treatment_dict[i[1]["group2"]])
                     for i in df_significant.iterrows()]

            annotator = Annotator(ax, pairs, data=data,
                                  x=x, y=y, order=order)
            # https://stackoverflow.com/questions/64081570/matplotlib-marker-annotation-fontsize-not-shrinking-below-1pt-in-pdf
            annotator.configure(text_format="star",
                                loc="inside", fontsize='xx-large')
            annotator.set_pvalues_and_annotate(p_values)

        return ax

    except ValueError:

        return ax


def plot_two_groups(treatment_dict, exist_dict, df_inc_ratios, palette_labeled, df_t_test_results,
                    name='tryplot', treatments_to_compare=[1, 2], p=0.05):

    pdf = PdfPages(name+'_histograms_.pdf')

    # count = 0
    for BR, BR_dict in exist_dict.items():
        # if count == 16:

        #                 break

        for comp, comp_dict in BR_dict.items():

            print('working on . . . ', comp, '  in  ', BR)
            indexes_to_keep = []

            # #create df to plot
            # df_to_plot = pd.DataFrame('NaN', index= np.arange(2*len( comp_dict['group_specific_ind'][1] )), columns=['group', 'dv'])

            for treatment in treatments_to_compare:

                # collect all relevant indicies
                indexes_to_keep.extend(
                    list(comp_dict['group_specific_ind'][treatment]))

            if len(indexes_to_keep) < 6:
                print('DICT ERROR / insufficient data for ', comp, 'in the ', BR)
                continue

            # create relevant df
            d = {'dv': df_inc_ratios[comp, BR][indexes_to_keep],
                 'group': df_inc_ratios['group', 'no'][indexes_to_keep]}
            df_working = pd.DataFrame(data=d)

            # print (df_working)
            # count = count +1

            df_working_labeled = df_working.replace(
                {'group': treatment_dict})  # label groups
            # df_working_labeled = df_working_labeled.dropna(axis=1)

            box_pairs = [('WT', 'SNORD_115116_KO')]

            fig, ax = plt.subplots(figsize=(20, 10))
            sns.barplot(x="group", y="dv", data=df_working_labeled, ci=95,
                        palette=palette_labeled, linewidth=1.5, errcolor=".2", edgecolor=".2")

            ax.set_title(comp + ' in ' + BR, y=1.04, fontsize=34)
            ax.tick_params(labelsize=24)
            ax.set_ylabel("pg/mg of tissue", fontsize=24)
            ax.set_xlabel(" ", fontsize=20)
            sns.despine(left=False)

            if df_t_test_results.loc['p_value', (comp, BR)] < p:

                test_results = add_stat_annotation(ax, data=df_working_labeled, x="group", y="dv",
                                                   box_pairs=box_pairs,
                                                   # 'simple' reports p value or 'star' for stars only
                                                   test='t-test_ind', text_format='star',
                                                   loc='inside', verbose=2)
                test_results
                plt.show()

            pdf.savefig(fig)

    pdf.close()
    plt.close('all')

    return

# %% CORELATION FUNCTIONS


def calculate_pvalues(df_corr_working):
    '''Calculates p values from a pearson correlation
        output is a square df of p values comparing column to column of origional df'''

    dfcols = pd.DataFrame(columns=df_corr_working.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')  # creating empty df

    for r in df_corr_working.columns:
        for c in df_corr_working.columns:

            # count Nans in list, if not 0 then cleaning needed
            Nan_r = df_corr_working[r].isnull().sum()
            Nan_c = df_corr_working[c].isnull().sum()

            if Nan_r or Nan_c > 0:
                print('unequal length data, trimming nans')
                print('NaN count: ', Nan_r, ' =  ',
                      r, ' and ', Nan_c, ' = ', c)

                index_nan_r = df_corr_working[df_corr_working[r].isnull(
                )].index.tolist()
                clean_r = df_corr_working[r].drop(index=index_nan_r)
                clean_c = df_corr_working[c].drop(index=index_nan_r)
                Nan_clean_c = clean_c.isnull().sum()

                if Nan_clean_c > 0:  # if we still have Nans in the other column
                    index_nan_c = df_corr_working[df_corr_working[c].isnull(
                    )].index.tolist()
                    # only deleting indexes that have not been removed
                    index_nan_c_no_duplicates = [
                        x for x in index_nan_c if x not in index_nan_r]
                    clean_r = clean_r.drop(index=index_nan_c_no_duplicates)
                    clean_c = clean_c.drop(index=index_nan_c_no_duplicates)
                    print('both trimmed, length = ', len(clean_r),)

                else:
                    print('trimed length = ', len(clean_r))

            else:
                clean_r = df_corr_working[r]
                clean_c = df_corr_working[c]

            df_corr_working[r]
            pvalues[r][c] = scipy.stats.pearsonr(clean_r, clean_c)[1]

    return pvalues


def pearson_correlations_within_BR(list_of_groups, exist_dict, df_inc_ratios, treatment_dict, name='compound', p_value=0.05,
                                   column_order=['NA', 'DA', 'DOPAC', 'HVA', '3MT', '5HT', '5HI'], n_minimum=8, method='pearson'):
    '''
    Input:
            dict, column_order (can be used to select specific subset of data) and df, 
            method = ‘pearson’, ‘kendall’, ‘spearman’
    Output:
            save pearson correlogram with all significant R vaklues displayed for each BR in pdf (ratios and compounds seperatly)

    '''

    pdf = PdfPages(name + '_correlograms_within_BR.pdf')

# SIVA this df will have some NaN values and i dont really understand how it is being delt with with pandas .corr
    # https://stackoverflow.com/questions/57155427/how-does-corr-remove-na-and-null-values

    for treatment in list_of_groups:  # create df for correlations of smetric data : each collumn is a factor to correlate maintining indicies

        # treatment = 1 # to test
        # count = 0

        for BR, BR_dict in exist_dict.items():
            print('BR = ', BR)
            # if count == 1:
            #         break

            # get all indexes of treatment group
            index_of_group = df_inc_ratios.index[df_inc_ratios['group', 'no'] == treatment].tolist(
            )

            if column_order == False:

                print('no column order')

                # include only columns(compounds) without missing groups sufficient n in group/treatment
                all_cols = [x for x in BR_dict.keys() if treatment not in list(
                    exist_dict[BR][x]['missing_groups'])]
                # transform cols to only include columns with n > n_minimum
                cols = [x for x in all_cols if len(
                    exist_dict[BR][x]['group_specific_ind'][treatment]) > n_minimum]

                # unordered but works for TCB2 and SNORD no column_order
                df_corr_working = pd.DataFrame(
                    'NaN', index=index_of_group, columns=cols)

                # df_corr_working = pd.DataFrame('NaN', index=index_of_group, columns=list(BR_dict.keys())) #unordered but works for TCB2 and SNORD no column_order

                # if you dont have a column order with all of your interest compounds
                if 'HSER' in list(BR_dict.keys()):

                    df_corr_working = df_corr_working.drop(
                        columns='HSER')  # remove HSER

            else:
                print('column ordered:  ', column_order)

                # include only columns without missing groups sufficient n in group/treatment
                all_cols = [x for x in BR_dict.keys() if treatment not in list(
                    exist_dict[BR][x]['missing_groups'])]
                # transform cols to only include columns with n > n_minimum
                viable_cols = [x for x in all_cols if len(
                    exist_dict[BR][x]['group_specific_ind'][treatment]) > n_minimum]

                # list of BR in order
                cols = [x for x in column_order if x in viable_cols]

                # create df to be filled for single BR and treatment
                df_corr_working = pd.DataFrame(
                    'NaN', index=index_of_group, columns=cols)

            for comp, comp_dict in BR_dict.items():
                # print(comp_dict)

                if comp == 'HSER':  # remove HSER
                    continue

                # if treatment in list(exist_dict[BR][comp]['missing_groups']): #is the treatment group  missing for that compound, skip
                #     continue
                # if len(exist_dict[BR][comp]['group_specific_ind'][treatment]) > n_minimum: #only fill df if n>n_minimum

                # only fill columns tat have been created
                if comp in list(df_corr_working.columns):

                    df_corr_working[comp] = df_inc_ratios[comp,
                                                          BR][comp_dict['group_specific_ind'][treatment]]     # fill df
                else:
                    continue
            # count = count +1

            corr_matrix_BR = df_corr_working.corr(method=method, min_periods=1)

            fig, ax = plt.subplots(figsize=(16, 10))

            # calculate p_values and set maske to P>0.05
            p_value_matrix = calculate_pvalues(df_corr_working)
            mask = np.invert(np.tril(p_value_matrix < p_value))

            heatmap = sns.heatmap(corr_matrix_BR, vmin=-1, vmax=1,
                                  annot=True, cmap='BrBG', mask=mask, annot_kws={"size": 8})
            heatmap.set_title(
                BR + ' in ' + treatment_dict[treatment], fontdict={'fontsize': 18}, pad=12)
            ax.set_xticklabels(
                ax.get_xticklabels(),
                rotation=45,
                horizontalalignment='right')
            ax.set_ylabel('')
            ax.set_xlabel('')
            pdf.savefig(fig)
            plt.close(fig)

    pdf.close()
    plt.close('all')
    return


def pearson_correlations_within_compound(list_of_groups, exist_dict, df_inc_ratios, list_of_compounds, treatment_dict, name='compound', p_value=0.05, n_minimum=8,
                                         column_order=['OF', 'PL', 'Cg', 'Nac', 'VM', 'DM', 'VL', 'DL', 'SP', 'AM', 'HD', 'HV', 'Th', 'DY',  'VY',  'SN', 'VTA', 'DR',  'MR'], method='pearson'):
    '''Input: 

        dict, df 
        method = ‘pearson’, ‘kendall’, ‘spearman’

        Output: 

            df for one correlogram (BR between compounds or compounds between BR), 
                    create correlogram with all significant R vaklues displayed for each compound/ratio 
            --> saved to pdf '''

    pdf = PdfPages(name + '_correlograms_between_BRs.pdf')

    for treatment in list_of_groups:  # create df for correlations of smetric data : each collumn is a factor to correlate maintining indicies
        print('treatment = ', treatment)

        for comp in list_of_compounds:
            print('comp =', comp)
            # get all indexes of treatment group
            index_of_group = df_inc_ratios.index[df_inc_ratios['group', 'no'] == treatment].tolist(
            )

            existing_BR_for_comp = [
                x for x in exist_dict.keys() if comp in list(exist_dict[x].keys())]
            print('existing BR for comp = ', existing_BR_for_comp)

            all_cols = [x for x in existing_BR_for_comp if treatment not in list(
                exist_dict[x][comp]['missing_groups'])]  # include only columns(BR's) without missing group/treatment
            # transform cols to only include columns with n > n_minimum
            cols = [x for x in all_cols if len(
                exist_dict[x][comp]['group_specific_ind'][treatment]) > n_minimum]
            # print('BRs used for correlation = ', cols)

            if column_order == False:
                print('no column order')

                df_corr_working = pd.DataFrame(
                    'NaN', index=index_of_group, columns=cols)

                # if you dont have a column order with all of your interest compounds
                if 'HSER' in list(exist_dict.keys()):

                    df_corr_working = df_corr_working.drop(
                        columns='HSER')  # remove HSER

            else:
                print('column ordered:    ', column_order)

                # list of viable BR in order (compounds not in order will not be included )
                cols_ordered = [x for x in column_order if x in list(cols)]

                # create df to be filled for single BR and treatment
                df_corr_working = pd.DataFrame(
                    'NaN', index=index_of_group, columns=cols_ordered)

            for BR, BR_dict in exist_dict.items():
                # for BR, BR_dict in cols:

                if BR in list(df_corr_working.columns):
                    df_corr_working[BR] = df_inc_ratios[comp,
                                                        BR][BR_dict[comp]['group_specific_ind'][treatment]]

                else:
                    print(BR, ' does not have ', comp)

            corr_matrix_BR = df_corr_working.corr(method=method, min_periods=1)

            fig, ax = plt.subplots(figsize=(16, 10))

            # calculate p_values and set maske to P>0.05
            p_value_matrix = calculate_pvalues(df_corr_working)
            # removing columns that have all 'NaN' vslues to match correlation matrix
            p_value_matrix = p_value_matrix.dropna(axis=1, how='all')
            # removing rows that have all 'NaN' vslues
            p_value_matrix = p_value_matrix.dropna(axis=0, how='all')

            mask = np.invert(np.tril(p_value_matrix < p_value))

            heatmap = sns.heatmap(corr_matrix_BR, vmin=-1, vmax=1,
                                  annot=True, cmap='BrBG', mask=mask, annot_kws={"size": 8})
            heatmap.set_title(
                comp+'  in  ' + treatment_dict[treatment], fontdict={'fontsize': 18}, pad=12)
            ax.set_xticklabels(
                ax.get_xticklabels(),
                rotation=45,
                horizontalalignment='right')
            ax.set_ylabel('')
            ax.set_xlabel('')
            pdf.savefig(fig)
            plt.close(fig)

    pdf.close()
    plt.close('all')
    return


def pearson_correlations_between_two_BR(treatment_dict, list_of_groups, exist_dict,
                                        df_inc_ratios, name='neurotransmitters_between_VL_VM', p_value=0.05,
                                        compounds_to_analise=[
                                            'DA', '5HT', 'NA'],
                                        BR_to_analise=['VL', 'VM']):
    '''input dict and df, 
            produce  df for one correlogram (BR between compounds or compounds between BR), 
                    create correlogram with all significant R vaklues displayed for each compound/ratio 
            --> saved to pdf '''

    pdf = PdfPages(name + '_custom_correlograms.pdf')

    for treatment in list_of_groups:  # create df for correlations of smetric data : each collumn is a factor to correlate maintining indicies

        # treatment = 1 # to test
        # count = 0
        print(treatment)

        index_of_group = df_inc_ratios.index[df_inc_ratios['group', 'no'] == treatment].tolist(
        )
        # get all indexes of treatment group

        coulmns = list(map("_".join, itertools.product(
            compounds_to_analise, BR_to_analise)))

        df_corr_working = pd.DataFrame(
            'NaN', index=index_of_group, columns=coulmns)

        for BR, BR_dict in exist_dict.items():
            if BR in BR_to_analise:

                # if count == 1:
                #         break
                # create df to be filled for single BR and treatment
                for comp, comp_dict in BR_dict.items():
                    if comp in compounds_to_analise:

                        df_corr_working[comp+'_'+BR] = df_inc_ratios[comp,
                                                                     BR][comp_dict['group_specific_ind'][treatment]]     # fill df
                    else:
                        continue

            # count = count +1

        corr_matrix_BR = df_corr_working.corr(method='pearson', min_periods=1)

        # calculate p_values and set maske to P>0.05
        p_value_matrix = calculate_pvalues(df_corr_working)

        # remove all of 1 BR from rows and other BR from columns

        # searching column names as aways identical but drops rows with value
        corr_matrix_BR = corr_matrix_BR.drop(
            list(corr_matrix_BR.filter(regex='VL').columns), axis=0)
        corr_matrix_BR_plot = corr_matrix_BR.drop(list(corr_matrix_BR.filter(
            regex='VM').columns), axis=1)  # dropping columns with value

        # searching column names as aways identical but drops rows with value
        p_value_matrix = p_value_matrix.drop(
            list(p_value_matrix.filter(regex='VL').columns), axis=0)
        p_value_matrix = p_value_matrix.drop(list(p_value_matrix.filter(
            regex='VM').columns), axis=1)  # dropping columns with value

        mask = np.invert((p_value_matrix < p_value))  # np.tril

        fig, ax = plt.subplots(figsize=(16, 10))
        heatmap = sns.heatmap(corr_matrix_BR_plot, vmin=-1, vmax=1,
                              annot=True, cmap='BrBG', mask=mask, annot_kws={"size": 16})
        # use matplotlib.colorbar.Colorbar object
        cbar = ax.collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=16)
        heatmap.set_title(BR_to_analise[0]+' and '+BR_to_analise[1]+'  '+str(
            treatment_dict[treatment]), fontdict={'fontsize': 36}, pad=12, y=1.03)
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            horizontalalignment='right', fontsize=18)
        ax.set_yticklabels(
            ax.get_yticklabels(), rotation=45, fontsize=18)
        ax.set_ylabel('')
        ax.set_xlabel('')
        pdf.savefig(fig)
        plt.close('all')

    pdf.close()
    plt.close('all')
    return


def find_group_specific_indices(df, exist_ind, list_of_groups):
    ''' segregate viable sample indices by  groups'''

    # exist_group = np.unique( df[ 'group', 'no'][exist_ind])   # missing groups are groups with n = 0

    groups = df['group', 'no'].values[exist_ind]

    sample_counts = np.bincount(groups)
    ii = np.nonzero(sample_counts)[0]
    group_samples = zip(ii, sample_counts[ii])  # groups that exist
    # groups exist when n < 3
    exist_group = np.array([gr_s[0] for gr_s in group_samples if gr_s[1] > 2])

    # groups present i.e. list of groups *not groups
    missing_groups = set(list_of_groups) - set(exist_group)

    # group_spec_dict = { gr : exist_ind [ df[ 'group', 'no'][exist_ind] == gr]      # including viabile indicies with n = 1 per group
    #                     for gr in groups}

    group_spec_dict = {gr: exist_ind[df['group', 'no'][exist_ind] == gr]
                       for gr in groups
                       if len(exist_ind[df['group', 'no'][exist_ind] == gr]) > 2}   # will not count groups with n < 3 as viable indicies

    return group_spec_dict, missing_groups


def get_summary_of_exsisting_data(df, list_of_brain_regions, list_of_groups, list_of_mice):

    headers = df.columns[2:]
    exist_dict = {}

    for BR in list_of_brain_regions:
        # print('BR = ', BR)
        corr_compunds = [h[0] for h in headers if h[1] == BR]

        exist_dict[BR] = {}

        for comp in corr_compunds:
            # print('compound = ', comp)
            # worked fine for TCB2 couldnt detect NaN in SNORD115116
            exist_ind = np.where(df[comp, BR] != 0)[0]
            # exist_ind = np.where(df[comp, BR].isnull() == False)[0]      # worked somewhat for SNORD115116 but didnt allow for correct dictionary of valid indicies to be created

            if len(exist_ind) > 0:

                group_spec_dict, missing_groups = find_group_specific_indices(
                    df, exist_ind, list_of_groups)

                exist_dict[BR][comp] = {
                    'viable_indices': exist_ind,
                    'missing_groups': missing_groups,
                    'group_specific_ind': group_spec_dict
                }
            # print(list_of_mice, exist_ind)
            if len(list_of_mice) - len(exist_ind) > 0:
                print(' BR : ', BR,
                      ' Comnpound : ', comp,
                      ', nb miss mice:', len(list_of_mice) - len(exist_ind),
                      ', nb miss groups', len(missing_groups))

    return exist_dict

# DO I EVEN NEED THIS FUNCTION


def get_summary_of_exsisting_data_comp_BR(df, list_of_brain_regions, list_of_groups, list_of_mice):

    headers = df.columns[2:]
    exist_dict = {}

    for BR in list_of_brain_regions:

        corr_compunds = [h[0] for h in headers if h[1] == BR]

        exist_dict[BR] = {}

        for comp in corr_compunds:

            exist_ind = np.where(df[comp, BR] != 0)[0]

            if len(exist_ind) > 0:

                group_spec_dict, missing_groups = find_group_specific_indices(
                    df, exist_ind, list_of_groups)

                exist_dict[BR][comp] = {
                    'viable_indices': exist_ind,
                    'missing_groups': missing_groups,
                    'group_specific_ind': group_spec_dict
                }

            print(' BR : ', BR,
                  'Comnpound : ', comp,
                  ', nb miss mice:', len(list_of_mice) - len(exist_ind),
                  ', nb miss groups', len(missing_groups))

    return exist_dict


def calculate_mean_SD_SEM(df_inc_ratios, exist_dict):
    ''' return the mean, SD, SEM of of the compound according to 
        the viable sample indices provided by <exist_dict>
    '''

    groups = np.unique(df_inc_ratios['group', 'no'])
    n_groups = len(groups)
    measure_list = ['mean', 'SD', 'SEM']

    summary = {measure: pd.DataFrame(groups,
                                     columns=pd.MultiIndex.from_tuples(
                                         [['group', 'no']])
                                     )
               for measure in measure_list
               }

    for BR, BR_dict in exist_dict.items():

        for comp, comp_dict in BR_dict.items():

            summary['mean'][comp, BR] = np.full(n_groups, np.nan)
            summary['SD'][comp, BR] = np.full(n_groups, np.nan)
            summary['SEM'][comp, BR] = np.full(n_groups, np.nan)

            for group, group_ind in comp_dict['group_specific_ind'].items():

                summary['mean'][comp, BR][int(
                    group) - 1] = df_inc_ratios[comp, BR][group_ind].mean()
                summary['SD'][comp, BR][int(
                    group) - 1] = df_inc_ratios[comp, BR][group_ind].std()
                summary['SEM'][comp, BR][int(
                    group) - 1] = df_inc_ratios[comp, BR][group_ind].sem()

    return summary['mean'], summary['SD'], summary['SEM']

 # %% SHIT I AM TOO SCARED TO

 # should be redundant now that ordered correlograms can be used to select subsection


def pearson_correlations_within_BR_custom(treatment_dict, list_of_groups, exist_dict,
                                          df_inc_ratios, name='neurotransmitters', p_value=0.05,
                                          compounds_to_analise=['DA', '5HT', 'NAD']):
    '''input dict and df, 
            produce  df for one correlogram (BR between compounds or compounds between BR), 
                    create correlogram with all significant R vaklues displayed for each compound/ratio 
            --> saved to pdf '''

    pdf = PdfPages(name + '_custom_correlograms_within_BR.pdf')

# SIVA this df will have some NaN values and i dont really understand how it is being delt with with pandas .corr
    # https://stackoverflow.com/questions/57155427/how-does-corr-remove-na-and-null-values

    for treatment in list_of_groups:  # create df for correlations of smetric data : each collumn is a factor to correlate maintining indicies

        # treatment = 1 # to test
        # count = 0

        for BR, BR_dict in exist_dict.items():

            # if count == 1:
            #         break

            # get all indexes of treatment group
            index_of_group = df_inc_ratios.index[df_inc_ratios['group', 'no'] == treatment].tolist(
            )

            df_corr_working = pd.DataFrame('NaN', index=index_of_group, columns=list(
                compounds_to_analise))      # create df to be filled for single BR and treatment

            for comp, comp_dict in BR_dict.items():
                if comp in compounds_to_analise:

                    df_corr_working[comp] = df_inc_ratios[comp,
                                                          BR][comp_dict['group_specific_ind'][treatment]]     # fill df
                else:
                    continue

            # count = count +1

            corr_matrix_BR = df_corr_working.corr(
                method='pearson', min_periods=1)

            fig, ax = plt.subplots(figsize=(16, 10))

            # calculate p_values and set maske to P>0.05
            p_value_matrix = calculate_pvalues(df_corr_working)
            mask = np.invert(np.tril(p_value_matrix < p_value))

            heatmap = sns.heatmap(corr_matrix_BR, vmin=-1, vmax=1,
                                  annot=True, cmap='BrBG', mask=mask, annot_kws={"size": 16})
            # use matplotlib.colorbar.Colorbar object
            cbar = ax.collections[0].colorbar
            # here set the labelsize by 20
            cbar.ax.tick_params(labelsize=16)
            heatmap.set_title(
                BR+' '+str(treatment_dict[treatment]), fontdict={'fontsize': 36}, pad=12, y=1.03)
            ax.set_xticklabels(
                ax.get_xticklabels(),
                rotation=45,
                horizontalalignment='right', fontsize=18)
            ax.set_yticklabels(
                ax.get_yticklabels(), rotation=45, fontsize=18)
            ax.set_ylabel('')
            ax.set_xlabel('')
            pdf.savefig(fig)
            plt.close(fig)

    pdf.close()
    plt.close('all')
    return

 # OLD FUNCTION PRIOR TO RATIO MATCH INDEX DICT


def pearson_correlations_within_compound_ratio(list_of_ratios, list_of_groups, ratio_match_ind_dict, df_inc_ratios, list_of_compounds, treatment_dict,
                                               name='ratio', p_value=0.05, column_order=['OF', 'PL', 'Cg', 'Nac', 'VM', 'DM', 'VL', 'DL', 'SP', 'AM', 'HD', 'HV', 'Th', 'DY',  'VY',  'SN', 'VTA', 'DR',  'MR']):
    '''input dict and df, 
            produce  df for one correlogram (BR between compounds or compounds between BR), 
                    create correlogram with all significant R vaklues displayed for each compound/ratio 
            --> saved to pdf '''

    pdf = PdfPages(name + '_correlograms_between_BRs.pdf')

# SIVA this df will have some NaN values and i dont really understand how it is being delt with with pandas .corr
    # https://stackoverflow.com/questions/57155427/how-does-corr-remove-na-and-null-values

    for treatment in list_of_groups:  # create df for correlations of smetric data : each collumn is a factor to correlate maintining indicies

        for comp in list_of_ratios:

            # get all indexes of treatment group
            index_of_group = df_inc_ratios.index[df_inc_ratios['group', 'no'] == treatment].tolist(
            )

            # df_corr_working = pd.DataFrame('NaN', index=index_of_group, columns=list(ratio_match_ind_dict.keys()))      # create df to be filled for single BR and treatment

            if column_order == False:
                print('no column order')

                df_corr_working = pd.DataFrame(
                    'NaN', index=index_of_group, columns=list(ratio_match_ind_dict.keys()))

                # if you dont have a column order with all of your interest compounds
                if 'HSER' in list(ratio_match_ind_dict.keys()):

                    df_corr_working = df_corr_working.drop(
                        columns='HSER')  # remove HSER

            else:
                print('column ordered:  ', column_order)
                cols = [x for x in column_order if x in list(
                    ratio_match_ind_dict.keys())]  # list of BR in order

                # create df to be filled for single BR and treatment
                df_corr_working = pd.DataFrame(
                    'NaN', index=index_of_group, columns=cols)

            for BR, BR_dict in ratio_match_ind_dict.items():
                compounds_in_this_BR = list(BR_dict.keys())

                if str(comp) in compounds_in_this_BR:
                    df_corr_working[BR] = df_inc_ratios[comp, BR][BR_dict[comp]
                                                                  ['group_specific_ind'][treatment]]     # fill df
                else:
                    print(BR, ' does not have ', comp)

            corr_matrix_BR = df_corr_working.corr(
                method='pearson', min_periods=1)

            fig, ax = plt.subplots(figsize=(16, 10))

            # calculate p_values and set maske to P>0.05
            p_value_matrix = calculate_pvalues(df_corr_working)
            # removing columns that have all 'NaN' vslues to match correlation matrix
            p_value_matrix = p_value_matrix.dropna(axis=1, how='all')
            # removing rows that have all 'NaN' vslues
            p_value_matrix = p_value_matrix.dropna(axis=0, how='all')

            mask = np.invert(np.tril(p_value_matrix < p_value))

            heatmap = sns.heatmap(corr_matrix_BR, vmin=-1, vmax=1,
                                  annot=True, cmap='BrBG', mask=mask, annot_kws={"size": 8})
            # heatmap.set_title(comp+' correlation treatment '+str(treatment), fontdict={'fontsize':18}, pad=12);
            heatmap.set_title(
                comp+'  in  ' + treatment_dict[treatment], fontdict={'fontsize': 18}, pad=12)
            ax.set_xticklabels(
                ax.get_xticklabels(),
                rotation=45,
                horizontalalignment='right')
            ax.set_ylabel('')
            ax.set_xlabel('')
            pdf.savefig(fig)
            plt.close(fig)

    pdf.close()
    plt.close('all')
    return
