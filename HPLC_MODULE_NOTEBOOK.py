#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# %% INSTALL
import os, shutil, itertools, json, time, functools, pickle #GENERIC UTILS

#STAT LIBS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import cm

import numpy as np

from outliers import smirnov_grubbs as grubbs

import pandas as pd

import pingouin as pg

import pprint

import scipy
from scipy.stats import pearsonr
import scipy.stats as stats
from scipy.stats import ttest_ind

import seaborn as sns

from sklearn.preprocessing import StandardScaler # mean = 0 vairance =1
from sklearn.decomposition import PCA

from statannotations.Annotator import Annotator
from statannot import add_stat_annotation

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison

# import weasyprint


######## CONSTANTS ######################
# Constants are like variables that should not be rewritten, they are declared in all caps by convention
ROOT = os.getcwd() #This gives terminal location (terminal working dir)
INPUT_DIR = f'{ROOT}/input'
OUTPUT_DIR = f'{ROOT}/output'
CACHE_DIR = f'{INPUT_DIR}/cache'

COLUMN_ORDER = {'region':['OF', 'PL', 'CC', 'IC', 'M', 'SJ', 'SL1', 'SR1', 'SL6', 'SR6', 'AC', 'V' , 
            'Am', 'dH', 'vH', 'NAc', 'VM', 'DM', 'VL', 'DL', 'MD', 'VPL', 'VPR', 'DG', 
            'Y', 'SC', 'SN', 'VTA', 'DR', 'MR', 'CE' ],
                   'compound': ['DA', 'DOPAC', 'HVA', '3MT', '5HT', '5HIAA', 'GLU', 'GLN']}

CORRELOGRAM_COLUMN_ORDER = { 'compound': COLUMN_ORDER['region'], 'region': COLUMN_ORDER['compound']}


########## UTILITARIES ############
#Check filesystem is set up for write operations
def saveMetadata(filename, treatment_mapping, experimental_info):
    subcache_dir = f"{CACHE_DIR}/{filename.split('.')[0]}"
    checkFilesystem(subcache_dir)
    saveJSON(f"{subcache_dir}/treatment_mapping.json", treatment_mapping)
    print(f"TREATMENT MAPPING {treatment_mapping} SAVED TO {subcache_dir} SUBCACHE")
    saveJSON(f"{subcache_dir}/experimental_info.json", experimental_info)
    print(f"EXPERIMENTAL INFO {experimental_info} SAVED TO {subcache_dir} SUBCACHE")

#This function saves dictionnaries, JSON is a dictionnary text format that you use to not have to reintroduce dictionnaries as variables 
def saveJSON(path, dict_to_save):
    with open(path, 'w', encoding ='utf8') as json_file:
        json.dump(dict_to_save, json_file)

def getMetadata(filename, metadata_type):
    return getJSON(f"{CACHE_DIR}/{filename.split('.')[0]}/{metadata_type}.json")

def getTreatmentMapping(filename):
    return getMetadata(filename, 'treatment_mapping')

def getExperimentalInfo(filename):
    return getMetadata(filename, 'experimental_info')
    
#This function gets JSON files and makes them into python dictionnaries
def getJSON(path):
    with open(path) as outfile:
        loaded_json = json.load(outfile)
    return loaded_json

def checkFileSystem(filepath):
    if not os.path.exists(filepath):
        os.mkdir(filepath)
        
#This checks that the filesystem has all the requisite folders (input, cache, etc..) and creates them if not
def initiateFileSystem():
    checkFileSystem(INPUT_DIR)
    checkFileSystem(CACHE_DIR)
    checkFileSystem(OUTPUT_DIR)

#This function deletes all cached files, it is used when you want to start from square one because all intermediary results will be cached
def resetCache():
    shutil.rmtree(CACHE_DIR)
    os.mkdir(CACHE_DIR)
    print('CACHE CLEARED')

#This function cahces (aka saves in a easily readable format) all dataframes used
def cache(filename, identifier, to_cache):
    filename = filename.split(".")[0]
    cache_subdir = f'{CACHE_DIR}/{filename}'
    checkFilesystem(cache_subdir)
    with open(f'{cache_subdir}/{identifier}.pkl','wb') as file:
        pickle.dump(to_cache, file)
    print(f'CREATED {cache_subdir}/{identifier}.pkl CACHE')

#This function gets the dataframes that are cached
def getCache(filename, identifier):
    filename = filename.split(".")[0]
    print(f'GETTING "{identifier}" FROM "{filename}" CACHE')
    with open(f'{CACHE_DIR}/{filename}/{identifier}.pkl','rb') as file:
        return pickle.load(file)

#This checks if a particulat dataframe/dataset is cached, return boolean
def isCached(filename, identifier):
    filename = filename.split(".")[0]
    return os.path.isfile(f'{CACHE_DIR}/{filename}/{identifier}.pkl')

def saveFigure(fig, identifier, fig_type):
    output_subdir = f"{OUTPUT_DIR}/{fig_type}"
    checkFilesystem(output_subdir)
    fig.savefig(f"{output_subdir}/{identifier}.svg")
    print(f'SAVED {output_subdir}/{identifier}.svg') 
    fig.savefig(f"{output_subdir}/{identifier}.png")
    print(f'SAVED {output_subdir}/{identifier}.png') 

def saveCorrelogram(fig, identifier):
    saveFigure(fig, identifier, 'correlograms')

def applyTreatmentMapping(df, filename):
    filename = filename.split(".")[0]
    treatment_mapping_path = f"{CACHE_DIR}/{filename}/treatment_mapping.json"
    if not os.path.isfile(treatment_mapping_path): #Check treatment mapping is present
        raise Exception("TREATMENT INFORMATION ABSENT, TO SAVE TREATMENT MAPPING RUN 'setTreatment(filename, treatment_mapping)'")
    treatment_mapping = getJSON((treatment_mapping_path))
    new_columns = list(list(treatment_mapping.values())[0].keys()) # Get the future column names from one of the treatments
    df[new_columns] = df.apply(lambda x: treatment_mapping[str(int(x['group_id']))].values(), axis=1, result_type='expand') # Get alll the values and assign to corresponding columns
    return df.explode('experiments').rename(columns={'experiments':'experiment'}) # Duplicate rows belonging to multiple experiment so that groupby can be done later
    
def dictToFilename(dict_to_stringify):
    result = str(dict_to_stringify)
    for replacement in [['{', ''], ['}', ''], [',', '_'], [':', ':'], ["'", ''], ['[', ''], [']', ''], [' ', '']]:
        result = result.replace(*replacement) #This syntaxt wil unpack the list as if I had written 'result.replace(replacement[0], replacement[1])'
    return result
        
########### GETTERS #################    

#Generic df getter
#First checks cache to see if the df already has been built and saved in cache
#If not it uses the builder callback to build the df appropriately
def getOrBuildDf(filename, df_identifier, builder_cb):
    filename_no_extension = filename.split(".")[0]
    if isCached(filename_no_extension, df_identifier): #Check cache to avoid recalcuating from scratch if alreasy done
        return getCache(filename_no_extension, df_identifier)
    print(f'BUILDING "{df_identifier}"')    #Build useing callback otherwise and cache result
    df = builder_cb(filename)
    cache(filename_no_extension, df_identifier, df)
    return df

#The three getters that follow just used the generic function to get the df if cached, injecting their own specific functions to build the df in the case where its not chached
def getRawDf(filename):
    return getOrBuildDf(filename, 'raw_df', buildRawDf)

def getCompoundDf(filename): #TODO: rename to compounds later
    return getOrBuildDf(filename, 'compound_df', buildCompoundDf)
    
def getCompoundAndRatiosDf(filename): #TODO: rename to compounds later
    return getOrBuildDf(filename, 'compound_and_ratios_df', buildCompoundAndRatiosDf)

def getRatiosDf(filename):
    return getOrBuildDf(filename, 'ratios_df', buildRatiosDf)
    
def getRatiosPerRegion(filename, ratios_mapping):
    cache_filename = dictToFilename(ratios_mapping)
    builder_callback_with_param = functools.partial(buildRatiosPerRegionDf, ratios_mapping=ratios_mapping)
    return getOrBuildDf(filename, cache_filename, builder_callback_with_param)

def getCompoundAggregateStatsDf(filename):
    return getAggregateStatsDf(filename, 'compound')

def getRatioAggregateStatsDf(filename):
    return getAggregateStatsDf(filename, 'ratio')

def getAggregateStatsDf(filename, df_type):
    if df_type not in ['ratio', 'compound']: raise Exception("DF TYPE MUST BE IN ['ratio', 'compound']")
    builder_callback_with_param = functools.partial(buildAggregateStatsDf, df_type=df_type)
    return getOrBuildDf(filename, f"{df_type}_aggregate_stats", builder_callback_with_param)


############ BUILDERS #########

#Contains the logic to build the raw df from the csv file
def buildRawDf(filename):
    file_name, file_type = filename.split('.')
    if not file_type == 'csv':
        raise Exception(f'METHOD TO DEAL WITH FILE TYPE "{file_type}" ABSENT')
    if not os.path.isfile(f'{INPUT_DIR}/{filename}'):
        raise Exception(f'FILE {filename} IS ABSENT IN "input/" DIRECTORY')
    return pd.read_csv(f'{INPUT_DIR}/{filename}', header=0).replace(np.nan, 0).rename(columns={'mouse_no':'mouse_id', 'group_no':'group_id'}) # to set all 0 to Nan

#contains the logic to build the df in the new format based on raw df
def buildCompoundDf(filename):
    raw_df = getRawDf(filename)
    new_format_df = raw_df.melt(id_vars=['mouse_id', 'group_id'], value_vars=raw_df.columns[2:])
    new_format_df[['compound', 'region']] = new_format_df.apply(lambda x: x.variable.split('_'), axis=1, result_type='expand')
    new_format_df = applyTreatmentMapping(new_format_df, filename)
    return new_format_df.drop(columns=['variable'])

#Contains the logic to build the ratios df based on the df with the new format
def buildCompoundAndRatiosDf(filename):
    compound_df = getCompoundDf(filename) #.iloc[0:100] #To speed up testing
    ratios_df = pd.merge(left=compound_df, right=compound_df, on=['mouse_id', 'group_id', 'region', 'experiment', 'color', 'treatment'], suffixes=['_1', '_2']) #merge every compound to every other for each mouse
    ratios_df = ratios_df[(ratios_df.compound_1 != ratios_df.compound_2)]
    ratios_df[['compound', 'value']] = ratios_df.apply(lambda x: [f'{x.compound_1}/{x.compound_2}', x.value_1 / x.value_2 if x.value_2 else np.NaN], axis=1, result_type='expand') #calculate ratio
    return pd.concat([compound_df, ratios_df.drop(columns=['compound_1', 'compound_2', 'value_1', 'value_2'])]) #Drop duplicate columns


#Contains the logic to build the ratios df based on the df with the new format
def buildRatiosDf(filename):
    compound_df = getCompoundDf(filename) #.iloc[0:100] #To speed up testing
    merged = pd.merge(left=compound_df, right=compound_df, on='mouse_id', suffixes=['_1', '_2']) #merge every region/compound combo to every other for each mouse
    merged = merged[~((merged.region_1 == merged.region_2) & (merged.compound_1 == merged.compound_2))] #eliminate duplicates (region1 == region2 & C1 == C2)
    # merged = merged[(merged.region_1 < merged.region_2) | (merged.compound_1 < merged.compound_2) #Eliminate duplicates (region1 == region2 & C1 == C2) and cross combinations (row1: region1=A, C1=B, region2=C, C2=D; row2: region1=C, C1=D, region2=A, C2=B))
    #         & ~((merged.region_1 > merged.region_2) & (merged.compound_1 < merged.compound_2))] #Uncomment code to only save half the ration (C1/C2) to divide time it takes by 2
    merged[['compounds', 'ratio']] = merged.apply(lambda x: [f'{x.compound_1}/{x.compound_2}', x.value_1 / x.value_2 if x.value_2 else np.NaN], axis=1, result_type='expand') #calculate the ratio
    return merged.rename(columns={'treatment_1': 'treatment'}).drop(columns=['treatment_2', 'value_1', 'value_2']) #Drop duplicate columns

#Function to get the specific rations (intra region) that you use based on a ratios dictionnary
def buildRatiosPerRegionDf(filename, ratios_mapping):
    ratios_df = getRatiosDf(filename)
    compound_ratios = []
    for compound_1, compound_2_list in ratios_mapping.items():
        for compound_2 in compound_2_list:
            compound_ratios.append(ratios_df[(ratios_df['compound_1']==compound_1) & (ratios_df['compound_2']==compound_2)])
    return pd.concat(compound_ratios)
    
#returns df columns = ['treatment', 'region', 'compound', 'F_value', 'p_value']
def buildAggregateStatsDf(filename, df_type):
    working_df = getCompoundDf(filename) if df_type == 'compound' else getRatiosDf(filename)
    result_ls = []
    for treat_region_comp, groupby_df in working_df.groupby(by =['treatment', 'region', 'compound']):
        F, p, is_valid = [*scipy.stats.shapiro(groupby_df['value']), True] if len(groupby_df) >= 3 else [np.NaN, np.NaN, False]
        mean, std, sem, values = [groupby_df.value.mean(), groupby_df.value.std(), groupby_df.value.sem(), list(groupby_df.value)]  
        result_ls.append([*treat_region_comp, F, p, is_valid, mean, std, sem, values]) #start unpacks the list of strings
    return pd.DataFrame(result_ls, columns= ['treatment', 'region', 'compound', 'shapiro_F', 'shapiro_p', 'is_valid', 'mean', 'std', 'sem', 'values'])

#returns df columns = ['treatment', 'region', 'compound', 'F_value', 'p_value']
def testBuildAggregateStatsDf(filename, df_type):
    working_df = getCompoundDf(filename) if df_type == 'compound' else getRatiosDf(filename)
    descriptive_stats_ls = [] #this one is just describing every region/compound/treament combo
    proper_stats_ls = [] #this one really comparing groups

    for region_comp, region_comp_groupby_df in working_df.groupby(by =['region', 'compound']):
        for treatment, treatment_groupby_df in region_comp_groupby_df.groupby(by =['treatment']):
            F, p, is_valid = [*scipy.stats.shapiro(treatment_groupby_df['value']), True] if len(treatment_groupby_df) >= 3 else [np.NaN, np.NaN, False]
            mean, std, sem, values = [treatment_groupby_df.value.mean(), treatment_groupby_df.value.std(), treatment_groupby_df.value.sem(), list(treatment_groupby_df.value)]  
            descriptive_stats_ls.append([*[treatment, *region_comp], F, p, is_valid, mean, std, sem, values]) #start unpacks the list of strings
            
        f_one, p_one = scipy.stats.f_oneway(g1, g2, g3) #pseudo code
        f_two, p_two = scipy.stats.f_twoway(g1, g2, g3) #pseudo code
        proper_stats_ls.append([treatment, region_comp, F, p]) #pseudo code 
        
    return pd.DataFrame(descriptive_stats_ls, columns= ['treatment', 'region', 'compound', 'shapiro_F', 'shapiro_p', 'is_valid', 'mean', 'std', 'sem', 'values'])


def editOutlier(subselect={'compound': 'DA', 'experiment': 'dose_response', 'region':'CC'}):
    query = '&'.join([f"{column}=='{value}'" for column, value in subselect.items()])
    df.query()

############### CALCULATE/STATISTICS ############

STAT_METHODS = {'pearson': []}


### The following functions are just here to be passed to the pd.corr() method, c and y are the two lists of values (columns) to be correlated
def isSignificant(stat_method_cb, pval_threshold = 0.05): #This is a classic design pattern of which i forget the name again. 
# As you can see here it will return a function. NOT CALL THE FUNCTION, but return it. The point here is to inject variables in the function that is returned.
# when isSignificant(callback, pval_threshold) is called, it declare the anonymous function (lambda) passing the variables, and returns this declaration
# this means that when the lambda function is called later on these will be 'harcoded' in the sense that they are no longer variables to be passed
# Why? because if you look at the usage I pass it to pd.corr() to generate the mask in getPearsonCorrStats(). The thing is that pd.corr() calls the method it is given with (x, y) arguments
# For the mask to be properly generated however, we also want to give a pval threshold to pass, as well as the statistical method that determines significance
# But there is no way to pass this once you are in th pd.corr() context, so this is why you use this design pattern. It is meant to import variable into a context where they are not normally available
    return lambda x, y: stat_method_cb(x, y) >= pval_threshold 

### The 
def getPearson(x,y):
        return scipy.stats.pearsonr(x,y)
        
def getPearsonR(x,y):
        return getPearson(x,y)[0]

def getPearsonPValue(x,y):
        return getPearson(x,y)[1]
#### Up to here ####


#### Here we actually do the proper stats, which include reorganizing the df so that its readily usable
def getPeasonCorrStats(df, pivot_column, p_value_threshold, n_minimum):
    methods = [getPearsonR, isSignificant(getPearsonPValue, p_value_threshold)] #This is the list of methods to pass to df.corr. I know you used to pass 'pearson' as a string but this way the R calculation and pvalue mask are 'linked' by being done on the same line
    pivot_df = df.pivot_table(values='value', index=df['mouse_id'], columns=pivot_column).filter(COLUMN_ORDER[pivot_column]) #Here we just picot our structure to one where each column is a region or compound and each line a mouse
    correlogram_df, p_value_mask = [pivot_df.corr(method=method, min_periods=n_minimum).dropna(axis=0, how='all').dropna(axis=1, how='all') for method in methods] # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
    return correlogram_df, p_value_mask.astype(bool) # Convert 1s and 0s to boolean for the plotting func


def getAndPlotMultipleCorrelograms(filename, selector, p_value_threshold=0.05, n_minimum=5, from_scratch=None): #TODO: Handle columns and more generality
    experiment = selector.pop('experiment', None)
    for correlogram_type, to_correlate_list in selector.items(): #Iterate through the selector dict
        for to_correlate in to_correlate_list:
            getAndPlotSingleCorrelogram(filename, experiment, correlogram_type, to_correlate, p_value_threshold, n_minimum, from_scratch)

def askMultipleChoice(question, choices):
    return choices[int(input(question + '\n' + '\n'.join([f'{i}: {choice}' for i, choice in choices.items()]) + '\n'))]

def getAndPlotSingleCorrelogram(filename, experiment=None, correlogram_type=None, to_correlate=None, p_value_threshold=0.05, n_minimum=5, columns=None, from_scratch=None):
    experiments = getExperimentalInfo(filename)
    experiment = experiment if experiment else askMultipleChoice("Which experiment?", {i: experiment for i, experiment in enumerate(experiments)})
    correlogram_type = correlogram_type if correlogram_type else askMultipleChoice("Which correlogram?", {0: 'compound', 1: 'region'})
    to_correlate = to_correlate if to_correlate else input(f"""Which {correlogram_type}?
                    (Enter simple {correlogram_type} or {correlogram_type}/{correlogram_type} for ratio,
                    Use {correlogram_type} or {correlogram_type}/{correlogram_type} for simple correlogram, or {correlogram_type}-{correlogram_type} or {correlogram_type}/{correlogram_type}-{correlogram_type} for square correlogram
                    Possibilities: {COLUMN_ORDER[correlogram_type]}
                    """)
    columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
    identifier = f"{experiment}_{correlogram_type}_{to_correlate.replace('/', ':')}_{(',').join(columns)}"
    from_scratch = from_scratch if from_scratch is not None else input("Recalculate figure even if previous version exists? (y/n)") == 'y'
    if from_scratch or not isCached(filename, identifier):
        fig = buildSingleCorrelogram(filename, experiment, correlogram_type, to_correlate.split('-'), p_value_threshold, n_minimum, columns)
        cache(filename, identifier, fig)
        saveCorrelogram(fig, identifier)
    else : fig = getCache(filename, identifier)
    fig.show()


def buildSingleCorrelogram(filename, experiment, correlogram_type, to_correlate, p_value_threshold, n_minimum, columns):
    compound_and_ratios_df = getCompoundAndRatiosDf(filename) #this is not the full ratios df, its only intra region compound ratios for nom
    value_type = {'compound': 'region', 'region': 'compound'}[correlogram_type]
    subselection_df = compound_and_ratios_df[(compound_and_ratios_df.experiment == experiment) & (compound_and_ratios_df[value_type].isin(columns))].query('|'.join([f"{correlogram_type}=='{value}'" for value in to_correlate]))
    pivot_columns = {'region':['region', 'compound'], 'compound': ['compound', 'region']}[correlogram_type]
    pivot_column_value = to_correlate if len(to_correlate) == 2 else to_correlate + to_correlate #Because table is îvoted on region and compound even in the case of simple correlogram I have to duplicate the selector in that case to avoid ugly labelling
    correlograms = []
    for treatment, group_df in subselection_df.groupby(by=['treatment']): 
        pivot_df = group_df.pivot_table(values='value', index=group_df['mouse_id'], columns=pivot_columns) #Here we just picot our structure to one where each column is a region or compound and each line a mouse    correlogram_df, p_value_mask = [pivot_df.corr(method=method, min_periods=n_minimum).dropna(axis=0, how='all').dropna(axis=1, how='all') for method in methods] # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
        methods = [getPearsonR, isSignificant(getPearsonPValue, p_value_threshold)] #This is the list of methods to pass to df.corr. I know you used to pass 'pearson' as a string but this way the R calculation and pvalue mask are 'linked' by being done on the same line
        correlogram_df, p_value_mask = [pivot_df.corr(method=method, min_periods=n_minimum).loc[tuple(pivot_column_value)].dropna(axis=0, how='all').dropna(axis=1, how='all') for method in methods] # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
        correlograms.append([correlogram_df, p_value_mask.astype(bool), treatment[0], to_correlate])
    fig = plotCorrelograms(correlograms)
    return fig


def plotCorrelograms(correlograms):
    fig, axs = plt.subplots(2, 2, figsize=(10,5))
    axs = list(itertools.chain.from_iterable(axs)) # Put all the axes into a list at the same level
    for (correlogram_df, p_value_mask, treatment, subvalues), ax in zip(correlograms, axs):
        plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax)
    return fig
    
def plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax):
    heatmap = sns.heatmap(correlogram_df, vmin=-1, vmax=1, annot=True, cmap='BrBG', mask=p_value_mask, annot_kws={"size": 8}, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    heatmap.set_title(f"{'-'.join(subvalues)} in {treatment}", fontdict={'fontsize': 10}, pad=20)
    if len(subvalues) == 1: 
        ax.set_ylabel('')
        ax.set_xlabel('')
    elif len(subvalues) == 2:
        ax.set_ylabel(subvalues[0])
        ax.set_xlabel(subvalues[1])

def plotAnything(anything):
    return plt(anything)

def showOutliers(raw_df):
    return plotAnything(getOutliers(raw_df))

def getOutliers(raw_df):
    return [doRawDfGrubbs(raw_df), doRawDfGrubbs(raw_df)]
    

def doRawDfGrubbs(raw_df):
    result_list = []
    for group in raw_df.groupby('treatment'):
        result_list.append(grubbsTest(raw_df)) #this func will loop through the df to feel grubbsTest
    return result_list
    
    

def grubbsTest(group_list): #include the vairable type in name i.e. group_list series
    return


######## INIT ##########
#Start by checking filesystem has all the folders necessary for read/write operations (cache) or create them otherwise
initiateFileSystem()














