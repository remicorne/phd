#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# %% INSTALL
import os
import shutil
import itertools
import json
import time
import functools
import pickle  # GENERIC UTILS

# STAT LIBS
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

from sklearn.preprocessing import StandardScaler  # mean = 0 vairance =1
from sklearn.decomposition import PCA

from statannotations.Annotator import Annotator
from statannot import add_stat_annotation

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison

# import weasyprint


######## CONSTANTS ######################
# Constants are like variables that should not be rewritten, they are declared in all caps by convention
ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
INPUT_DIR = f"{ROOT}/input"
OUTPUT_DIR = f"{ROOT}/output"
CACHE_DIR = f"{INPUT_DIR}/cache"

COLUMN_ORDER = {
    "region": [
        "OF",
        "PL",
        "CC",
        "IC",
        "M",
        "SJ",
        "SL1",
        "SL6",
        "SR6",
        "SR1",
        "AC",
        "V",
        "Am",
        "dH",
        "vH",
        "NAc",
        "VM",
        "DM",
        "VL",
        "DL",
        "MD",
        "VPL",
        "VPR",
        "DG",
        "Y",
        "SC",
        "SN",
        "VTA",
        "DR",
        "MR",
        "CE",
    ],
    "compound": ["DA", "DOPAC", "HVA", "3MT", "5HT", "5HIAA", "GLU", "GLN"],
}


CORRELOGRAM_COLUMN_ORDER = {
    "compound": COLUMN_ORDER["region"],
    "region": COLUMN_ORDER["compound"],
}


########## UTILITARIES ############
# Check filesystem is set up for write operations
def saveMetadata(filename, treatment_mapping, experimental_info):
    subcache_dir = f"{CACHE_DIR}/{filename.split('.')[0]}"
    checkFileSystem(subcache_dir)
    saveJSON(f"{subcache_dir}/treatment_mapping.json", treatment_mapping)
    print(f"TREATMENT MAPPING {treatment_mapping} SAVED TO {subcache_dir} SUBCACHE")
    saveJSON(f"{subcache_dir}/experimental_info.json", experimental_info)
    print(f"EXPERIMENTAL INFO {experimental_info} SAVED TO {subcache_dir} SUBCACHE")


# This function saves dictionnaries, JSON is a dictionnary text format that you use to not have to reintroduce dictionnaries as variables


def saveJSON(path, dict_to_save):
    with open(path, "w", encoding="utf8") as json_file:
        json.dump(dict_to_save, json_file)


def getMetadata(filename, metadata_type):
    return getJSON(f"{CACHE_DIR}/{filename.split('.')[0]}/{metadata_type}.json")


def getTreatmentMapping(filename):
    return getMetadata(filename, "treatment_mapping")


def getExperimentalInfo(filename):
    return getMetadata(filename, "experimental_info")


# This function gets JSON files and makes them into python dictionnaries


def getJSON(path):
    with open(path) as outfile:
        loaded_json = json.load(outfile)
    return loaded_json


def checkFileSystem(filepath):
    if not os.path.exists(filepath):
        os.mkdir(filepath)


# This checks that the filesystem has all the requisite folders (input, cache, etc..) and creates them if not


def initiateFileSystem():
    checkFileSystem(INPUT_DIR)
    checkFileSystem(CACHE_DIR)
    checkFileSystem(OUTPUT_DIR)


# This function deletes all cached files, it is used when you want to start from square one because all intermediary results will be cached


def resetCache():
    shutil.rmtree(CACHE_DIR)
    os.mkdir(CACHE_DIR)
    print("CACHE CLEARED")


# This function cahces (aka saves in a easily readable format) all dataframes used


def cache(filename, identifier, to_cache):
    filename = filename.split(".")[0]
    cache_subdir = f"{CACHE_DIR}/{filename}"
    checkFileSystem(cache_subdir)
    with open(f"{cache_subdir}/{identifier}.pkl", "wb") as file:
        pickle.dump(to_cache, file)
    print(f"CREATED {cache_subdir}/{identifier}.pkl CACHE")


# This function gets the dataframes that are cached


def getCache(filename, identifier):
    filename = filename.split(".")[0]
    print(f'GETTING "{identifier}" FROM "{filename}" CACHE')
    with open(f"{CACHE_DIR}/{filename}/{identifier}.pkl", "rb") as file:
        return pickle.load(file)


# This checks if a particulat dataframe/dataset is cached, return boolean


def isCached(filename, identifier):
    filename = filename.split(".")[0]
    return os.path.isfile(f"{CACHE_DIR}/{filename}/{identifier}.pkl")


def saveFigure(fig, identifier, fig_type):
    output_subdir = f"{OUTPUT_DIR}/{fig_type}"
    checkFileSystem(output_subdir)
    fig.savefig(f"{output_subdir}/{identifier}.svg")  # dpi also?
    print(f"SAVED {output_subdir}/{identifier}.svg")
    # https://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
    fig.savefig(f"{output_subdir}/{identifier}.png", dpi=fig.dpi)
    print(f"SAVED {output_subdir}/{identifier}.png")


def saveCorrelogram(fig, identifier):
    saveFigure(fig, identifier, "correlograms")


def saveHistogram(fig, identifier):
    saveFigure(fig, identifier, "histograms")


def saveHTHistogram(fig, identifier):
    saveFigure(fig, identifier, "HT_histograms")


def applyTreatmentMapping(df, filename):
    filename = filename.split(".")[0]
    treatment_mapping_path = f"{CACHE_DIR}/{filename}/treatment_mapping.json"
    # Check treatment mapping is present
    if not os.path.isfile(treatment_mapping_path):
        raise Exception(
            "TREATMENT INFORMATION ABSENT, TO SAVE TREATMENT MAPPING RUN 'setTreatment(filename, treatment_mapping)'"
        )
    treatment_mapping = getJSON((treatment_mapping_path))
    # Get the future column names from one of the treatments
    new_columns = list(list(treatment_mapping.values())[0].keys())
    df[new_columns] = df.apply(
        lambda x: treatment_mapping[str(int(x["group_id"]))].values(),
        axis=1,
        result_type="expand",
    )  # Get alll the values and assign to corresponding columns
    # Duplicate rows belonging to multiple experiment so that groupby can be done later
    return df.explode("experiments").rename(columns={"experiments": "experiment"})


def dictToFilename(dict_to_stringify):
    result = str(dict_to_stringify)
    for replacement in [
        ["{", ""],
        ["}", ""],
        [",", "_"],
        [":", ":"],
        ["'", ""],
        ["[", ""],
        ["]", ""],
        [" ", ""],
    ]:
        # This syntaxt wil unpack the list as if I had written 'result.replace(replacement[0], replacement[1])'
        result = result.replace(*replacement)
    return result


########### GETTERS #################

# Generic df getter
# First checks cache to see if the df already has been built and saved in cache
# If not it uses the builder callback to build the df appropriately


def getOrBuildDf(filename, df_identifier, builder_cb):
    filename_no_extension = filename.split(".")[0]
    # Check cache to avoid recalcuating from scratch if alreasy done
    if isCached(filename_no_extension, df_identifier):
        return getCache(filename_no_extension, df_identifier)
    # Build useing callback otherwise and cache result
    print(f'BUILDING "{df_identifier}"')
    df = builder_cb(filename)
    cache(filename_no_extension, df_identifier, df)
    return df


# The three getters that follow just used the generic function to get the df if cached, injecting their own specific functions to build the df in the case where its not chached


def getRawDf(filename):
    return getOrBuildDf(filename, "raw_df", buildRawDf)


def getCompoundDf(filename):  # TODO: rename to compounds later
    return getOrBuildDf(filename, "compound_df", buildCompoundDf)


def getCompoundAndRatiosDf(filename):  # TODO: rename to compounds later
    return getOrBuildDf(filename, "compound_and_ratios_df", buildCompoundAndRatiosDf)


def getRatiosDf(filename):
    return getOrBuildDf(filename, "ratios_df", buildRatiosDf)


def getRatiosPerRegion(filename, ratios_mapping):
    cache_filename = dictToFilename(ratios_mapping)
    builder_callback_with_param = functools.partial(
        buildRatiosPerRegionDf, ratios_mapping=ratios_mapping
    )
    return getOrBuildDf(filename, cache_filename, builder_callback_with_param)


def getCompoundAggregateStatsDf(filename):
    return getAggregateStatsDf(filename, "compound")


def getRatioAggregateStatsDf(filename):
    return getAggregateStatsDf(filename, "ratio")


def getAggregateStatsDf(filename, df_type):
    if df_type not in ["ratio", "compound"]:
        raise Exception("DF TYPE MUST BE IN ['ratio', 'compound']")
    builder_callback_with_param = functools.partial(
        buildAggregateStatsDf, df_type=df_type
    )
    return getOrBuildDf(
        filename, f"{df_type}_aggregate_stats", builder_callback_with_param
    )


############ BUILDERS #########


# Contains the logic to build the raw df from the csv file
def buildRawDf(filename):
    file_name, file_type = filename.split(".")
    if not file_type == "csv":
        raise Exception(f'METHOD TO DEAL WITH FILE TYPE "{file_type}" ABSENT')
    if not os.path.isfile(f"{INPUT_DIR}/{filename}"):
        raise Exception(f'FILE {filename} IS ABSENT IN "input/" DIRECTORY')
    # to set all 0 to Nan
    return pd.read_csv(f"{INPUT_DIR}/{filename}", header=0).replace(np.nan, 0)


# contains the logic to build the df in the new format based on raw df


def buildCompoundDf(filename):
    raw_df = getRawDf(filename)
    new_format_df = raw_df.melt(
        id_vars=["mouse_id", "group_id"], value_vars=raw_df.columns[2:]
    )
    new_format_df[["compound", "region"]] = new_format_df.apply(
        lambda x: x.variable.split("_"), axis=1, result_type="expand"
    )
    new_format_df = applyTreatmentMapping(new_format_df, filename)
    return new_format_df.drop(columns=["variable"])


# Contains the logic to build the ratios df based on the df with the new format


def buildCompoundAndRatiosDf(filename):
    compound_df = getCompoundDf(filename)  # .iloc[0:100] #To speed up testing
    ratios_df = pd.merge(
        left=compound_df,
        right=compound_df,
        on=["mouse_id", "group_id", "region", "experiment", "color", "treatment"],
        suffixes=["_1", "_2"],
    )  # merge every compound to every other for each mouse
    ratios_df = ratios_df[(ratios_df.compound_1 != ratios_df.compound_2)]
    ratios_df[["compound", "value"]] = ratios_df.apply(
        lambda x: [
            f"{x.compound_1}/{x.compound_2}",
            x.value_1 / x.value_2 if x.value_2 else np.NaN,
        ],
        axis=1,
        result_type="expand",
    )  # calculate ratio
    # Drop duplicate columns
    return pd.concat(
        [
            compound_df,
            ratios_df.drop(columns=["compound_1", "compound_2", "value_1", "value_2"]),
        ]
    )


# Contains the logic to build the ratios df based on the df with the new format
def buildRatiosDf(filename):
    compound_df = getCompoundDf(filename)  # .iloc[0:100] #To speed up testing
    merged = pd.merge(
        left=compound_df, right=compound_df, on="mouse_id", suffixes=["_1", "_2"]
    )  # merge every region/compound combo to every other for each mouse
    # eliminate duplicates (region1 == region2 & C1 == C2)
    merged = merged[
        ~(
            (merged.region_1 == merged.region_2)
            & (merged.compound_1 == merged.compound_2)
        )
    ]
    # merged = merged[(merged.region_1 < merged.region_2) | (merged.compound_1 < merged.compound_2) #Eliminate duplicates (region1 == region2 & C1 == C2) and cross combinations (row1: region1=A, C1=B, region2=C, C2=D; row2: region1=C, C1=D, region2=A, C2=B))
    #         & ~((merged.region_1 > merged.region_2) & (merged.compound_1 < merged.compound_2))] #Uncomment code to only save half the ration (C1/C2) to divide time it takes by 2
    merged[["compounds", "ratio"]] = merged.apply(
        lambda x: [
            f"{x.compound_1}/{x.compound_2}",
            x.value_1 / x.value_2 if x.value_2 else np.NaN,
        ],
        axis=1,
        result_type="expand",
    )  # calculate the ratio
    # Drop duplicate columns
    return merged.rename(columns={"treatment_1": "treatment"}).drop(
        columns=["treatment_2", "value_1", "value_2"]
    )


# Function to get the specific rations (intra region) that you use based on a ratios dictionnary


def buildRatiosPerRegionDf(filename, ratios_mapping):
    ratios_df = getRatiosDf(filename)
    compound_ratios = []
    for compound_1, compound_2_list in ratios_mapping.items():
        for compound_2 in compound_2_list:
            compound_ratios.append(
                ratios_df[
                    (ratios_df["compound_1"] == compound_1)
                    & (ratios_df["compound_2"] == compound_2)
                ]
            )
    return pd.concat(compound_ratios)


# returns df columns = ['treatment', 'region', 'compound', 'F_value', 'p_value']


def buildAggregateStatsDf(filename, df_type):
    working_df = (
        getCompoundDf(filename) if df_type == "compound" else getRatiosDf(filename)
    )
    result_ls = []
    for treat_region_comp, groupby_df in working_df.groupby(
        by=["treatment", "region", "compound", "experiment"]
    ):
        F, p, is_valid = (
            [*scipy.stats.shapiro(groupby_df["value"]), True]
            if len(groupby_df) >= 3
            else [np.NaN, np.NaN, False]
        )
        mean, std, sem, values = [
            groupby_df.value.mean(),
            groupby_df.value.std(),
            groupby_df.value.sem(),
            list(groupby_df.value),
        ]
        # start unpacks the list of strings
        result_ls.append([*treat_region_comp, F, p, is_valid, mean, std, sem, values])
    return pd.DataFrame(
        result_ls,
        columns=[
            "treatment",
            "region",
            "compound",
            "experiment",
            "shapiro_F",
            "shapiro_p",
            "is_valid",
            "mean",
            "std",
            "sem",
            "values",
        ],
    )


# returns df columns = ['treatment', 'region', 'compound', 'F_value', 'p_value']


def testBuildAggregateStatsDf(filename, df_type):
    working_df = (
        getCompoundDf(filename) if df_type == "compound" else getRatiosDf(filename)
    )
    # this one is just describing every region/compound/treament combo
    descriptive_stats_ls = []
    proper_stats_ls = []  # this one really comparing groups

    for region_comp, region_comp_groupby_df in working_df.groupby(
        by=["region", "compound"]
    ):
        for treatment, treatment_groupby_df in region_comp_groupby_df.groupby(
            by=["treatment"]
        ):
            F, p, is_valid = (
                [*scipy.stats.shapiro(treatment_groupby_df["value"]), True]
                if len(treatment_groupby_df) >= 3
                else [np.NaN, np.NaN, False]
            )
            mean, std, sem, values = [
                treatment_groupby_df.value.mean(),
                treatment_groupby_df.value.std(),
                treatment_groupby_df.value.sem(),
                list(treatment_groupby_df.value),
            ]
            # start unpacks the list of strings
            descriptive_stats_ls.append(
                [*[treatment, *region_comp], F, p, is_valid, mean, std, sem, values]
            )

        f_one, p_one = scipy.stats.f_oneway(g1, g2, g3)  # pseudo code
        f_two, p_two = scipy.stats.f_twoway(g1, g2, g3)  # pseudo code
        proper_stats_ls.append([treatment, region_comp, F, p])  # pseudo code

    return pd.DataFrame(
        descriptive_stats_ls,
        columns=[
            "treatment",
            "region",
            "compound",
            "shapiro_F",
            "shapiro_p",
            "is_valid",
            "mean",
            "std",
            "sem",
            "values",
        ],
    )


def editOutlier(
    subselect={"compound": "DA", "experiment": "dose_response", "region": "CC"}
):
    query = "&".join([f"{column}=='{value}'" for column, value in subselect.items()])
    df.query()


############### CALCULATE/STATISTICS ############
import warnings


# This is a decorator design pattern. Its basically a function that wraps another function and does some operations before
# Calling the wrapped function and returning the result. In this cas it extracts the names of the parameter for the stat function called
# This makes is possible to generically call stat functions and pass them every variable we have all the while not knowing exactly
# Which one(s) they need. If we don't do this, we will be rejected when passing unexpected arguments to the function
# The decorator will look at the function definition to get only the relevant ones from the kwargs variable (keyword args)
# We could also have used **kwargs in all the stat functions and done the parameter selection after that but the decorator saves codelines and also its cool
def select_params(stat_function):
    def wrapper(*args, **kwargs):
        warnings.filterwarnings(
            "error"
        )  # This is to have warning raise exceptions as if they were errors in order to capture and store them
        # This here selects the name of the arguments the function need (it uses internal python variable, rtm if necessary)
        stat_function_args = stat_function.__code__.co_varnames[
            : stat_function.__code__.co_argcount
        ]
        # Select the argument based on their name from kwargs, kwargs is all the keyword arguments passed to the function including the useless ones
        selected_args = [kwargs[arg_name] for arg_name in stat_function_args]
        try:
            # call the function with the args it needs
            return ["OK", stat_function(*selected_args)]
        except Warning as w:
            warnings.resetwarnings()
            return ["WARNING", stat_function(*selected_args)]
        except Exception as e:
            return ["ERROR", str(e)]

    return wrapper


# The following functions are just here to be passed to the pd.corr() method, c and y are the two lists of values (columns) to be correlated
# This is a classic design pattern of which i forget the name again.
def isSignificant(stat_method_cb, pval_threshold=0.05):
    # As you can see here it will return a function. NOT CALL THE FUNCTION, but return it. The point here is to inject variables in the function that is returned.
    # when isSignificant(callback, pval_threshold) is called, it declare the anonymous function (lambda) passing the variables, and returns this declaration
    # this means that when the lambda function is called later on these will be 'harcoded' in the sense that they are no longer variables to be passed
    # Why? because if you look at the usage I pass it to pd.corr() to generate the mask in getPearsonCorrStats(). The thing is that pd.corr() calls the method it is given with (x, y) arguments
    # For the mask to be properly generated however, we also want to give a pval threshold to pass, as well as the statistical method that determines significance
    # But there is no way to pass this once you are in th pd.corr() context, so this is why you use this design pattern. It is meant to import variable into a context where they are not normally available
    return lambda x, y: stat_method_cb(x, y) >= pval_threshold


# The


def getPearson(x, y):
    return scipy.stats.pearsonr(x, y)


def getPearsonR(x, y):
    return getPearson(x, y).statistic


def getPearsonPValue(x, y):
    return getPearson(x, y).pvalue


@select_params
def getTukey(data, p_value_threshold):
    columns, *stats_data = pairwise_tukeyhsd(
        endog=data["value"], groups=data["treatment"], alpha=p_value_threshold
    )._results_table.data
    return pd.DataFrame(stats_data, columns=columns)


@select_params
def getOneWayAnova(data):
    F_value, p_value = scipy.stats.f_oneway(
        *[list(group_df["value"]) for treatment, group_df in data.groupby("treatment")]
    )
    # print(f'oneWAY_ANOVA F_value: {F_value}, p_value: {p_value}')
    return pd.DataFrame([[F_value, p_value]], columns=["F", "p_value"])


@select_params
def getTwoWayAnova(data, independant_vars):
    data[independant_vars] = data.apply(
        lambda x: [var in x["treatment"] for var in independant_vars],
        axis=1,
        result_type="expand",
    )
    return pg.anova(
        data=data,
        dv="value",
        between=independant_vars,
        detailed=True,
    ).round(3)


def getQuantitativeStats(filename, p_value_threshold, from_scratch=False):
    identifier = "quantitative_stats"
    if from_scratch or not isCached(filename, identifier):
        quantitative_stats_df = buildQuantitativeStatsDf(filename, p_value_threshold)
        cache(filename, identifier, quantitative_stats_df)
    else:
        quantitative_stats_df = getCache(filename, identifier)
    return quantitative_stats_df


def buildQuantitativeStatsDf(filename, p_value_threshold):
    experimental_info = getExperimentalInfo(filename)
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    data_keys = ["experiment", "compound", "region"]
    quantitative_stats = []
    for experiment_compound_region, data in compound_and_ratios_df.groupby(
        by=data_keys
    ):
        experiment, compound, region = experiment_compound_region
        if (
            experiment == "dose_response"
        ):  ## TODO: Remove when info mapping is completed
            continue
        experiment_info = experimental_info[experiment]
        current_result = [
            experiment,
            compound,
            region,
        ]
        for test in experiment_info["quantitative_statistics"]:
            current_result.extend(
                QUANTITATIVE_STAT_METHODS[test](
                    data=data,
                    independant_vars=experiment_info["independant_vars"],
                    p_value_threshold=p_value_threshold,
                )
            )
        quantitative_stats.append(current_result)
    return pd.DataFrame(
        quantitative_stats,
        columns=[
            data_keys
            + flatten(
                [
                    [f"{test}_status", f"{test}_result"]
                    for test in experiment_info["quantitative_statistics"]
                ]
            )
        ],
    )


QUANTITATIVE_STAT_METHODS = {
    "twoway_anova": getTwoWayAnova,
    "oneway_anova": getOneWayAnova,
    "tukey": getTukey,
}


QUALITATIVE_STAT_METHODS = {"pearson": getPearson}


# STAT_METHODS = QUANTITATIVE_STAT_METHODS + QUALITATIVE_STAT_METHODS
#### Up to here ####


# #### Here we actually do the proper stats, which include reorganizing the df so that its readily usable
# def getPeasonCorrStats(df, pivot_column, p_value_threshold, n_minimum):
#     methods = [getPearsonR, isSignificant(getPearsonPValue, p_value_threshold)] #This is the list of methods to pass to df.corr. I know you used to pass 'pearson' as a string but this way the R calculation and pvalue mask are 'linked' by being done on the same line
#     pivot_df = df.pivot_table(values='value', index=df['mouse_id'], columns=pivot_column).filter(COLUMN_ORDER[pivot_column]) #Here we just picot our structure to one where each column is a region or compound and each line a mouse
#     correlogram_df, p_value_mask = [pivot_df.corr(method=method, min_periods=n_minimum).dropna(axis=0, how='all').dropna(axis=1, how='all') for method in methods] # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
#     return correlogram_df, p_value_mask.astype(bool) # Convert 1s and 0s to boolean for the plotting func


# TODO: Handle columns and more generality
def getAndPlotMultipleCorrelograms(
    filename, selector, p_value_threshold=0.05, n_minimum=5, from_scratch=None
):
    experiment = selector.pop("experiment", None)
    for (
        correlogram_type,
        to_correlate_list,
    ) in selector.items():  # Iterate through the selector dict
        for to_correlate in to_correlate_list:
            getAndPlotSingleCorrelogram(
                filename,
                experiment,
                correlogram_type,
                to_correlate,
                p_value_threshold,
                n_minimum,
                from_scratch,
            )


def askMultipleChoice(question, choices):
    return choices[
        int(
            input(
                question
                + "\n"
                + "\n".join([f"{i}: {choice}" for i, choice in choices.items()])
                + "\n"
            )
        )
    ]


def getAndPlotSingleCorrelogram(
    filename,
    experiment=None,
    correlogram_type=None,
    to_correlate=None,
    p_value_threshold=0.05,
    n_minimum=5,
    columns=None,
    from_scratch=None,
):
    experiments = getExperimentalInfo(filename)
    experiment = (
        experiment
        if experiment
        else askMultipleChoice(
            "Which experiment?",
            {i: experiment for i, experiment in enumerate(experiments)},
        )
    )
    correlogram_type = (
        correlogram_type
        if correlogram_type
        else askMultipleChoice("Which correlogram?", {0: "compound", 1: "region"})
    )
    to_correlate = (
        to_correlate
        if to_correlate
        else input(
            f"""Which {correlogram_type}?
                    (Enter simple {correlogram_type} or {correlogram_type}/{correlogram_type} for ratio,
                    Use {correlogram_type} or {correlogram_type}/{correlogram_type} for simple correlogram, or {correlogram_type}-{correlogram_type} or {correlogram_type}/{correlogram_type}-{correlogram_type} for square correlogram
                    Possibilities: {COLUMN_ORDER[correlogram_type]}
                    """
        )
    )
    columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
    identifier = f"{experiment}_{correlogram_type}_{to_correlate.replace('/', ':')}_{(',').join(columns)}"
    from_scratch = (
        from_scratch
        if from_scratch is not None
        else input("Recalculate figure even if previous version exists? (y/n)") == "y"
    )
    if from_scratch or not isCached(filename, identifier):
        fig = buildSingleCorrelogram(
            filename,
            experiment,
            correlogram_type,
            to_correlate.split("-"),
            p_value_threshold,
            n_minimum,
            columns,
        )
        cache(filename, identifier, fig)
        saveCorrelogram(fig, identifier)
    else:
        fig = getCache(filename, identifier)
    fig.show()


def buildSingleCorrelogram(
    filename,
    experiment,
    correlogram_type,
    to_correlate,
    p_value_threshold,
    n_minimum,
    columns,
):
    # this is not the full ratios df, its only intra region compound ratios for nom
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    value_type = {"compound": "region", "region": "compound"}[correlogram_type]
    subselection_df = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df[value_type].isin(columns))
    ].query("|".join([f"{correlogram_type}=='{value}'" for value in to_correlate]))
    columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
    pivot_columns = {
        "region": ["region", "compound"],
        "compound": ["compound", "region"],
    }[correlogram_type]
    # Because table is îvoted on region and compound even in the case of simple correlogram I have to duplicate the selector in that case to avoid ugly labelling
    pivot_column_value = (
        to_correlate if len(to_correlate) == 2 else to_correlate + to_correlate
    )
    correlograms = []
    treatment_mapping = getTreatmentMapping(filename)
    order = [
        treatment_mapping[str(group)]["treatment"]
        for group in getExperimentalInfo(filename)[experiment]["groups"]
    ]
    subselection_df_ordered = subselection_df.iloc[
        subselection_df["treatment"]
        .astype(pd.CategoricalDtype(categories=order))
        .argsort()
    ]  # order df in order to plot, then sort=False in .groupby
    for treatment, group_df in subselection_df_ordered.groupby(
        by=["treatment"], sort=False
    ):
        # Here we just picot our structure to one where each column is a region or compound and each line a mouse    correlogram_df, p_value_mask = [pivot_df.corr(method=method, min_periods=n_minimum).dropna(axis=0, how='all').dropna(axis=1, how='all') for method in methods] # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
        pivot_df = group_df.pivot_table(
            values="value", index=group_df["mouse_id"], columns=pivot_columns
        )
        # This is the list of methods to pass to df.corr. I know you used to pass 'pearson' as a string but this way the R calculation and pvalue mask are 'linked' by being done on the same line
        methods = [getPearsonR, isSignificant(getPearsonPValue, p_value_threshold)]
        # order columns in desired plotting order
        columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
        pivot_columns_ordered = sorted(
            pivot_df.columns,
            key=lambda x: columns.index(x[1]) if x[1] in columns else float("inf"),
        )
        pivot_df_ordered = pivot_df[pivot_columns_ordered]
        correlogram_df, p_value_mask = [
            pivot_df_ordered.corr(method=method, min_periods=n_minimum)
            .loc[tuple(pivot_column_value)]
            .dropna(axis=0, how="all")
            .dropna(axis=1, how="all")
            for method in methods
        ]  # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
        correlograms.append(
            [correlogram_df, p_value_mask.astype(bool), treatment[0], to_correlate]
        )
    fig = plotCorrelograms(correlograms)
    return fig


def getSingleHistogram(
    filename, experiment, compound, region, p_value_threshold, from_scratch=False
):
    identifier = f"{experiment}_for_{compound}_in_{region}"
    from_scratch = (
        from_scratch
        if from_scratch is not None
        else input("Recalculate figure even if previous version exists? (y/n)") == "y"
    )
    if from_scratch or not isCached(filename, identifier):
        fig = buildSingleHistogram(
            filename, experiment, compound, region, p_value_threshold
        )
        cache(filename, identifier, fig)
        saveHistogram(fig, identifier)
    else:
        fig = getCache(filename, identifier)
    fig.show()


def buildSingleHistogram(filename, experiment, compound, region, p_value_threshold):
    # this is not the full ratios df, its only intra region compound ratios for nom
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    subselection_df = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df.compound == compound)
        & (compound_and_ratios_df.region == region)
    ]
    subselection_df = subselection_df[["value", "mouse_id", "treatment"]]
    treatment_mapping = getTreatmentMapping(filename)
    experimental_info = getExperimentalInfo(filename)[experiment]
    palette = {
        info["treatment"]: info["color"] for number, info in treatment_mapping.items()
    }
    order = [
        treatment_mapping[str(group)]["treatment"]
        for group in experimental_info["groups"]
    ]
    # REMI: i commented this as its missing a : but idk where - i just need to work on plotters for correlograms
    # QUANTITATIVE_STAT_METHODS[stat_name](subselection_df, experimental_info) for stat_name, necessary_for_diplay in experimental_info['quantitative_statistics'].items()}
    fig, ax = plt.subplots(figsize=(20, 10))
    ax = sns.barplot(
        x="treatment",
        y="value",
        data=subselection_df,
        palette=palette,
        ci=68,
        order=order,
        capsize=0.1,
        alpha=0.8,
        errcolor=".2",
        edgecolor=".2",
    )
    ax = sns.swarmplot(
        x="treatment",
        y="value",
        palette=palette,
        order=order,
        data=subselection_df,
        edgecolor="k",
        linewidth=1,
        linestyle="-",
    )
    ax.tick_params(labelsize=24)
    ax.set_ylabel("ng/mg of tissue", fontsize=24)
    if "/" in compound:
        ax.set_ylabel(" ", fontsize=24)
    ax.set_xlabel(" ", fontsize=20)  # treatments
    ax.set_title(compound + " in " + region, y=1.04, fontsize=34)  # '+/- 68%CI'
    sns.despine(left=False)
    return fig


def put_significnce_stars(
    stat_data,
    ax,
    treatment_dict,
    test_path,
    data=None,
    x=None,
    y=None,
    order=None,
    sheet="5HT_DL",
):  # , p_values):
    if len(df_significant.index) > 0:
        print(df_significant)
        p_values = df_significant["p-adj"].values
        pairs = [
            (treatment_dict[i[1]["group1"]], treatment_dict[i[1]["group2"]])
            for i in df_significant.iterrows()
        ]

        annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=order)
        # https://stackoverflow.com/questions/64081570/matplotlib-marker-annotation-fontsize-not-shrinking-below-1pt-in-pdf
        annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
        annotator.set_pvalues_and_annotate(p_values)

    return ax


def plotCorrelograms(correlograms):
    fig, axs = plt.subplots(2, 2, figsize=(22, 22))
    axs = flatten(axs)  # Put all the axes into a list at the same level
    for (correlogram_df, p_value_mask, treatment, subvalues), ax in zip(
        correlograms, axs
    ):
        plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax)
    fig.tight_layout()
    return fig


def plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax):
    if np.array_equal(correlogram_df, correlogram_df.T):  # remove duplicate data
        mask = np.triu(np.ones(p_value_mask.shape, dtype=bool), k=1)
        p_value_mask[mask] = True

    heatmap = sns.heatmap(
        correlogram_df,
        vmin=-1,
        vmax=1,
        square=True,
        annot=True,
        cmap="BrBG",
        mask=p_value_mask,
        annot_kws={"size": 6},
        ax=ax,
    )
    ax.set_xticklabels(
        ax.get_xticklabels()
    )  # rotation=45, horizontalalignment='right',
    heatmap.set_title(
        f"{'-'.join(subvalues)} in {treatment}", fontdict={"fontsize": 20}, pad=20
    )
    if len(subvalues) == 1:
        ax.set_ylabel("")
        ax.set_xlabel("")
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
    for group in raw_df.groupby("treatment"):
        result_list.append(
            grubbsTest(raw_df)
        )  # this func will loop through the df to feel grubbsTest
    return result_list


def grubbsTest(group_list):  # include the vairable type in name i.e. group_list series
    # from outliers import smirnov_grubbs as grubbs
    # x_ = grubbs.test(x, alpha=p_value)
    return


######### HEAD TWITCH BEHAVIOR FUNCS ###########
def getRawHTDf(filename):
    return getOrBuildDf(filename, "raw_HT_df", buildRawHTDf)


def buildRawHTDf(filename):
    file_name, file_type = filename.split(".")
    if not file_type == "csv":
        raise Exception(f'METHOD TO DEAL WITH FILE TYPE "{file_type}" ABSENT')
    if not os.path.isfile(f"{INPUT_DIR}/{filename}"):
        raise Exception(f'FILE {filename} IS ABSENT IN "input/" DIRECTORY')
    return pd.read_csv(f"{INPUT_DIR}/{filename}", header=0).replace(
        np.nan, 0
    )  # to set all 0 to Nan


def buildHTHistogram(
    filename,
    HT_filename,
    experiment="agonist_antagonist",
    p_value_threshold=0.05,
    to_plot=["HT_20"],
):
    HT_df = getRawHTDf(HT_filename)
    applyTreatmentMapping(HT_df, filename)

    treatment_mapping = getTreatmentMapping(filename)
    treatment_palette = {
        info["treatment"]: info["color"] for number, info in treatment_mapping.items()
    }
    treatments = [
        treatment_mapping[str(group)]["treatment"]
        for group in getExperimentalInfo(filename)[experiment]["groups"]
    ]
    experimental_df = HT_df[
        HT_df["treatment"].isin(treatments)
    ]  # select only relevent treatments
    columns = to_plot + ["treatment"]  # slice df for plotting
    experimental_df = experimental_df[columns]
    experimental_df = pd.melt(
        experimental_df, id_vars=["treatment"], value_vars=to_plot
    )

    fig, ax = plt.subplots(figsize=(20, 10))
    if len(to_plot) == 1:
        sns.barplot(
            data=experimental_df,
            x="treatment",
            y="value",
            ci=68,
            order=treatments,
            capsize=0.1,
            alpha=0.8,
            palette=treatment_palette,
            errcolor=".2",
            edgecolor=".2",
        )
    else:
        sns.barplot(
            data=experimental_df,
            x="treatment",
            y="value",
            hue="variable",
            ci=68,
            order=treatments,
            capsize=0.1,
            alpha=0.8,
            errcolor=".2",
            edgecolor=".2",
        )

    ax.set_title(f"Head Twitch", y=1.04, fontsize=34)
    ax.set_ylabel("twitches / min", fontsize=24)
    ax.set_xlabel(" ", fontsize=24)
    ax.tick_params(labelsize=24)
    sns.despine(left=False)
    return fig


def getHTHistogram(
    filename,
    HT_filename,
    experiment="dose_responce",
    p_value_threshold=0.05,
    to_plot=[],
    from_scratch=None,
):
    identifier = f"HT_Histogram_{experiment}_for_{to_plot}"
    from_scratch = (
        from_scratch
        if from_scratch is not None
        else input("Recalculate figure even if previous version exists? (y/n)") == "y"
    )
    if from_scratch or not isCached(filename, identifier):
        fig = buildHTHistogram(
            filename,
            HT_filename,
            experiment=experiment,
            p_value_threshold=p_value_threshold,
            to_plot=to_plot,
        )

        cache(filename, identifier, fig)
        saveHTHistogram(fig, identifier)
    else:
        fig = getCache(filename, identifier)
    fig.show()


###### UTILS ########
def flatten(two_d_list):
    return list(itertools.chain.from_iterable(two_d_list))


######## INIT ##########
# Start by checking filesystem has all the folders necessary for read/write operations (cache) or create them otherwise
initiateFileSystem()
