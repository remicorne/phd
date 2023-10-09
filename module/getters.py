import functools
from matplotlib import pyplot as plt

import scipy
import os
from module.statistics import QUANTITATIVE_STAT_METHODS
from module.utils import (
    flatten,
    isCached,
    getCache,
    cache,
    getJSON,
    dictToFilename,
    askYesorNo,
    subselectDf,
)
from module.constants import CACHE_DIR, INPUT_DIR, COLUMN_ORDER
from module.metadata import applyTreatmentMapping
import pandas as pd
import numpy as np


######### GETTERS THAT DEAL WITH THE RAW DATA #########
## TODO Seperate and finalize aggregate stats functions


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


def getMetadata(filename, metadata_type):
    return getJSON(f"{CACHE_DIR}/{filename.split('.')[0]}/{metadata_type}.json")


def getTreatmentMapping(filename):
    return getMetadata(filename, "treatment_mapping")


def getExperimentalInfo(filename):
    return getMetadata(filename, "experimental_info")


def getRegionSubclassification(filename):
    return getMetadata(filename, "region_subclassification")


def getRawDf(filename):
    return getOrBuildDf(filename, "raw_df", buildRawDf)


def getCompoundDf(filename):  # TODO: rename to compounds later
    return getOrBuildDf(filename, "compound_df", buildCompoundDf)


def getCompoundAndRatiosDf(filename):  # TODO: rename to compounds later
    if "compound_and_ratios_df" not in globals():
        global compound_and_ratios_df
        compound_and_ratios_df = getOrBuildDf(
            filename, "compound_and_ratios_df", buildCompoundAndRatiosDf
        )
    return compound_and_ratios_df


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


def getAggregateStatsDf(filename):
    builder_callback_with_param = functools.partial(
        buildAggregateStatsDf, df_type="compound_and_ratios_df"
    )
    return getOrBuildDf(
        filename, "compound_and_ratios_df_aggregate_stats", builder_callback_with_param
    )


def getQuantitativeStats(filename):
    return getOrBuildDf(filename, "quantitative_stats", buildQuantitativeStatsDf)


def getHeadTwitchDf(filename):
    return getOrBuildDf(filename, "headtwitch_df", buildHeadTwitchDf)


########### SETTERS ######


def updateQuantitativeStats(filename, row):
    quantitative_stats_df = getQuantitativeStats(filename)
    quantitative_stats_df = pd.concat(
        [quantitative_stats_df, pd.DataFrame(row)], ignore_index=True
    )
    cache(filename, "quantitative_stats", quantitative_stats_df)
    print("QUANTITATIVE STATS UPDATED")


######### BUILDERS #######


def buildRawDf(filename):
    file_name, file_type = filename.split(".")
    if not file_type == "csv":
        raise Exception(f'METHOD TO DEAL WITH FILE TYPE "{file_type}" ABSENT')
    if not os.path.isfile(f"{INPUT_DIR}/{filename}"):
        raise Exception(f'FILE {filename} IS ABSENT IN "input/" DIRECTORY')
    # to set all 0 to Nan
    return pd.read_csv(f"{INPUT_DIR}/{filename}", header=0).replace(0, np.nan)


# contains the logic to build the df in the new format based on raw df


def buildHeadTwitchDf(HT_filename):
    new_format_df = applyTreatmentMapping(getRawDf(HT_filename), HT_filename)
    return new_format_df


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


# Only compound ratios included here
def buildCompoundAndRatiosDf(filename):
    compound_df = getCompoundDf(filename)  # .iloc[0:100] #To speed up testing
    ratios_df = pd.merge(
        left=compound_df,
        right=compound_df,
        on=[
            "mouse_id",
            "group_id",
            "region",
            "experiment",
            "color",
            "treatment",
        ],
        suffixes=["_1", "_2"],
    ).reset_index(
        drop=True
    )  # merge every compound to every other for each mouse, we want to reset the index (ie make sure it has no duplicates) becaus many update operations will use it

    ratios_df = ratios_df[(ratios_df.compound_1 != ratios_df.compound_2)]

    def calculateRatio(row):
        ratio_name = f"{row.compound_1}/{row.compound_2}"
        print("CALCULATING", row.name, "OF", len(ratios_df), "RATIOS")
        return [
            ratio_name,
            row.value_1 / row.value_2 if row.value_2 else np.NaN,
        ]

    ratios_df[["compound", "value"]] = ratios_df.apply(
        calculateRatio,
        axis=1,
        result_type="expand",
    )  # calculate ratio
    # Drop duplicate columns
    compound_and_ratios_df = pd.concat(
        [
            compound_df,
            ratios_df.drop(columns=["compound_1", "compound_2", "value_1", "value_2"]),
        ]
    )
    return compound_and_ratios_df


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
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    result_ls = []
    for treat_region_comp, groupby_df in compound_and_ratios_df.groupby(
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
            groupby_df.value.values,
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


def buildQuantitativeStatsDf(filename):
    return pd.DataFrame(
        columns=[
            "data_type",
            "experiment",
            "region",
            "compound",
            "test",
            "is_significant",
            "result",
        ]
    )

