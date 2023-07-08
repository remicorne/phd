import functools

import scipy
from module.metadata import applyTreatmentMapping
from module.statistics import *
from module.utils import *
import pandas as pd
import numpy as np

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


def getQuantitativeStats(filename, p_value_threshold, from_scratch=False):
    identifier = "quantitative_stats"
    if from_scratch or not isCached(filename, identifier):
        quantitative_stats_df = buildQuantitativeStatsDf(filename, p_value_threshold)
        cache(filename, identifier, quantitative_stats_df)
    else:
        quantitative_stats_df = getCache(filename, identifier)
    return quantitative_stats_df

def getRawHeadTwitchDf(filename):
    return getOrBuildDf(filename, "raw_head_twitch_df", buildRawHeadTwitchDf)

######### BUILDERS #######


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
    compound_df = getCompoundDf(filename) #.iloc[0:100] #To speed up testing
    ratios_df = pd.merge(left=compound_df, right=compound_df, on=['mouse_id', 'group_id', 'region', 'experiment', 'color', 'treatment'], suffixes=['_1', '_2']) #merge every compound to every other for each mouse
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
    working_df = getCompoundDf(filename) if df_type == 'compound' else getRatiosDf(filename)
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

def buildRawHeadTwitchDf(filename):
    file_name, file_type = filename.split(".")
    if not file_type == "csv":
        raise Exception(f'METHOD TO DEAL WITH FILE TYPE "{file_type}" ABSENT')
    if not os.path.isfile(f"{INPUT_DIR}/{filename}"):
        raise Exception(f'FILE {filename} IS ABSENT IN "input/" DIRECTORY')
    return pd.read_csv(f"{INPUT_DIR}/{filename}", header=0).replace(
        np.nan, 0
    )  # to set all 0 to Nan

