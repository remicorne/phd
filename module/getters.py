import functools
from matplotlib import pyplot as plt
import seaborn as sns
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
    saveQuantitativeSummaryFig,
)
from module.constants import CACHE_DIR, INPUT_DIR, COLUMN_ORDER
from module.metadata import applyTreatmentMapping
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches

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
        buildAggregateStatsDf, df_type='compound_and_ratios_df'
    )
    return getOrBuildDf(
        filename, "compound_and_ratios_df_aggregate_stats", builder_callback_with_param
    )


def getQuantitativeStats(filename, p_value_threshold, from_scratch=False):
    identifier = "quantitative_stats"
    if from_scratch or not isCached(filename, identifier):
        quantitative_stats_df = buildQuantitativeStatsDf(filename, p_value_threshold)
        cache(filename, identifier, quantitative_stats_df)
    else:
        quantitative_stats_df = getCache(filename, identifier)
    return quantitative_stats_df


def getHeadTwitchDf(filename):
    return getOrBuildDf(filename, "headtwitch_df", buildHeadTwitchDf)

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

def buildHeadTwitchDf(HT_filename):
    new_format_df=applyTreatmentMapping(getRawDf(HT_filename), HT_filename)
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


# Change experiment param to none when experiment info mapping is complete
def buildQuantitativeStatsDf(
    filename, subselect=None, experiment="agonist_antagonist", p_value_threshold=0.05
):
    experiment_info = getExperimentalInfo(filename)[experiment]
    compound_and_ratios_df = (
        subselectDf(getCompoundAndRatiosDf(filename), subselect)
        if subselect
        else getCompoundAndRatiosDf(filename)
    )
    compound_and_ratios_df = compound_and_ratios_df[
        compound_and_ratios_df.experiment == experiment
    ]  # TODO: Remove when info mapping is completed
    data_keys = ["experiment", "compound", "region"]
    groupbys = compound_and_ratios_df.groupby(by=data_keys)

    return pd.DataFrame(
        list(
            map(
                buildGroupbyQuantitativeStatsMapper(
                    experiment_info, len(groupbys), p_value_threshold
                ),
                list(enumerate(groupbys)),
            )
        ),
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


def buildGroupbyQuantitativeStatsMapper(
    experiment_info, groupby_length, p_value_threshold
):
    def executor(i_groupby):
        i, groupby = i_groupby
        print("CALULATING", i + 1, "OF", groupby_length, "GROUPBYS")
        experiment_compound_region, data = groupby
        experiment, compound, region = experiment_compound_region
        return buildGroupbyQuantitativeStats(
            experiment_info, experiment, compound, region, data, p_value_threshold
        )

    return executor


def buildGroupbyQuantitativeStats(
    experiment_info, experiment, compound, region, data, p_value_threshold
):
    return [
        experiment,
        compound,
        region,
    ] + flatten(
        [
            QUANTITATIVE_STAT_METHODS[test](
                data=data,
                independant_vars=experiment_info["independant_vars"],
                p_value_threshold=p_value_threshold,
            )
            for test in experiment_info["quantitative_statistics"]
        ]
    )

#replacing in the quantitative.py as only builder should be here no?
def getQuantitativeSummaryFig(
    filename,
    experiment="dose_response",
    value_type="ratio",
    value="5HIAA/5HT",
    regions_to_plot=COLUMN_ORDER,
    from_scratch=None,
):
    identifier = (
        f"{experiment}_for_{value.replace('/', ':')}_{(',').join(regions_to_plot)}"
    )
    from_scratch = (
        from_scratch
        if from_scratch is not None
        else askYesorNo("Recalculate figure even if previous version exists?")
    )
    if from_scratch or not isCached(filename, identifier):
        fig = buildQuantitativeSummaryFig(
            filename,
            experiment=experiment,
            compound=value,
            regions_to_plot=regions_to_plot,
        )
        cache(filename, identifier, fig)
        saveQuantitativeSummaryFig(fig, identifier)
    else:
        fig = getCache(filename, identifier)
    fig.show()



# @get_or_add("quantitative_summary")
# def singleQuantitativeSummaryFig(

def buildQuantitativeSummaryFig(
    filename,
    experiment="dose_response",
    compound="5HIAA/5HT",
    regions=COLUMN_ORDER,
):
    # FIX ME: build in scafolding so if experiment is not in experiments in treatmentmapping then raise error: use 'dose_response' or 'agonist_antagonist'
    # REMI: i guess you will also want me to modularise this, and we need to unify naming for different plots
    # get aggstats df
    # compound_type= 'ratio' if '/' in compound else 'compound'#redundant with new agg stats getter?
    values_df = getAggregateStatsDf(filename)

    # slice df to experiment and vale
    treatment_mapping = getTreatmentMapping(filename)
    experimental_df = values_df[
        (values_df["region"].isin(regions))
        & (values_df["experiment"] == experiment)
        & (values_df["compound"] == compound)
    ]  # select only relevent treatments AND REGIONS

    # create a new column % of control mean/control_mean * 100
    experimental_df.loc[:, "percentage_of_vehicles"] = experimental_df.groupby(
        "region"
    )["mean"].transform(
        lambda x: (x / x.loc[experimental_df["treatment"] == "vehicles"].values[0])
        * 100
    )

    # order and reshape df to plot
    experimental_df = (
        experimental_df.loc[experimental_df["region"].isin(regions)]
        .assign(
            region=lambda x: pd.Categorical(
                x["region"], categories=regions, ordered=True
            )
        )
        .sort_values("region")
    )

    plot_experimental_df = pd.melt(
        experimental_df,
        id_vars=["region", "treatment"],
        value_vars=["percentage_of_vehicles"],
    )

    # load pallette and open fig
    treatment_palette = {
        info["treatment"]: info["color"] for number, info in treatment_mapping.items()
    }


    fig, ax = plt.subplots(figsize=(12, 9))
    # sns.set_style("whitegrid")
    sns.set_style("white")
    sns.set_context("notebook")

    # plot lines
    sns.lineplot(
        data=plot_experimental_df,
        x="region",
        y="value",
        hue="treatment",
        palette=treatment_palette,
    )

    # add markers for each region
    marker_mapping = {
        value["treatment"]: value["markers"]
        for value in treatment_mapping.values()
        if experiment in value["experiments"]
    }

    # Define the minimum and maximum marker sizes
    min_marker_size = 20
    max_marker_size = 100

    # Calculate the marker sizes based on the difference from 100
    plot_experimental_df["marker_size"] = abs(plot_experimental_df["value"] - 100)

    # Normalize marker sizes between the minimum and maximum sizes
    plot_experimental_df["marker_size"] = (
        (
            plot_experimental_df["marker_size"]
            - plot_experimental_df["marker_size"].min()
        )
        / (
            plot_experimental_df["marker_size"].max()
            - plot_experimental_df["marker_size"].min()
        )
    ) * (max_marker_size - min_marker_size) + min_marker_size

    # Plot with adjusted marker sizes
    sns.scatterplot(
        data=plot_experimental_df,
        x="region",
        y="value",
        hue="treatment",
        palette=treatment_palette,
        style="treatment",
        markers=marker_mapping,
        legend=False,
        size="marker_size",
        sizes=(min_marker_size, max_marker_size),  # Set the range of marker sizes
        alpha=0.7,  # Adjust the marker transparency if desired
    )

    # Add axvspan to each grouping
    region_subclassification = getRegionSubclassification(filename)
    region_subclassification = {
        group.replace("_", " "): properties
        for group, properties in region_subclassification.items()
    }
    for group, properties in region_subclassification.items():
        regions = properties["regions"]
        color = properties["color"]
        label = group.replace("_", " ")
        if regions:
            start_idx = regions_to_plot.index(regions[0]) - 0.5
            end_idx = regions_to_plot.index(regions[-1]) + 0.5
            ax.axvspan(start_idx, end_idx, facecolor=color, alpha=0.2, label=label)
            # Adjust x-axis limits
            ax.set_xlim(-0.5, len(regions_to_plot) - 0.5)  # Set the x-axis limits

            # Add axvspan label

            label_x = (start_idx + end_idx) / 2
            label_y = ax.get_ylim()[1] - 0.05 * (
                ax.get_ylim()[1] - ax.get_ylim()[0]
            )  # Adjust the label_y coordinate
            words = label.split()
            lines = [
                words[i : i + 2] for i in range(0, len(words), 2)
            ]  # Split words into pairs for multiline display
            line_height = 18  # Adjust the height between lines

            for i, line in enumerate(lines):
                line_label = "\n".join(line)  # Join words with a line break
                current_label_y = label_y - i * line_height
                ax.annotate(
                    line_label,
                    xy=(label_x, current_label_y),
                    xycoords="data",
                    xytext=(0, -12),
                    textcoords="offset points",
                    ha="center",
                    va="top",
                    color=color,
                    fontsize=18,
                )

    # Set x-ticks and labels
    #    compound_type= 'ratio' if '/' in compound else 'compound'#redundant with new agg stats getter?

    if value_type == "ratio":
        y_label = "ratio"
    elif value_type == "compound":
        y_label = "ng/mg"
    # Set x-tick positions and labels
    x_ticks = range(len(regions_to_plot))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(regions_to_plot, rotation=45)

    # ax.set_xticklabels(regions_to_plot, rotation=45)
    # ax.tick_params(axis='x', rotation=45)
    ax.set_xlabel("Region")
    ax.set_ylabel(f"{y_label} {value}")
    ax.set_title(
        f'{experiment.replace("_", " ").capitalize()} for {value} (% of vehicles)',
        fontsize=20,
        pad=20,
    )

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Create custom legend handles for 'treatment' with black outline in the specified order
    order = [
        treatment_mapping[str(group)]["treatment"]
        for group in getExperimentalInfo(filename)[experiment]["groups"]
    ]
    legend_handles = [
        mpatches.Patch(
            facecolor=treatment_palette[treatment], edgecolor="black", label=treatment
        )
        for treatment in order
    ]

    # Add a legend with black outline
    plt.legend(handles=legend_handles, loc="upper left", frameon=False)

    # Adjust the spacing
    plt.tight_layout()
    return fig
