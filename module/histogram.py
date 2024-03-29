from matplotlib import pyplot as plt
import numpy as np 
from module.utils import figure_cache
from module.getters import getCompoundAndRatiosDf
import seaborn as sns
from statannotations.Annotator import Annotator
from module.core.Constants import REGIONS, COMPOUNDS

# from brokenaxis import brokenaxis #REMI or future JAS for quantitaiveSummary() modulenotfound in pip list tho

########## GENERIC HISTOGRAM FUNCTIONS MEANT TO BE USED TO BUILD ANY HISTOGRAM #########

########## HISTOGRAM DATA BUILDERS


def buildHistogramData(  # REMI THIS IS NOT SO GENERIC its better in quantitative - I can not use for behavior at all as no compound or region #TODC
    filename,
    experiment,
    compound,
    region,
):
    compound_and_ratios_df = getCompoundAndRatiosDf(
        filename
    )  # this is not the full ratios df, its only intra region compound ratios for nom
    data = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df.compound == compound)
        & (compound_and_ratios_df.region == region)
    ]

    order = data.sort_values(by="group_id", ascending=True).treatment.unique()
    palette = {
        treatment: color
        for treatment, color in data.groupby(by=["treatment", "color"]).groups.keys()
    }

    return data, order, palette




def buildQuantitativeSummaryHistogramData(
    filename, experiment, histogram_type, to_plot, columns
):
    # subselect and transorm to long format
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    value_type = {"compound": "region", "region": "compound"}[histogram_type]

    COLUMN_ORDER = {'region': REGIONS.list, 'compound': COMPOUNDS.list}

    data = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df[value_type].isin(columns))
        & (compound_and_ratios_df[histogram_type] == to_plot)
        & (compound_and_ratios_df.eliminated_grubbs_outlier != True)
    ]
    order = sorted(
        columns, key=lambda x: COLUMN_ORDER[value_type].index(x)
    )  # orders regions / compounds as in constants
    hue_order = data.sort_values(by="group_id", ascending=True).treatment.unique()

    hue_palette = {
        treatment: color
        for treatment, color in data.groupby(by=["treatment", "color"]).groups.keys()
    }

    return data, order, hue_order, hue_palette, value_type


########## HISTOGRAM BUILDERS - can it be made more generic and just one? REMI


def buildHistogram(
    title,
    ylabel,
    data,
    order,
    hue=None,
    palette=None,
    swarm_hue=None,
    swarm_palette=None,
    significance_infos=None,  # x='treatment',y='value'
):
    # JASMINE: in what case would the x and y be variables? #REMI we need to talk about this func as it should be more general
    x = "treatment"
    y = "value"

    fig, ax = plt.subplots(figsize=(20, 10))
    ax = sns.barplot(
        x=x,
        y=y,
        data=data,
        hue=hue,
        palette=palette,
        errorbar=("ci", 68),
        order=order,
        capsize=0.1,
        alpha=0.8,
        errcolor=".2",
        edgecolor=".2",
        dodge=False,
    )
    # #REMI so thiis for the outliers! I was trying to have this function work for my other histogram needs but i cant with this
    ax = sns.swarmplot(
        x=x,
        y=y,
        hue=swarm_hue or hue,
        palette=swarm_palette or palette,
        order=order,
        data=data,
        edgecolor="k",
        linewidth=1,
        linestyle="-",
        dodge=False,
        legend=True if swarm_palette else False,
    )

    if significance_infos:
        ax = labelStats(ax, data, x, y, order, significance_infos)

    ax.tick_params(labelsize=24)
    # ax.set_ylabel(ylabel, fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_xlabel(" ", fontsize=20)  # treatments
    ax.set_title(title, y=1.04, fontsize=34)  # '+/- 68%CI'
    sns.despine(left=False)
    return fig


def buildHueHistogram(
    title,
    ylabel,
    data,
    order,
    x=None,
    y=None,
    hue=None,
    palette=None,
    hue_order=None,
    significance_infos=None,
):
    '''
    Plots histogram of quanatative data for each treatment across region subset. #TODO need to make reversible 
    Args:
        title (string): 
        ylable (string):
        data (pd.DataFrame):                    subset of data - single compound - region subset - experiment
        order (list):                           list of strings corrisponding to x-ticks ['', '', ...]
        x (string):                             column in data of x values
        y (string):                             column in data of y values
        hue (string):                           column name for hue in data  i.e. 'treatment'
        palette (dict):                         dict mapping hue colors {'treatment1':'color1', ... }
        hue_order ( np.array(string) ):         string of each hue in order
        significance_infos (pd.DataFrame):      QuantativeStats for data * post hoc test *  optional pram

    Retuns:
        fig (figure):                           histograms/barplot 
    '''
    fig, ax = plt.subplots(figsize=(20, 10))
    ax = sns.barplot(
        x=x,
        y=y,
        data=data,
        hue=hue,
        palette=palette,
        errorbar=("ci", 68),
        errwidth=1,
        order=order,
        hue_order=hue_order,
        capsize=0.1,
        alpha=0.8,
        errcolor=".2",
        edgecolor=".2",
    )
    # Set the size of x and y ticks
    ax.tick_params(labelsize=16)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.yaxis.set_label_coords(-0.035, 0.5)
    ax.set_xlabel(" ", fontsize=20)  # remove x title
    ax.set_title(title, y=1.04, fontsize=34)
    ax.legend(loc="upper right")  # , bbox_to_anchor=(0.1, 1))
    plt.tight_layout()

    # #handel significance info #already filtered to post hoc (check higher test passes?) # single compound across regions
    
    if significance_infos is not None:
        comparison_hue = hue_order[0] # stats compare against 'vehicles' 

        #loop seaborn bars : xtick1, hue1 --> xtick2, hue1 --> xtick3, hue1   #seaborn v0.11.2
        for i, bar in enumerate(ax.patches):


            # Calculate the index of the group (region) and the hue for this bar
            group_index = i % len(order)
            hue_index = i // len(order)
            # Get the region and hue names based on their indices
            region = order[group_index]
            hue = hue_order[hue_index]

            # IF LOOPING BARS FROM LEFT TO RIGHT THROUGH ALL HUES #if seaborn updates have different ordering of ax.patches
            # # Calculate the number of bars per group (number of hues)
            # bars_per_group = len(hue_order)
            # # Find the index of the group (region) and the hue for this bar
            # group_index = i // bars_per_group
            # hue_index = i % bars_per_group
            # # Get the region and hue names based on their indices
            # region = order[group_index]
            # hue = hue_order[hue_index]
            
            if hue == comparison_hue:
                continue # do not plot stats on control/vehicle

            if region not in significance_infos['region'].unique():
                continue #continue if no stats for that region

            # Check if this hue is in a significant pair for this region
            significant_posthoc_pairs = significance_infos[significance_infos['region'] == region]['p_value'].values[0][0]  #  [ [(hue, hue), (hue,hue)] ,   [p_val, p_val] ]
            significant_posthoc_p_values = significance_infos[significance_infos['region'] == region]['p_value'].values[0][1]
            for pair, p_value in zip(significant_posthoc_pairs, significant_posthoc_p_values):
                if comparison_hue in pair and hue in pair:
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.06, '*', ha='center', va='bottom', fontsize=14)

                    break  

    sns.despine(left=False)
    return fig


def labelStats(ax, data, x, y, order, significance_infos):
    pairs, p_values = significance_infos
    annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=order)
    annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
    annotator.set_pvalues_and_annotate(p_values)

    return ax
