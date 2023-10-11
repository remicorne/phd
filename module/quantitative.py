from matplotlib import pyplot as plt
import pandas as pd
import scipy
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from module.getters import (
    getCompoundAndRatiosDf,
    getAggregateStatsDf,
    getExperimentalInfo,
    getTreatmentMapping,
    getRegionSubclassification,
    updateQuantitativeStats,
    getQuantitativeStats,
)
from module.histogram import (
    buildHistogram,
    buildHistogramData,
    buildHueHistogram,
    buildQuantitativeSummaryHistogramData,
)
from module.outliers import OUTLIER_TESTS, processOutliers
from module.statistics import QUANTITATIVE_STAT_METHODS
from module.utils import (
    askMultipleChoice,
    askSelectParameter,
    askYesorNo,
    get_or_add,
    select_params,
)
import pingouin as pg

import matplotlib.patches as mpatches
from module.utils import subselectDf
import seaborn as sns
from module.correlogram import askColumnsToUser


def doQuantitativeStatLogic(multiple_factors, multiple_treatments, paired, parametric):
    return {
        (False, False, False, True): ["ttest"],
        (False, False, True, True): ["paired_ttest"],
        (False, True, False, True): ["one_way_anova", "tukey"],
        (False, True, True, True): ["repeated_measures_anova", "paired_ttest"],
        (True, True, False, True): ["two_way_anova", "one_way_anova", "tukey"],
    }[(multiple_factors, multiple_treatments, paired, parametric)]

def justStats(filename,
    experiments: list,
    compounds: list,
    regions: list,
    p_value_threshold = 0.05,
):
    """User MUST provide all parameters they must be list even if single possibility (iteration on list)
    calculates the stats for every combination of parameters and stores them in the stats df


    Args:
        filename (_type_): name
        experiment (_type_, optional): _description_. Defaults to None.
        compound (_type_, optional): _description_. Defaults to None.
        region (_type_, optional): _description_. Defaults to None.
        p_value_threshold (_type_, optional): _description_. Defaults to None.
    """
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    
    ## Iterat over every combination of parameters
    for experiment in experiments:
        for compound in compounds:
            for region in regions:
                print(f"Calculating stats for {experiment}, {compound}, {region}")
                experiment_info = getExperimentalInfo(filename)[experiment]
                print(compound_and_ratios_df.columns)
                data = compound_and_ratios_df[
                    (compound_and_ratios_df.experiment == experiment)
                    & (compound_and_ratios_df.compound == compound)
                    & (compound_and_ratios_df.region == region)
                ]
                
                if 'eliminated_grubbs_outlier' in data.columns:
                    data = data[data.eliminated_grubbs_outlier != True]
                # the last quantitative test is coded to return the labels directly, thus the need for the bool
                p_value, is_significant, test_results = processQuantitativeStats(
                    experiment_info, data, p_value_threshold
                )

                updateQuantitativeStats(
                    filename,
                    [
                        {
                            "data_type": "HPLC",
                            "experiment": experiment,
                            "compound": compound,
                            "region": region,
                            **test_result,
                        }
                        for test_result in test_results
                    ],
                )
                
    print('DONE')

    
    
            

def quantitativeHistogram(
    filename,
    experiment=None,
    compound=None,
    region=None,
    outlier_test=None,
    p_value_threshold=None,
    from_scratch=None,
):
    """
    Ask user histogram parameters and go through the process
    outliers -> stats -> display and save fig

    """
    data = getCompoundAndRatiosDf(filename)
    experimental_info = getExperimentalInfo(filename)
    exit_loop = False
    while not exit_loop:
        # We use keyword params here even though they actually are mandatory for the decorator
        singleQuantitativeHistogram(
            filename,
            experiment=experiment
            or askMultipleChoice("Select experiment", experimental_info.keys()),
            compound=compound or askSelectParameter(data, "compound"),
            region=region or askSelectParameter(data, "region"),
            outlier_test=outlier_test
            or askMultipleChoice("Select outlier test", OUTLIER_TESTS.keys()),
            p_value_threshold=p_value_threshold
            or askMultipleChoice(
                "Select p value threshold", [0.05, 0.01, 0.001, 0.0001]
            ),
            from_scratch=from_scratch,
        )
        exit_loop = askYesorNo("Exit?")


# We use keyword params here even though they actually are mandatory for the decorator
@get_or_add("histogram")
def singleQuantitativeHistogram(
    filename,
    experiment=None,
    compound=None,
    region=None,
    outlier_test=None,
    p_value_threshold=None,
    from_scratch=None,
):
    """
    Does the whole process of eliminating outliers for a given experiment,
    compound and region, and calculating the stats and figure based on this
    This wil ask you to redo basically the whole process, finalize figure storing
    and handling of multiple deffently parametered figures is the users responsability
    """
    confirmed = False
    experiment_info = getExperimentalInfo(filename)[experiment]
    eliminated_outlier_col_name = f"eliminated_{outlier_test}_outlier"

    while not confirmed:
        data, order, palette = buildHistogramData(
            filename, experiment, compound, region
        )

        if (
            eliminated_outlier_col_name not in data
            or data[eliminated_outlier_col_name].isna().any()
            or askYesorNo("Redo outlier selection?")
        ):
            data = processOutliers(
                filename,
                experiment,
                compound,
                region,
                outlier_test,
                p_value_threshold,
            )

        data = data[data[eliminated_outlier_col_name] == False]
        # the last quantitative test is coded to return the labels directly, thus the need for the bool
        is_significant, significance_infos, test_results = processQuantitativeStats(
            experiment_info, data, p_value_threshold
        )

        updateQuantitativeStats(
            filename,
            [
                {
                    "data_type": "HPLC",
                    "experiment": experiment,
                    "compound": compound,
                    "region": region,
                    **test_result,
                }
                for test_result in test_results
            ],
        )

        # JJB ok for the title would like to have either:  " passes: twowayANOVA, onewayANOVA " oder "failed: two-way-anova"
        title = f"{compound} in {region}"
        ylabel = " " if "/" in compound else "ng/mg of tissue"

        fig = buildHistogram(
            title,
            ylabel,
            data,
            order,
            palette,
            significance_infos=significance_infos if is_significant else None,
        )

        plt.show()

        return fig


def processQuantitativeStats(experiment_info, data, p_value_threshold):
    """_summary_

    Args:
        experiment_info (dict): experiment info mapping
        data (df): data for signle experiment
        p_value_threshold (float): threshold applied to all tests

    Returns:
        is_significant (bool), significance_infos (list, list): [treatment parings, pvalues]
    """
    multiple_factors = len(experiment_info["independant_vars"]) > 1
    multiple_treatments = len(experiment_info["groups"]) > 1
    paired = experiment_info["paired"]
    parametric = experiment_info["parametric"]

    tests = doQuantitativeStatLogic(
        multiple_factors, multiple_treatments, paired, parametric
    )
    test_results = []
    for test in tests:
        p_value, is_significant, stats_results = QUANTITATIVE_STAT_METHODS[
            test
        ](
            data=data[data.value.notna()],
            independant_vars=experiment_info["independant_vars"],
            p_value_threshold=p_value_threshold,
        )
        print()
        print(test.upper(), "SIGNIFICANT" if is_significant else "NOT SIGNIFICANT")
        print(stats_results)
        test_results.append(
            {"test": test, "is_significant": is_significant, "result": stats_results, "p_value": p_value}
        )

    return is_significant, p_value, test_results


def quantitativeSummary(
    filename,
    experiment=None,
    histogram_type=None,  #  compound or region
    to_plot=None,  # chosen compound or region 
    p_value_threshold=None,
    columns=None, # x values to plot 
    from_scratch=None,
):
    # get df
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    experiments = getExperimentalInfo(filename)

    # JASI input loop here if you like see quantitativeHistogram()
    # prompts for inputs if not provided
    experiment = (
        experiment
        if experiment
        else askMultipleChoice(
            "Which experiment?",
            experiments.keys(),
        )
    )
    histogram_type = (
        histogram_type
        if histogram_type
        else askMultipleChoice("Which histogram?", ["compound", "region"])
    )
    to_plot = (
        to_plot
        if to_plot
        else input(
            f"""Which {histogram_type}?
                    Possibilities: {set(compound_and_ratios_df[histogram_type])}
                    """
        ).upper()
    )
    columns = (
        columns if columns else askColumnsToUser(histogram_type, compound_and_ratios_df)
    )

    buildSingleQuantitativeSummary(
        filename,
        experiment=experiment,
        histogram_type=histogram_type,
        to_plot=to_plot,
        columns=columns,
        from_scratch=from_scratch,  # Used in decorator
    )


@get_or_add("quantitative_summary")
def buildSingleQuantitativeSummary(
    filename,
    experiment=None,
    histogram_type=None, # e.g. compound 
    to_plot=None, # DA
    columns=None, # region list
    from_scratch=None,  # Used in decorator
):
    title = ""  # EMPTY to feed generic plotter
    ylabel = f"{to_plot} {'ng/mg' if '/' not in to_plot else ''} +/-98CI"

    (
        data,
        order,
        hue_order,
        hue_palette,
        value_type,
    ) = buildQuantitativeSummaryHistogramData(
        filename, experiment, histogram_type, to_plot, columns
    )

    #create a significance mapping for the x values #RN this is T/F but will be corrected to pvalues and stars later 
    if histogram_type == 'compound':
        compound = to_plot
        region = columns
    else:
        compound = columns
        region = to_plot #this allows the use of the function both ways but feels fucking stupid tbh

    #TO DO pretty sure this should be its own function at some point statsTests(filename, experiment) also in processQuantitativeStats
    experiment_info = getExperimentalInfo(filename)[experiment]
    multiple_factors = len(experiment_info["independant_vars"]) > 1
    multiple_treatments = len(experiment_info["groups"]) > 1
    paired = experiment_info["paired"]
    parametric = experiment_info["parametric"]

    tests = doQuantitativeStatLogic(
        multiple_factors, multiple_treatments, paired, parametric
    )
    test = tests[0]

    #get stats #REMI TODO missing p_value col for non  post hoc testing
    quant_stats_df = subselectDf(getQuantitativeStats(filename), {'experiment':experiment, 'compound':compound, 'region': region, 'test':test })

    #  buildHistogram() IS NOT GENERAL enough becaue of hue/outlier stuff - REMI can we combine? I tried and couldnt manage ... yet
    #hue hard coded to treatment for quantitativeSummary()
    fig = buildHueHistogram(
        title,
        ylabel,
        data,
        order,
        palette=hue_palette,
        x=value_type,
        y="value",
        hue="treatment",
        hue_order=hue_order,
        significance_infos=quant_stats_df,
    )


    return fig


@get_or_add(
    "percentage_vehicles"
)  # REMI this should probebly be singlePercentageVehiclesFig
def percentageVehiclesFig(
    filename,
    experiment=None,
    compound=None,
    regions=None,  # REMI i would ike this to work the same way it does for correlograms i.e. also specifying the order
    from_scratch=True,
):
    # slice aggstats df
    experimental_df = subselectDf(
        getAggregateStatsDf(filename), {"compound": compound, "experiment": experiment}
    )  # REMI can use completly omnce has list comprehension
    experimental_df = experimental_df[experimental_df["region"].isin(regions)]

    # check that COLUMN_ORDER is consistent with region_subclasification i.e. CORTEX then SUBCORT...
    region_order = (
        regions  # TODO check that regions inputted are sorted by subclasification
    )
    experimental_df = experimental_df.assign(
        region=lambda x: pd.Categorical(
            x["region"], categories=region_order, ordered=True
        )
    ).sort_values("region")

    # create a new column % of control mean/control_mean * 100
    experimental_df.loc[:, "percentage_of_vehicles"] = experimental_df.groupby(
        "region"
    )["mean"].transform(
        lambda x: (x / x.loc[experimental_df["treatment"] == "vehicles"].values[0])
        * 100
    )

    # melt df
    plot_experimental_df = pd.melt(
        experimental_df,
        id_vars=["region", "treatment"],
        value_vars=["percentage_of_vehicles"],
    )

    # palete
    treatment_palette = {
        info["treatment"]: info["color"]
        for number, info in getTreatmentMapping(filename).items()
    }

    # BUILD
    fig, ax = plt.subplots(figsize=(12, 9))
    sns.set_style("white")
    sns.set_context("notebook")

    # plot lines for each treatment
    for treatment_value in plot_experimental_df["treatment"].unique():
        treatment_data = plot_experimental_df[
            plot_experimental_df["treatment"] == treatment_value
        ]
        ax.plot(
            treatment_data["region"],
            treatment_data["value"],
            marker="o",
            markersize=2,
            linestyle="-",
            label=treatment_value,
            color=treatment_palette.get(
                treatment_value, None
            ),  # Use the palette if needed
        )

    # Create custom legend handles for 'treatment' with black outline
    order = [
        getTreatmentMapping(filename)[str(group)]["treatment"]
        for group in getExperimentalInfo(filename)[experiment]["groups"]
    ]
    legend_handles = [
        mpatches.Patch(
            facecolor=treatment_palette[treatment], edgecolor="black", label=treatment
        )
        for treatment in order
    ]
    # Add a legend with custom handles and formatting
    ax.legend(handles=legend_handles, loc="upper left", frameon=False)

    # Add spans for region sub-classifications
    for subcls, info in getRegionSubclassification(filename).items():
        regions_to_highlight = info["regions"]
        color = info["color"]

        # Find the indices of regions to highlight
        indices_to_highlight = [
            region_order.index(region) for region in regions_to_highlight
        ]

        # Create a single span for the highlighted regions
        start_index = indices_to_highlight[0]
        end_index = indices_to_highlight[-1]
        ax.axvspan(
            start_index - 0.5,
            end_index + 0.5,
            facecolor=color,
            alpha=0.1,
            edgecolor="none",
            label=subcls,
        )

        # Add a label for the region sub-classification at the midpoint of the span
        label = subcls.replace("_", " ")  # Replace underscores with spaces
        label_x = (start_index + end_index) / 2
        label_y = ax.get_ylim()[1] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
        # Split label into lines if it contains multiple words
        words = label.split()
        lines = ["\n".join(words[i : i + 2]) for i in range(0, len(words), 2)]
        label_text = "\n".join(lines)
        ax.text(
            label_x,
            label_y,
            label_text,
            ha="center",
            va="top",
            fontsize=12,
            color=color,
            bbox=None,  # {'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.8}  # Add a white background for visibility
        )

    ax.set_xlim(-0.5, len(region_order) - 0.5)  # Set x-axis limits to edge of axspan
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ax.set_xlabel("Regions")
    ax.set_ylabel(f"{compound} as a % of vehicles")

    return fig
