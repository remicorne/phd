from matplotlib import pyplot as plt
import pandas as pd
import scipy
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from module.getters import getCompoundAndRatiosDf, getExperimentalInfo, getTreatmentMapping 
from module.histogram import buildHistogram, buildHistogramData
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
        singleQuantitativeHistogram( #REMI i think this is where the endless prompts come fom re hist 
            filename,
            experiment=askMultipleChoice("Select experiment", experimental_info.keys()),
            compound=askSelectParameter(data, "compound"),
            region=askSelectParameter(data, "region"),
            outlier_test=askMultipleChoice("Select outlier test", OUTLIER_TESTS.keys()),
            p_value_threshold=askMultipleChoice(
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
        (
            is_significant,
            significance_infos,
        ) = processQuantitativeStats(experiment_info, data, p_value_threshold)

        #JJB ok for the title its nice to show the stats that are significant/shown by stars not all tests done as it was
        #i.e. if its two factor (agonist antagonist) and two way anova is calculated you have this choice to continue even if you failed the two way in this case perhaps we add '(failed two-way-anova)' to title
        title = f"{compound} in {region} " #{experiment_info['quantitative_statistics'].keys()}
        ylabel = " " if "/" in compound else "ng/mm of tissue"

        
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
    # TODO: replace with independant var and factor logic
    # Implement test passing choice logic

    for test in experiment_info["quantitative_statistics"]:
        is_significant, stats_results, *significance_infos = QUANTITATIVE_STAT_METHODS[
            test
        ](
            data=data,
            independant_vars=experiment_info["independant_vars"],
            p_value_threshold=p_value_threshold,
        )
        print()
        print(test.upper(), "SIGNIFICANT" if is_significant else "NOT SIGNIFICANT")
        print(stats_results)
        if not is_significant:
            if not askYesorNo(f"{test} FAILED, proceed?"):
                return False, None

    return is_significant, significance_infos[0]
