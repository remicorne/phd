from matplotlib import pyplot as plt
import pandas as pd
import scipy
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from module.getters import (
    getCompoundAndRatiosDf,
    getExperimentalInfo,
    getTreatmentMapping,
)
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


def doQuantitativeStatLogic(multiple_factors, multiple_treatments, paired, parametric):
    return {
        (False, False, False, True): ["ttest"],
        (False, False, True, True): ["paired_ttest"],
        (False, True, False, True): ["one_way_anova", "tukey"],
        (False, True, True, True): ["repeated_measures_anova", "paired_ttest"],
        (True, True, False, True): ["two_way_anova", "one_way_anova", "tukey"],
    }[(multiple_factors, multiple_treatments, paired, parametric)]


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
        (is_significant, significance_infos, passed_tests) = processQuantitativeStats(
            experiment_info, data, p_value_threshold
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
    test_results = {}
    for test in tests:
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
        if is_significant in test_results:
            test_results[is_significant].append(test)
        else:
            test_results[is_significant] = [test]
    return is_significant, significance_infos[0], test_results
