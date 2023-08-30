from module.imports import *


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
):
    """
    Ask user histogram parameters and go through the process
    outliers -> stats -> display and save fig

    """
    data = getCompoundAndRatiosDf(filename)
    experimental_info = getExperimentalInfo(filename)
    edit_outliers = True
    while edit_outliers:
        compound = compound or askSelectParameter(data, "compound")
        region = region or askSelectParameter(data, "region")
        experiment = experiment or askMultipleChoice(
            "Select experiment", experimental_info.keys()
        )
        outlier_test = outlier_test or askMultipleChoice(
            "Select outlier test", OUTLIER_TESTS.keys()
        )
        p_value_threshold = p_value_threshold or askMultipleChoice(
            "Select p value threshold", [0.05, 0.01, 0.001, 0.0001]
        )

        processQuantitativeData(
            filename, experiment, compound, region, outlier_test, p_value_threshold
        )
        edit_outliers = askYesorNo("Edit more outliers?")


@get_or_add("histogram")
def processQuantitativeData(
    filename,
    experiment,
    compound,
    region,
    outlier_test,
    p_value_threshold,
):
    """
    Does the whole process of eliminating outliers for a given experiment,
    compound and region, and calculating the stats and figure based on this
    This wil ask you to redo basically the whole process, finalize figure storing
    and handling of multiple deffently parametered figures is the users responsability
    """
    confirmed = False
    while not confirmed:
        data, order, palette = buildHistogramData(
            filename, experiment, compound, region
        )

        if data[f"eliminated_{outlier_test}_outlier"].notna().all():
            if askYesorNo("Redo outlier selection?"):
                data = processOutliers(
                    filename,
                    experiment,
                    compound,
                    region,
                    outlier_test,
                    p_value_threshold,
                )

            is_significant, stats_results = processQuantitativeStats(
                filename, data, p_value_threshold
            )

            if is_significant:
                title = f"{compound} in {region}"
                ylabel = " " if "/" in compound else "ng/mm of tissue"
                fig = buildHistogram(title, ylabel, data, order, palette, stats_results)

            plt.show()

            return fig


def processQuantitativeStats(
    filename, data, experiment, compound, region, p_value_threshold
):
    experiment_info = getExperimentalInfo(filename)[experiment]
    # TODO: replace with independant var and factor logic
    # Implement test passing choice logic

    for test in experiment_info["quantitative_statistics"]:
        is_significant, stats_results = QUANTITATIVE_STAT_METHODS[test](
            data=data,
            independant_vars=experiment_info["independant_vars"],
            p_value_threshold=p_value_threshold,
        )
        if not is_significant:
            if not askYesorNo(f"{test} FAILED, proceed?"):
                break

    return is_significant, stats_results
