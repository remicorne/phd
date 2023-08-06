from module.imports import *


def grubbsTest(values, p_value_threshold):
    return grubbs.test(values, alpha=p_value_threshold)


def getOutlierLabels(data, test_name, p_value_threshold):
    normal_values = OUTLIER_METHODS[test_name](data.value.values, p_value_threshold)
    return data.apply(lambda row: row.value not in normal_values, axis=1)


def labelTreatmentGroups(subselection_df, test_name, p_value_threshold):
    """calculate outliers for each treatment group and label them in the dataframe"""
    outlier_col_name = f"{test_name}_outlier"
    for treatment in subselection_df.treatment.unique():
        treatment_indices = subselection_df[
            subselection_df.treatment == treatment
        ].index
        outlier_labels = getOutlierLabels(
            subselection_df.loc[treatment_indices], test_name, p_value_threshold
        )
        subselection_df.loc[treatment_indices, outlier_col_name] = outlier_labels

    return subselection_df


def editOutliers(
    filename,
    p_value_threshold,
):
    """Ask user to edit outliers for a given experiment, compound and region"""
    data = getCompoundAndRatiosDf(filename)
    experimental_info = getExperimentalInfo(filename)
    edit_outliers = True
    while edit_outliers:
        compound = askSelectParameter(data, "compound")
        region = askSelectParameter(data, "region")
        experiment = askMultipleChoice("Which experiment", experimental_info.keys())
        experiment_info = experimental_info[experiment]
        test = askMultipleChoice("Which outliers", experiment_info["outliers"])
        processOutliers(filename, experiment, compound, region, test, p_value_threshold)
        edit_outliers = input("Edit more outliers? (y/n)\n") == "y"


def processOutliers(
    filename,
    experiment,
    compound,
    region,
    test_name,
    p_value_threshold,
):
    """Does the whole process of eliminating outliers for a given experiment, compound and region"""
    confirmed = False
    while not confirmed:
        # Get data
        outlier_col_name = f"{test_name}_outlier"
        eliminated_col_name = f"{test_name}_eliminated_outlier"

        data, order, palette = buildHistogramData(
            filename, experiment, compound, region, p_value_threshold
        )
        # Label outliers
        data = labelTreatmentGroups(
            data,
            test_name,
            p_value_threshold,
        )
        # If there are no outliers, exit
        if not data[outlier_col_name].any():
            print("NO OUTLIERS FOUND")
            data[eliminated_col_name] = False
            break

        hue = {f"{test_name}_outlier": {True: "red", False: "black"}}
        # Ask user to codanfirm outliers
        data[eliminated_col_name] = data.apply(
            lambda row: decideOutlier(data, row, order, palette, hue)
            if row[outlier_col_name]
            else False,
            axis=1,
        )
        # Display figure with outliers highlighted
        hue = {eliminated_col_name: {True: "red", False: "black"}}
        fig = buildHistogram(compound, region, data, order, palette, hue=hue)
        plt.show()
        # Ask user to confirm outliers
        confirmed = (
            input(f"Confirm {len(data[data[eliminated_col_name]])} outliers? (y/n)\n")
            == "y"
        )
    # Update cache
    updateCompoundAndRatiosDf(
        filename,
        data,
        [outlier_col_name, eliminated_col_name],
        experiment,
        compound,
        region,
    )
    return fig


def decideOutlier(data, row, order, palette, hue):
    """Ask user to confirm outlier"""
    data[data.mouse_id == row.mouse_id].is_outlier = f"current ({row.mouse_id})"
    title = f"{row.compound} in {row.region}"
    ylabel = " " if "/" in row.compound else "ng/mm of tissue"
    fig = buildHistogram(
        title, ylabel, data, order, palette, hue=hue
    )  # Display figure with current outlier highlighted
    plt.show()
    remove_outliers = input("Remove outlier? (y/n)\n") == "y"
    return remove_outliers


# Get compound and ratios dataframe, update outliers column in the rows that have been edited
# and save it to the cache


def updateCompoundAndRatiosDf(
    filename,
    data,
    columns,
    experiment,
    compound,
    region,
):
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    compound_and_ratios_df.loc[data.index, columns] = data[columns]
    cache(filename, "compound_and_ratios_df", compound_and_ratios_df)
    print(f"{columns} UPDATED FOR {compound} in {region} of {experiment} experiment")


OUTLIER_METHODS = {
    "grubbs": grubbsTest,
}
