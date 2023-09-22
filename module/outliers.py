from outliers import smirnov_grubbs as grubbs
import matplotlib.pyplot as plt
from module.getters import getCompoundAndRatiosDf
from module.histogram import buildHistogram, buildHistogramData
from module.utils import askYesorNo, cache

########### FUNCTIONS THAT DEAL WITH THE OUTLIER PROCESSING LOGIC ###########


def grubbsTest(values, p_value_threshold):
    """
    Takes a list of values on which to perform the test and returns normal values
    """
    return grubbs.test(values, alpha=p_value_threshold)


def getOutliers(data, outlier_test, p_value_threshold):
    """
    Takes a series containing the values and returns a series with T/F outlier labels
    The outlier test called must return the list of normal values for this to work
    """
    normal_values = OUTLIER_TESTS[outlier_test](data.value.values, p_value_threshold)
    return data.apply(lambda row: row.value not in normal_values, axis=1)


def labelOutliers(subselection_df, outlier_test, p_value_threshold):
    """
    calculate outliers for each treatment group and label them in the dataset
    Return labelled dataset
    """
    outlier_col_name = f"{outlier_test}_outlier"
    for treatment in subselection_df.treatment.unique():
        treatment_indices = subselection_df[
            subselection_df.treatment == treatment
        ].index
        outlier_labels = getOutliers(
            subselection_df.loc[treatment_indices], outlier_test, p_value_threshold
        )
        subselection_df.loc[treatment_indices, outlier_col_name] = outlier_labels

    return subselection_df


def processOutliers(
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
    """
    confirmed = False
    while not confirmed:
        outlier_col_name = f"{outlier_test}_outlier"
        eliminated_outlier_col_name = f"eliminated_{outlier_col_name}"

        data, order, palette = buildHistogramData(
            filename, experiment, compound, region
        )
        # Label outliers
        data = labelOutliers(
            data,
            outlier_test,
            p_value_threshold,
        )
        # If there are no outliers, exit
        if not data[outlier_col_name].any():
            print("NO OUTLIERS FOUND")
            data.loc[data.index, eliminated_outlier_col_name] = False
            break

        title = f"{compound} in {region}"
        ylabel = " " if "/" in compound else "ng/mm of tissue"
        hue = {eliminated_outlier_col_name: {True: "red", False: "white"}}

        # Ask user to codanfirm outliers
        data.loc[:, eliminated_outlier_col_name] = False
        data[data[outlier_col_name]].apply(
            lambda row: decideOutlier(
                data,
                row,
                order,
                outlier_col_name,
                eliminated_outlier_col_name,
                title,
                ylabel,
                palette,
            )
            if row[outlier_col_name]
            else False,
            axis=1,
        )
        # Display figure with outliers highlighted
        fig = buildHistogram(title, ylabel, data, order, palette, hue=hue)
        plt.show()
        # Ask user to confirym outliers
        confirmed = askYesorNo(
            f"Confirm {len(data[data[eliminated_outlier_col_name]])} outliers?"
        )
    # Update cache
    updateCompoundAndRatiosDf(
        filename,
        data,
        [outlier_col_name, eliminated_outlier_col_name],
        experiment,
        compound,
        region,
    )
    return data


def decideOutlier(
    data,
    row,
    order,
    outlier_col_name,
    eliminated_outlier_col_name,
    title,
    ylabel,
    palette,
):
    """Ask user to confirm outlier"""
    current_label = f"current (id={row.mouse_id})"
    hue = {outlier_col_name: {True: "white", False: "black", current_label: "red"}}
    data.loc[row.name, outlier_col_name] = current_label
    fig = buildHistogram(
        title, ylabel, data, order, palette, hue=hue
    )  # Display figure with current outlier highlighted
    data.loc[row.name, outlier_col_name] = True
    plt.show()
    data.loc[row.name, eliminated_outlier_col_name] = askYesorNo("Remove outlier?")


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


OUTLIER_TESTS = {
    "grubbs": grubbsTest,
}
