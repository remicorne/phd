from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
from module.getters import (
    getExperimentalInfo,
    getTreatmentMapping,
    getCompoundAndRatiosDf,
)
from module.utils import askMultipleChoice, flatten, get_or_add, inputEscape
from module.statistics import getPearsonR, getPearsonPValue, isSignificant
from module.constants import getCorrelogramColumns, CORRELOGRAM_TYPE_CONVERTER
import pandas as pd

# TODO: Handle columns and more generality


def getAndPlotMultipleCorrelograms(
    filename, selector, p_value_threshold=0.05, n_minimum=5, from_scratch=None
):
    experiment = selector.pop("experiment", None)
    for (
        correlogram_type,
        to_correlate_list,
    ) in selector.items():  # Iterate through the selector dict
        for to_correlate in to_correlate_list:
            correlogram()(
                filename,
                experiment,
                correlogram_type,
                to_correlate,
                p_value_threshold,
                n_minimum,
                from_scratch,
            )


##TOOD: make more interactive with loop for user to be able to calibrate correlogram as desired
def correlogram(
    filename,
    experiment=None,
    correlogram_type=None,
    to_correlate=None,
    p_value_threshold=0.05,
    n_minimum=5,
    columns=None,
    from_scratch=None,
):
    """
    This is the function that is called by the user, it will call buildSingleCorrelogram that builds a signle correlogram
    The function can be called with parameters, user input will be required if not
    """
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    experiments = getExperimentalInfo(filename)
    experiment = (
        experiment
        if experiment
        else askMultipleChoice(
            "Which experiment?",
            experiments.keys(),
        )
    )
    correlogram_type = (
        correlogram_type
        if correlogram_type
        else askMultipleChoice("Which correlogram?", ["compound", "region"])
    )
    to_correlate = (
        to_correlate
        if to_correlate
        else input(
            f"""Which {correlogram_type}?
                    (Enter simple {correlogram_type} or {correlogram_type}/{correlogram_type} for ratio,
                    Use {correlogram_type} or {correlogram_type}/{correlogram_type} for simple correlogram, or {correlogram_type}-{correlogram_type} or {correlogram_type}/{correlogram_type}-{correlogram_type} for square correlogram
                    Possibilities: {set(compound_and_ratios_df[correlogram_type])}
                    """
        ).upper()
    )
    columns = (
        columns if columns 
        else askColumnsToUser(correlogram_type, compound_and_ratios_df) #REMI: for me I would like to be able to manuely code these in pls ass the if not assigned, danke
    )

    buildSingleCorrelogram(
        filename,
        experiment=experiment,
        correlogram_type=correlogram_type,
        to_correlate=to_correlate,
        p_value_threshold=p_value_threshold,
        n_minimum=n_minimum,
        columns=columns,
        from_scratch=from_scratch,
    )


@get_or_add("correlogram")
def buildSingleCorrelogram(
    filename,
    experiment,
    correlogram_type,
    to_correlate,
    p_value_threshold,
    n_minimum,
    columns,
    from_scratch,  # Used in decorator
):
    # this is not the full ratios df, its only intra region compound ratios for nom
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    value_type = {"compound": "region", "region": "compound"}[correlogram_type]
    to_correlate = to_correlate.split("-")
    subselection_df = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df[value_type].isin(columns))
    ].query("|".join([f"{correlogram_type}=='{value}'" for value in to_correlate]))
    columns = columns #if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]#REMI this doesnt exist anywhere?
    pivot_columns = {
        "region": ["region", "compound"],
        "compound": ["compound", "region"],
    }[correlogram_type]
    # Because table is înveted on region and compound even in the case of simple correlogram I have to duplicate the selector in that case to avoid ugly labelling
    pivot_column_value = (
        to_correlate if len(to_correlate) == 2 else to_correlate + to_correlate
    )
    correlograms = []
    treatment_mapping = getTreatmentMapping(filename)
    order = [
        treatment_mapping[str(group)]["treatment"]
        for group in getExperimentalInfo(filename)[experiment]["groups"]
    ]
    subselection_df_ordered = subselection_df.iloc[
        subselection_df["treatment"]
        .astype(pd.CategoricalDtype(categories=order))
        .argsort()
    ]  # order df in order to plot, then sort=False in .groupby
    for treatment, group_df in subselection_df_ordered.groupby(
        by=["treatment"], sort=False
    ):
        # Here we just picot our structure to one where each column is a region or compound and each line a mouse    correlogram_df, p_value_mask = [pivot_df.corr(method=method, min_periods=n_minimum).dropna(axis=0, how='all').dropna(axis=1, how='all') for method in methods] # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
        pivot_df = group_df.pivot_table(
            values="value", index=group_df["mouse_id"], columns=pivot_columns
        )
        # This is the list of methods to pass to df.corr. I know you used to pass 'pearson' as a string but this way the R calculation and pvalue mask are 'linked' by being done on the same line
        methods = [getPearsonR, isSignificant(getPearsonPValue, p_value_threshold)]
        # order columns in desired plotting order
        columns = columns #if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type] #REMI i am pretty sure this is old code no?
        pivot_columns_ordered = sorted(
            pivot_df.columns,
            key=lambda x: columns.index(x[1]) if x[1] in columns else float("inf"),
        )
        pivot_df_ordered = pivot_df[pivot_columns_ordered]
        correlogram_df, p_value_mask = [
            pivot_df_ordered.corr(method=method, min_periods=n_minimum)
            .loc[tuple(pivot_column_value)]
            .dropna(axis=0, how="all")
            .dropna(axis=1, how="all")
            for method in methods
        ]  # ANd finally we calculate the R values and the pvalue mask in a for loop becaus the go through the exact same treatment
        correlograms.append(
            [correlogram_df, p_value_mask.astype(bool), treatment[0], to_correlate]
        )
    fig = plotCorrelograms(correlograms)
    return fig


def plotCorrelograms(correlograms):
    #JJB might be nice to have one color bar for all figures
    fig, axs = plt.subplots(2, 2, figsize=(22, 22))
    axs = flatten(axs)  # Put all the axes into a list at the same level
    for (correlogram_df, p_value_mask, treatment, subvalues), ax in zip(
        correlograms, axs
    ):
        plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax)
    fig.tight_layout()
    return fig


def plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax):
    

    if np.array_equal(correlogram_df, correlogram_df.T):  # TRIANGLE CORRELOGRAMS remove duplicate data  
        title = ax.set_title(f"{'-'.join(subvalues)} in {treatment}", fontsize=28, pad=20, y=0.9)  # Adjust the y position of the title manually #JJB set

        mask = np.triu(np.ones(p_value_mask.shape, dtype=bool), k=1)
        p_value_mask[mask] = True
        np.fill_diagonal(
            p_value_mask.values, False
        )  # this makes the diagonal correlations of 1 visible
    else:
        title = ax.set_title(f"{'-'.join(subvalues)} in {treatment}", fontsize=28, pad=20, y=1)  # Adjust the y position of the title manually for square correlogram

    heatmap = sns.heatmap(
        correlogram_df,
        vmin=-1,
        vmax=1,
        square=True,
        annot=True,
        cmap=matplotlib.colormaps['BrBG'],
        mask=p_value_mask,
        annot_kws={"size": 8}, #, 'fontweight':'bold'
        ax=ax,
        cbar_kws={"shrink": 0.8} #adj color bar size
    )
    ax.set_xticklabels(
        ax.get_xticklabels()
    )  # rotation=45, horizontalalignment='right',

    if len(subvalues) == 1:
        ax.set_ylabel("")
        ax.set_xlabel("")
    elif len(subvalues) == 2:
        ax.set_ylabel(subvalues[0])
        ax.set_xlabel(subvalues[1])


def askColumnsToUser(correlogram_type, compound_and_ratios_df):
    """ "
    Asks user for columns to display and correlogram
    """
    col_name = CORRELOGRAM_TYPE_CONVERTER[correlogram_type]
    columns = getCorrelogramColumns(correlogram_type)
    answer = input(
        f"""Default {col_name} are: {','.join(columns)}. 
                   Press enter to confirm or write new column list in the same format"""
    )
    if answer:
        columns = answer.replace(" ", "").split(",") #.upper() was this the problem a problem but not THE problem
        # errors = filter(
        #     lambda i, col: col not in set(compound_and_ratios_df[col_name]),
        #     enumerate(columns), 
        # )
        errors = [col for col in columns if col not in set(compound_and_ratios_df[col_name])] #JJB this fixed TyepError: missing 1 required positional argument: 'col' when trying to actualy select cols

        while errors:
            columns = (
                inputEscape(
                    f"""Unknown {col_name}s {','.join(errors)}
                Please retype"""
                )
                .upper()
                .replace(" ", "")
                .split(",")
            )
            errors = filter(
                lambda i, col: col not in set(compound_and_ratios_df[col_name]),
                enumerate(columns),
            )
    return columns
