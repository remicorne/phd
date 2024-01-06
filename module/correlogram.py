from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
from module.getters import (
    getExperimentalInfo,
    getTreatmentMapping,
    getCompoundAndRatiosDf,
)
from module.utils import (
    askMultipleChoice,
    flatten,
    get_or_add,
    inputEscape,
    plotExperiment,
)
from module.statistics import (
    CORR_STAT_METHODS,
    isSignificant,
    getPearsonR,
    getPearsonPValue,
)  ##FIX ME REMIs OLD REDUNDANT SYSTEM lol just without diff corr methids
from module.constants import getCorrelogramColumns, CORRELOGRAM_TYPE_CONVERTER
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from module.Matrix import Matrix

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
    corr_method="pearson"  # HARD CODE FIXME
    # hierarchical_clustering=None,
):
    """
    This is the function that is called by the user, it will call buildExperimentalCorrelogram that builds a signle correlogram
    The function can be called with parameters, user input will be required if not
    """
    (
        filename,
        experiment,
        correlogram_type,
        to_correlate,
        p_value_threshold,
        n_minimum,
        columns,
        corr_method,
        from_scratch,
    ) = corrSelector(  # generic prompter for selecting corr matrices, probably need to add pearson/spearman
        filename,
        experiment=experiment,
        correlogram_type=correlogram_type,
        to_correlate=to_correlate,
        p_value_threshold=p_value_threshold,
        n_minimum=n_minimum,
        columns=columns,
        corr_method=corr_method,
        from_scratch=from_scratch,
    )
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate

    buildExperimentalCorrelogram(
        filename,
        experiment=experiment,
        correlogram_type=correlogram_type,
        to_correlate=to_correlate,
        p_value_threshold=p_value_threshold,
        n_minimum=n_minimum,
        columns=columns,
        corr_method=corr_method,
        from_scratch=from_scratch,
        # hierarchical_clustering=hierarchical_clustering,
    )


@get_or_add("correlogram")
def buildExperimentalCorrelogram(
    filename,
    experiment,
    correlogram_type,
    to_correlate,
    p_value_threshold,
    n_minimum,
    columns,
    from_scratch,  # Used in decorator
    corr_method,
):
    matrices = buildExperimentCorrmatrices(
        filename,
        experiment,
        correlogram_type,  # compound / ratio
        to_correlate,  # whihc compound/s / ratio/s
        p_value_threshold,
        n_minimum,
        columns,  # ordered to plot
        corr_method,  # 'pearson' 'spearman' 'kendall'
        from_scratch,  # Used in decorator
    )

    fig = plotExperiment(matrices, plotCorrelogram)

    return fig


# prompts for correlation matrices
# done by experiment
def corrSelector(  # generic prompter for selecting corr matrices
    filename,
    experiment=None,
    correlogram_type=None,
    to_correlate=None,
    p_value_threshold=0.05,
    n_minimum=5,
    columns=None,
    corr_method=None,
    from_scratch=None,
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate
):
    """
    This is the function that is called by the user, it will generate the information sto be passed to buildExperimentCorrmatrices
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
        columns
        if columns
        else askColumnsToUser(
            correlogram_type, compound_and_ratios_df
        )  # REMI: for me I would like to be able to manuely code these in pls ass the if not assigned, danke
    )
    corr_method = (
        corr_method
        if corr_method
        else askMultipleChoice(
            "Which correlation method?",
            CORR_STAT_METHODS.keys(),
        )
    )
    return (
        filename,
        experiment,
        correlogram_type,
        to_correlate,
        p_value_threshold,
        n_minimum,
        columns,
        corr_method,
        from_scratch,
    )


# for an experiment and given set of correlations
# we return a list of lists called matrices:
# [[df_to_corr , correlation_matrix , T_F_mask_matrix , treatment , to_correlate],
# [ other treatment ],
# [ other treatment ],
# [ other treatment ]]
def buildExperimentCorrmatrices(
    filename,
    experiment,
    correlogram_type,  # compound / ratio
    to_correlate,  # whihc compound/s / ratio/s
    p_value_threshold,
    n_minimum,
    columns,  # ordered to plot
    corr_method,  # 'pearson' 'spearman' 'kendall'
    from_scratch,  # Used in decorator
):
    # get and subselect df in long format
    compound_and_ratios_df = getCompoundAndRatiosDf(filename).sort_values("group_id")
    experiment_df = compound_and_ratios_df[compound_and_ratios_df.experiment == experiment]
    matrices = []
    conserved_rows = []
    conserved_cols = []
    for treatment, group_df in experiment_df.groupby(
        by=["treatment"], sort=False  # preserve order in subselection_df_ordered
    ):
        matrix = Matrix(
            data=group_df,
            experiment=experiment,
            treatment=treatment[0],
            pvalue_threshold=p_value_threshold,
            n_minimum=n_minimum,
            between=correlogram_type,
            variables=to_correlate.split('-'),
            method=corr_method,
            accross="compound" if correlogram_type == "region" else "region",
            columns=columns,
        )
        conserved_rows = (
            set(matrix.corr.index)
            if not conserved_rows
            else conserved_rows.intersection(matrix.corr.index)
        )
        conserved_cols = (
            set(matrix.corr.columns)
            if not conserved_cols
            else conserved_cols.intersection(matrix.corr.columns)
        )
        matrices.append(matrix)

    # Homogenize all datasets
    for matrix in matrices:
        rows_to_drop = [row for row in matrix.corr.index if row not in conserved_rows]
        cols_to_drop = [col for col in matrix.corr.columns if col not in conserved_cols]
        matrix.corr = matrix.corr.drop(index=rows_to_drop, columns=cols_to_drop)
    return matrices


# I would think that like build correlogram there may be build hirichal correlogram whihc would use this for now its all under construction
def perform_hierarchical_clustering(correlograms):
    """_summary_

    Args:
        correlograms (_type_): _description_

    Returns:
        _type_: _description_
    """
    hierarchical_correlograms = []
    for ordered_corr, p_corr, treatment, tocolorlate in correlograms:
        print(f"Hierarchical clustering for {treatment} correlating {tocolorlate}")
        # Linkage matrix
        Z = linkage(squareform(1 - abs(ordered_corr)), "complete")

        # Plot dendrogram based on col of ordered_corr
        hierarchical_labels = dendrogram(
            Z, labels=ordered_corr.columns, orientation="top", leaf_rotation=90
        )["ivl"]
        plt.title(f"{treatment} correlating {tocolorlate}", pad=20)
        plt.tight_layout()
        plt.show()

        # not used - check inout prams
        threshold = 0.8  # Adjust as needed
        labels = fcluster(Z, threshold, criterion="distance")

        # Create mapping from original labels to sorted indices
        label_indices_map = {label: i for i, label in enumerate(hierarchical_labels)}
        # Reordering columns based on clustering using label_indices_map
        hierarchical_labels_order = [
            label_indices_map[label] for label in ordered_corr.columns
        ]
        hierarchical_corr = ordered_corr.iloc[
            hierarchical_labels_order, hierarchical_labels_order
        ]
        hierarchical_p_corr = p_corr.iloc[
            hierarchical_labels_order, hierarchical_labels_order
        ]

        hierarchical_correlograms.append(
            [hierarchical_corr, hierarchical_p_corr, treatment, tocolorlate]
        )
    return hierarchical_correlograms


def plotCorrelogram(matrix, ax):
    """
    Correlogram plotter for single correlation matrix ~ to be fed to plotExperiment()
    input: a single element from matricies i.e. for one treatment
    output:  ax with graph plotted
    """

    ax.set_title(
        matrix.get_title(), fontsize=28, pad=20, y=1
    )  # Adjust the y position of the title manually for square correlogram

    sns.heatmap(
        matrix.corr,
        vmin=-1,
        vmax=1,
        square=True,
        annot=True,
        cmap="coolwarm",
        annot_kws={"size": 8},
        ax=ax,
        cbar_kws={"shrink": 0.7},  # adj color bar size
    )
    ax.set_xticklabels(
        ax.get_xticklabels()
    )  # rotation=45, horizontalalignment='right',

    ax.set_ylabel(matrix.var1, fontsize=28)
    ax.set_xlabel(matrix.var2, fontsize=28)


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
        columns = answer.replace(" ", "").split(
            ","
        )  # .upper() was this the problem a problem but not THE problem
        # errors = filter(
        #     lambda i, col: col not in set(compound_and_ratios_df[col_name]),
        #     enumerate(columns),
        # )
        errors = [
            col for col in columns if col not in set(compound_and_ratios_df[col_name])
        ]  # JJB this fixed TyepError: missing 1 required positional argument: 'col' when trying to actualy select cols

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
