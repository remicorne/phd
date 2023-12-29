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
from module.statistics import CORR_STAT_METHODS, isSignificant, getPearsonR, getPearsonPValue ##FIX ME REMIs OLD REDUNDANT SYSTEM lol just without diff corr methids
from module.constants import getCorrelogramColumns, CORRELOGRAM_TYPE_CONVERTER
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

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
    corr_method = 'pearson' #HARD CODE FIXME
    # hierarchical_clustering=None,
):
    """
    This is the function that is called by the user, it will call buildSingleCorrelogram that builds a signle correlogram
    The function can be called with parameters, user input will be required if not
    """
    filename,experiment,correlogram_type,to_correlate,p_value_threshold,n_minimum,columns,corr_method,from_scratch=corrSelector( # generic prompter for selecting corr matricies, probably need to add pearson/spearman
                                                                                                                                filename,
                                                                                                                                experiment=experiment,
                                                                                                                                correlogram_type=correlogram_type,
                                                                                                                                to_correlate=to_correlate,
                                                                                                                                p_value_threshold=p_value_threshold,
                                                                                                                                n_minimum=n_minimum,
                                                                                                                                columns=columns,
                                                                                                                                corr_method = corr_method,
                                                                                                                                from_scratch=from_scratch )
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate 

    buildSingleCorrelogram(
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
def buildSingleCorrelogram(
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
    matricies=buildExperimentCorrMatricies(
                                            filename,
                                            experiment,
                                            correlogram_type, #compound / ratio
                                            to_correlate, # whihc compound/s / ratio/s
                                            p_value_threshold,
                                            n_minimum,
                                            columns, #ordered to plot
                                            corr_method, # 'pearson' 'spearman' 'kendall'
                                            from_scratch,  # Used in decorator
                                            )

    fig = plotCorrelograms(matricies)

    return fig


#prompts for correlation matricies
#done by experiment
def corrSelector( # generic prompter for selecting corr matricies
    filename,
    experiment=None,
    correlogram_type=None,
    to_correlate=None,
    p_value_threshold=0.05,
    n_minimum=5,
    columns=None,
    corr_method = None,
    from_scratch=None,
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate 
):
    """
    This is the function that is called by the user, it will generate the information sto be passed to buildExperimentCorrMatricies 
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
    corr_method = (
        corr_method
        if corr_method
        else askMultipleChoice(
            "Which correlation method?",
            CORR_STAT_METHODS.keys(), )
        )
    return filename,experiment,correlogram_type,to_correlate,p_value_threshold,n_minimum,columns,corr_method,from_scratch
        

#for an experiment and given set of correlations 
# we return a list of lists called matricies: 
# [[df_to_corr , correlation_matrix , T_F_mask_matrix , treatment , correlated],
# [ other treatment ],
# [ other treatment ],
# [ other treatment ]] 
def buildExperimentCorrMatricies(
    filename,
    experiment,
    correlogram_type, #compound / ratio
    to_correlate, # whihc compound/s / ratio/s
    p_value_threshold,
    n_minimum,
    columns, #ordered to plot
    corr_method, # 'pearson' 'spearman' 'kendall'
    from_scratch,  # Used in decorator
):

    #get and subselect df in long format
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    value_type = {"compound": "region", "region": "compound"}[correlogram_type]
    to_correlate = to_correlate.split("-")
    subselection_df = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df[value_type].isin(columns))
    ].query("|".join([f"{correlogram_type}=='{value}'" for value in to_correlate]))

    #columns to pivot the df - if correlogram_type is compound (i.e. GLU-GLU) then pivot_columns ["compound", "region"]
    pivot_columns = {
        "region": ["region", "compound"],
        "compound": ["compound", "region"],
    }[correlogram_type]

    # Because table is înveted on region and compound even in the case of simple correlogram I have to duplicate the selector in that case to avoid ugly labelling
    pivot_column_value = (
        to_correlate if len(to_correlate) == 2 else to_correlate + to_correlate
    )
    treatment_mapping = getTreatmentMapping(filename)

    #order treatments within experiment
    order = [
        treatment_mapping[str(group)]["treatment"]
        for group in getExperimentalInfo(filename)[experiment]["groups"]
    ]
    subselection_df_ordered = subselection_df.iloc[
        subselection_df["treatment"]
        .astype(pd.CategoricalDtype(categories=order))
        .argsort()
    ]  

    matricies = []
    for treatment, group_df in subselection_df_ordered.groupby(
        by=["treatment"], sort=False #preserve order in subselection_df_ordered
    ):
        # Here we just pivot our structure to one where each column is a region or compound and each line a mouse    
        pivot_df = group_df.pivot_table(
            values="value", index=group_df["mouse_id"], columns=pivot_columns
        )
        # This is the list of methods to pass to df.corr. I know you used to pass 'pearson' as a string but this way the R calculation and pvalue mask are 'linked' by being done on the same line
        
        print(f"correlating using {corr_method}")
        base_scipy_test, get_corr, get_p_val = CORR_STAT_METHODS[corr_method]#fetchnig functions for stats 
        methods =[get_corr, isSignificant(get_p_val, p_value_threshold)]

        # order columns in desired plotting order
        columns = columns 
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

        matricies.append(
            [pivot_df_ordered, correlogram_df, p_value_mask.astype(bool), treatment[0], to_correlate]
            #df_to_corr , correlation_matrix , T_F_mask_matrix , treatment , correlated
        )

    return matricies 

# I would think that like build correlogram there may be build hirichal correlogram whihc would use this for now its all under construction 
def perform_hierarchical_clustering(correlograms):
    hierarchical_correlograms = []
    for ordered_corr, p_corr, treatment, tocolorlate in correlograms:
        print(f"Hierarchical clustering for {treatment} correlating {tocolorlate}")
        # Linkage matrix
        Z = linkage(squareform(1 - abs(ordered_corr)), 'complete')
        
        # Plot dendrogram based on col of ordered_corr
        hierarchical_labels = dendrogram(Z, labels=ordered_corr.columns, orientation='top', leaf_rotation=90)['ivl']
        plt.title(f"{treatment} correlating {tocolorlate}", pad=20) 
        plt.tight_layout()
        plt.show()

        #not used - check inout prams 
        threshold = 0.8  # Adjust as needed
        labels = fcluster(Z, threshold, criterion='distance') 
        
        # Create mapping from original labels to sorted indices
        label_indices_map = {label: i for i, label in enumerate(hierarchical_labels)}
        # Reordering columns based on clustering using label_indices_map
        hierarchical_labels_order = [label_indices_map[label] for label in ordered_corr.columns]
        hierarchical_corr = ordered_corr.iloc[hierarchical_labels_order, hierarchical_labels_order]
        hierarchical_p_corr = p_corr.iloc[hierarchical_labels_order, hierarchical_labels_order]
        
        hierarchical_correlograms.append([hierarchical_corr, hierarchical_p_corr, treatment, tocolorlate])
    return hierarchical_correlograms


########## replaced 
def plotCorrelograms(matricies):
    #JJB might be nice to have one color bar for all figures
    fig, axs = plt.subplots(2, 2, figsize=(22, 22))
    axs = flatten(axs)  # Put all the axes into a list at the same level
    for (df_to_corr,correlogram_df, p_value_mask, treatment, correlated), ax in zip(
        matricies, axs
    ):
        plotCorrelogram(correlogram_df, p_value_mask, treatment, correlated, ax)
    fig.tight_layout()
    return fig

def plotCorrelogram(correlogram_df, p_value_mask, treatment, correlated, ax):
    
    if np.array_equal(correlogram_df, correlogram_df.T):  # TRIANGLE CORRELOGRAMS remove duplicate data  
        title = ax.set_title(f"{'-'.join(correlated)} in {treatment}", fontsize=28, pad=20, y=0.9)  # Adjust the y position of the title manually #JJB set

        mask = np.triu(np.ones(p_value_mask.shape, dtype=bool), k=1)
        p_value_mask[mask] = True
        np.fill_diagonal(
            p_value_mask.values, False
        )  # this makes the diagonal correlations of 1 visible
    else:
        title = ax.set_title(f"{'-'.join(correlated)} in {treatment}", fontsize=28, pad=20, y=1)  # Adjust the y position of the title manually for square correlogram

    heatmap = sns.heatmap(
        correlogram_df,
        vmin=-1,
        vmax=1,
        square=True,
        annot=True,
        cmap='coolwarm', #matplotlib.colormaps['coolwarm'], #  Spectral_r
        mask=p_value_mask,
        annot_kws={"size": 8}, #, 'fontweight':'bold'
        ax=ax,
        # fmt=".2f",
        cbar_kws={"shrink": 0.8} #adj color bar size
    )
    ax.set_xticklabels(
        ax.get_xticklabels()
    )  # rotation=45, horizontalalignment='right',

    if len(correlated) == 1:
        # ax.set_ylabel("")
        # ax.set_xlabel("")
        ax.set_ylabel(correlated[0], fontsize=28)
        ax.set_xlabel(correlated[0], fontsize=28)
    elif len(correlated) == 2:
        ax.set_ylabel(correlated[0], fontsize=28)
        ax.set_xlabel(correlated[1], fontsize=28)


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
