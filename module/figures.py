from matplotlib import pyplot as plt
import numpy as np
from module.utils import *
from module.getters import *
from module.statistics import *
from module.metadata import *
from module.getters import *
import seaborn as sns
from statannotations.Annotator import Annotator
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
                    Use {correlogram_type} or {correlogram_type}/q{correlogram_type} for simple correlogram, or {correlogram_type}-{correlogram_type} or {correlogram_type}/{correlogram_type}-{correlogram_type} for square correlogram
                    Possibilities: {COLUMN_ORDER[correlogram_type]}
                    """
        )
    )
    columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
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
    columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
    pivot_columns = {
        "region": ["region", "compound"],
        "compound": ["compound", "region"],
    }[correlogram_type]
    # Because table is îvoted on region and compound even in the case of simple correlogram I have to duplicate the selector in that case to avoid ugly labelling
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
        columns = columns if columns else CORRELOGRAM_COLUMN_ORDER[correlogram_type]
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
    fig, axs = plt.subplots(2, 2, figsize=(22, 22))
    axs = flatten(axs)  # Put all the axes into a list at the same level
    for (correlogram_df, p_value_mask, treatment, subvalues), ax in zip(
        correlograms, axs
    ):
        plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax)
    fig.tight_layout()
    return fig


def plotCorrelogram(correlogram_df, p_value_mask, treatment, subvalues, ax):
    if np.array_equal(correlogram_df, correlogram_df.T):  # remove duplicate data
        mask = np.triu(np.ones(p_value_mask.shape, dtype=bool), k=1)
        p_value_mask[mask] = True
        np.fill_diagonal(
            p_value_mask.values, False
        )  # this maked the diagonal correlations of 1 visible

    heatmap = sns.heatmap(
        correlogram_df,
        vmin=-1,
        vmax=1,
        square=True,
        annot=True,
        cmap="BrBG",
        mask=p_value_mask,
        annot_kws={"size": 6},
        ax=ax,
    )
    ax.set_xticklabels(
        ax.get_xticklabels()
    )  # rotation=45, horizontalalignment='right',
    title = ax.set_title(
        f"{'-'.join(subvalues)} in {treatment}", fontsize=28, pad=20, y=0.9
    )  # Adjust the y position of the title manually

    if len(subvalues) == 1:
        ax.set_ylabel("")
        ax.set_xlabel("")
    elif len(subvalues) == 2:
        ax.set_ylabel(subvalues[0])
        ax.set_xlabel(subvalues[1])


# TODO: Plug stat methods in here
@get_or_add("histogram")
def histogram(
    filename,
    experiment,
    compound,
    region,
    p_value_threshold,
    from_scratch=False,  # Used in decorator
    select_outliers=False,
):
    compound_and_ratios_df = getCompoundAndRatiosDf(
        filename
    )  # this is not the full ratios df, its only intra region compound ratios for nom
    subselection_df = compound_and_ratios_df[
        (compound_and_ratios_df.experiment == experiment)
        & (compound_and_ratios_df.compound == compound)
        & (compound_and_ratios_df.region == region)
    ]
    subselection_df = subselection_df[["value", "mouse_id", "treatment"]]
    treatment_mapping = getTreatmentMapping(filename)
    experimental_info = getExperimentalInfo(filename)[experiment]
    subselection_df = applyHistogramOutlierMapping(
        pd.concat(
            [
                subselection_df,
                # pd.DataFrame(
                #     {"treatment": ["vehicles"], "value": [50], "mouse_id": [203]}
                # ),
            ]
        ),
        experimental_info["outliers"],
        p_value_threshold,
    )
    order = [
        treatment_mapping[group]["treatment"]
        for group in treatment_mapping
        if experiment in treatment_mapping[group]["experiments"]
    ]
    palette = {info["treatment"]: info["color"] for info in treatment_mapping.values()}

    if select_outliers:
        hue_or_palette = {"hue": "is_outlier"}
        fig = buildHistogram(
            compound, region, subselection_df, order, palette, hue_or_palette
        )
        plt.show(fig)
        remove_outliers = input("Remove outliers?") == "y"
        if remove_outliers:
            subselection_df = subselection_df[subselection_df.is_outlier == False]

    hue_or_palette = {"palette": palette}
    fig = buildHistogram(
        compound, region, subselection_df, order, palette, hue_or_palette
    )
    return fig


def buildHistogram(title, ylabel, data, order, palette, hue=None):
    fig, ax = plt.subplots(figsize=(20, 10))
    ax = sns.barplot(
        x="treatment",
        y="value",
        data=data,
        palette=palette,
        ci=68,
        order=order,
        capsize=0.1,
        alpha=0.8,
        errcolor=".2",
        edgecolor=".2",
    )
    #
    hue, palette = list(hue.items())[0] if hue else (None, palette)
    ax = sns.swarmplot(
        x="treatment",
        y="value",
        hue=hue,
        palette=palette,
        order=order,
        data=data,
        edgecolor="k",
        linewidth=1,
        linestyle="-",
    )
    ax.tick_params(labelsize=24)
    ax.set_ylabel("ng/mg of tissue", fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_xlabel(" ", fontsize=20)  # treatments
    ax.set_title(title, y=1.04, fontsize=34)  # '+/- 68%CI'
    sns.despine(left=False)
    return fig


def put_significnce_stars(
    stat_data,
    ax,
    treatment_dict,
    test_path,
    data=None,
    x=None,
    y=None,
    order=None,
    sheet="5HT_DL",
):  # , p_values):
    if len(df_significant.index) > 0:
        print(df_significant)
        p_values = df_significant["p-adj"].values
        pairs = [
            (treatment_dict[i[1]["group1"]], treatment_dict[i[1]["group2"]])
            for i in df_significant.iterrows()
        ]

        annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=order)
        # https://stackoverflow.com/questions/64081570/matplotlib-marker-annotation-fontsize-not-shrinking-below-1pt-in-pdf
        annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
        annotator.set_pvalues_and_annotate(p_values)

    return ax


@get_or_add("head_twitch_histogram")
def headTwitchHistogram(
    filename,
    HT_filename,
    experiment="agonist_antagonist",
    p_value_threshold=0.05,
    to_plot=["HT_20"],
):
    HT_df = getRawHeadTwitchDf(HT_filename)
    applyTreatmentMapping(HT_df, filename)
    treatment_mapping = getTreatmentMapping(filename)
    treatment_palette = {
        info["treatment"]: info["color"] for number, info in treatment_mapping.items()
    }
    treatments = [
        treatment_mapping[group]["treatment"]
        for group in treatment_mapping
        if experiment in treatment_mapping[group]["experiments"]
    ]
    # treatments = [treatment_mapping[str(group)]['treatment'] for group in getExperimentalInfo(filename)[experiment]['groups']]
    experimental_df = HT_df[
        HT_df["treatment"].isin(treatments)
    ]  # select only relevent treatments
    columns = to_plot + ["treatment"]  # slice df for plotting
    experimental_df = experimental_df[columns]
    experimental_df = pd.melt(
        experimental_df, id_vars=["treatment"], value_vars=to_plot
    )

    fig, ax = plt.subplots(figsize=(20, 10))
    if len(to_plot) == 1:
        time = to_plot[0].split("_")[1]
        sns.barplot(
            data=experimental_df,
            x="treatment",
            y="value",
            ci=68,
            order=treatments,
            capsize=0.1,
            alpha=0.8,
            palette=treatment_palette,
            errcolor=".2",
            edgecolor=".2",
        )
        sns.swarmplot(
            data=experimental_df,
            x="treatment",
            y="value",
            order=treatments,
            palette=treatment_palette,
            edgecolor="k",
            linewidth=1,
            linestyle="-",
        )
        ax.set_title(f"Head Twitch at {time} minutes", y=1.04, fontsize=34)
    else:
        sns.barplot(
            data=experimental_df,
            x="treatment",
            y="value",
            hue="variable",
            ci=68,
            order=treatments,
            capsize=0.1,
            alpha=0.8,
            errcolor=".2",
            edgecolor=".2",
        )
        sns.swarmplot(
            data=experimental_df,
            x="treatment",
            y="value",
            hue="variable",
            order=treatments,
            edgecolor="k",
            linewidth=1,
            linestyle="-",
            dodge=True,
            marker="o",
        )
        ax.set_title(f"Head Twitch", y=1.04, fontsize=34)

    ax.set_ylabel("twitches / min", fontsize=24)
    ax.set_xlabel(" ", fontsize=24)
    ax.tick_params(labelsize=24)
    sns.despine(left=False)
    return fig