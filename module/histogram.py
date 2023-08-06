from module.imports import *

# TODO: Handle columns and more generality


# TODO: Plug stat methods in here
@get_or_add("histogram")
def histogram(
    filename,
    experiment,
    compound,
    region,
    p_value_threshold,
    from_scratch=False,  # Used in decorator
):
    data, order, palette = buildHistogramData(
        filename, experiment, compound, region, p_value_threshold
    )
    fig = buildHistogram(
        compound,
        region,
        data,
        order,
        palette,
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
