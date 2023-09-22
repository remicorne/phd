from matplotlib import pyplot as plt
import pandas as pd
from module.getters import getRawHeadTwitchDf
from module.getters import getTreatmentMapping
from module.metadata import applyTreatmentMapping
from module.utils import get_or_add
import seaborn as sns


## TODO: modularize with generic histogram funcs


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
