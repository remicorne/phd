from main import * 
from module.headtwitch import *
from module.getters import getTreatmentMapping

treatment_mapping = (
    {  # TO DO : change name to treatment_info and add columns in df REMI
        1: {"treatment": "vehicles",
            "color": "white",
            "experiments": ["dose_response", "agonist_antagonist"],
            },
        2: {"treatment": "0.3mg/kgTCB",
            "color": "firebrick",
            "experiments": ["dose_response"],
            },
        3: {"treatment": "3mg/kgTCB",
            "color": "red",
            "experiments": ["dose_response", "agonist_antagonist"],
            },
        4: {"treatment": "10mg/kgTCB",
            "color": "salmon",
            "experiments": ["dose_response"],
            },
        5: {"treatment": "0.3mg/kgMDL",
            "color": "black",
            "experiments": ["agonist_antagonist"],
            },
        6: {"treatment": "TCB+MDL",
            "color": "grey",
            "experiments": ["agonist_antagonist"],
            },
    }
)
region_subclassification = {
    'cortex': {'regions': ['OF', 'PL', 'CC', 'IC', 'M', 'SL1', 'SR1', 'AC', 'V'], 'color': 'mediumblue'},
    'subcortical_telencephalon': {'regions': ['Am', 'dH', 'vH', 'NAc', 'VM', 'DM', 'VL', 'DL'], 'color': 'orange'},
    'diencephalon': {'regions': ['MD', 'VPL', 'VPR', 'DG', 'Y'], 'color': 'darkorchid'},
    'mesencephalon': {'regions': ['SC', 'SN', 'VTA', 'DR', 'MR'], 'color': 'forestgreen'},
    'cerebellum': {'regions': ['CE'], 'color': 'peru'}}

filename = "TCB2_data_HPLC.csv"  # TCB2 #using current working directory plus file name
HT_filename = "TCB2_data_HT.csv"

experimental_info = {
    "dose_response": {"groups": [1, 2, 3, 4], "independant_vars": ["TCB2"]},
    "agonist_antagonist": {
        "groups": [1,3,5,6 ],  # JASMINETODO: the order here will be the one used by histogram and correlogram
        "independant_vars": [
            "TCB2",
            "MDL",
        ],
        "outliers": ["grubbs"],
        "correlation_statistics": ["pearson"],
        "quantitative_statistics": {
            "twoway_anova": True,
            "oneway_anova": True,
            "tukey": True,
        },
    },
}


saveMetadata(
    filename, treatment_mapping=treatment_mapping, experimental_info=experimental_info
)

saveMetadata(
    HT_filename, treatment_mapping=treatment_mapping, experimental_info=experimental_info
)

#26/9/23 running

# Exception has occurred: OptionError
# "No such keys(s): 'mode.use_inf_as_null'"
#   File "/Users/jasminebutler/Desktop/phd/debug.py", line 73, in <module>
#     sns.lineplot(data=rando_df, x='x',y='y')
# pandas._config.config.OptionError: "No such keys(s): 'mode.use_inf_as_null'"

# rando_df= pd.DataFrame({'x':[1,2,1,3], 'y': [3,3,3,5]})

# sns.lineplot(data=rando_df, x='x',y='y')

print('ja, genau')

#building 
# @get_or_add("quantitative_summary")
# def singleQuantitativeSummaryFig(
from module.constants import COLUMN_ORDER

filename
experiment="dose_response"
compound="5HIAA/5HT"
regions=COLUMN_ORDER['region'] #REMI i would ike this to work the same way it does for correlograms i.e. also specifying the order 
from_scratch=True

#slice aggstats df 
experimental_df = subselectDf(getAggregateStatsDf(filename), {"compound": compound, "experiment": experiment })   #REMI can use completly omnce has list comprehension
experimental_df = experimental_df[experimental_df['region'].isin(regions)]

#check that COLUMN_ORDER is consistent with region_subclasification i.e. CORTEX then SUBCORT...
region_order = regions #TODO check that regions inputted are sorted by subclasification
experimental_df = experimental_df.assign(region=lambda x: pd.Categorical(x["region"], categories=region_order, ordered=True)).sort_values("region")


# create a new column % of control mean/control_mean * 100
experimental_df.loc[:, "percentage_of_vehicles"] = experimental_df.groupby(
    "region"
)["mean"].transform(
    lambda x: (x / x.loc[experimental_df["treatment"] == "vehicles"].values[0])
    * 100
)

#melt df
plot_experimental_df = pd.melt(
    experimental_df,
    id_vars=["region", "treatment"],
    value_vars=["percentage_of_vehicles"],
)

#palete 
treatment_palette = {
        info["treatment"]: info["color"] for number, info in getTreatmentMapping(filename).items()
    }

#BUILD

fig, ax = plt.subplots(figsize=(12, 9))
sns.set_style("white")
sns.set_context("notebook")

# plot lines
sns.lineplot(data=plot_experimental_df,
    x="region",
    y="value",
    hue="treatment",
    palette=treatment_palette,
)

#MARKER SIZE mapping #seperate func

# add markers for each region
marker_mapping = {
    value["treatment"]: value["markers"]
    for value in treatment_mapping.values()
    if experiment in value["experiments"]
}
# Define the minimum and maximum marker sizes
min_marker_size = 20
max_marker_size = 100

# Calculate the marker sizes based on the difference from 100
plot_experimental_df["marker_size"] = abs(plot_experimental_df["value"] - 100)

# Normalize marker sizes between the minimum and maximum sizes
plot_experimental_df["marker_size"] = (
    (
        plot_experimental_df["marker_size"]
        - plot_experimental_df["marker_size"].min()
    )
    / (
        plot_experimental_df["marker_size"].max()
        - plot_experimental_df["marker_size"].min()
    )
) * (max_marker_size - min_marker_size) + min_marker_size

# Plot with adjusted marker sizes
sns.scatterplot(
    data=plot_experimental_df,
    x="region",
    y="value",
    hue="treatment",
    palette=treatment_palette,
    style="treatment",
    markers=marker_mapping,
    legend=False,
    size="marker_size",
    sizes=(min_marker_size, max_marker_size),  # Set the range of marker sizes
    alpha=0.7,  # Adjust the marker transparency if desired
)

# Add AXSPAN to plot another funct ?

region_subclassification = getRegionSubclassification(filename)
region_subclassification = {
    group.replace("_", " "): properties
    for group, properties in region_subclassification.items()
}
for group, properties in region_subclassification.items():
    what_regions = properties["regions"]
    color = properties["color"]
    label = group.replace("_", " ")
    if regions:
        start_idx = region_order .index(regions[0]) - 0.5
        end_idx = region_order .index(regions[-1]) + 0.5
        ax.axvspan(start_idx, end_idx, facecolor=color, alpha=0.2, label=label)
        # Adjust x-axis limits
        ax.set_xlim(-0.5, len(region_order ) - 0.5)  # Set the x-axis limits

        # Add axvspan label

        label_x = (start_idx + end_idx) / 2
        label_y = ax.get_ylim()[1] - 0.05 * (
            ax.get_ylim()[1] - ax.get_ylim()[0]
        )  # Adjust the label_y coordinate
        words = label.split()
        lines = [
            words[i : i + 2] for i in range(0, len(words), 2)
        ]  # Split words into pairs for multiline display
        line_height = 18  # Adjust the height between lines

        for i, line in enumerate(lines):
            line_label = "\n".join(line)  # Join words with a line break
            current_label_y = label_y - i * line_height
            ax.annotate(
                line_label,
                xy=(label_x, current_label_y),
                xycoords="data",
                xytext=(0, -12),
                textcoords="offset points",
                ha="center",
                va="top",
                color=color,
                fontsize=18,
            )





print('kitty cat')








