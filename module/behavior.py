### OBSERVED BEHAVIOR ### 

#HT : head twitch 


######### HEAD TWITCH BEHAVIOR FUNCS ###########
def getRawHTDf(filename):
    return getOrBuildDf(filename, "raw_HT_df", buildRawHTDf)


def buildRawHTDf(filename):
    file_name, file_type = filename.split(".")
    if not file_type == "csv":
        raise Exception(f'METHOD TO DEAL WITH FILE TYPE "{file_type}" ABSENT')
    if not os.path.isfile(f"{INPUT_DIR}/{filename}"):
        raise Exception(f'FILE {filename} IS ABSENT IN "input/" DIRECTORY')
    return pd.read_csv(f"{INPUT_DIR}/{filename}", header=0).replace(
        np.nan, 0
    )  # to set all 0 to Nan


def buildHTHistogram(
    filename,
    HT_filename,
    experiment="agonist_antagonist",
    p_value_threshold=0.05,
    to_plot=["HT_20"],
):
    HT_df = getRawHTDf(HT_filename)
    applyTreatmentMapping(HT_df, filename)
    treatment_mapping = getTreatmentMapping(filename)
    treatment_palette = {info['treatment']:info['color'] for number, info in treatment_mapping.items()}
    treatments = [treatment_mapping[str(group)]['treatment'] for group in getExperimentalInfo(filename)[experiment]['groups']]
    experimental_df = HT_df[HT_df['treatment'].isin(treatments)] #select only relevent treatments
    columns = to_plot +['treatment'] # slice df for plotting
    experimental_df = experimental_df[columns]
    experimental_df = pd.melt(
        experimental_df, id_vars=["treatment"], value_vars=to_plot
    )

    fig, ax = plt.subplots(figsize=(20, 10))
    if len(to_plot)==1:
        sns.barplot(data = experimental_df, x = 'treatment', y='value', 
                    ci=68, order=treatments,capsize=.1, alpha=0.8, palette=treatment_palette,
                    errcolor=".2", edgecolor=".2")
    else:
        sns.barplot(data = experimental_df, x = 'treatment', y='value', hue='variable', 
                    ci=68, order=treatments, capsize=.1, alpha=0.8, errcolor=".2", edgecolor=".2")

    ax.set_title(f'Head Twitch', y=1.04, fontsize=34)
    ax.set_ylabel("twitches / min",fontsize=24)
    ax.set_xlabel(" ",fontsize=24)
    ax.tick_params(labelsize=24)
    sns.despine(left=False)
    return fig


def getHTHistogram(
    filename,
    HT_filename,
    experiment="dose_responce",
    p_value_threshold=0.05,
    to_plot=[],
    from_scratch=None,
):
    identifier = f"HT_Histogram_{experiment}_for_{to_plot}"
    from_scratch = (
        from_scratch
        if from_scratch is not None
        else input("Recalculate figure even if previous version exists? (y/n)") == "y"
    )
    if from_scratch or not isCached(filename, identifier):
        fig = buildHTHistogram(
            filename,
            HT_filename,
            experiment=experiment,
            p_value_threshold=p_value_threshold,
            to_plot=to_plot,
        )

        cache(filename, identifier, fig)
        saveHTHistogram(fig, identifier)
    else:
        fig = getCache(filename, identifier)
    fig.show()
