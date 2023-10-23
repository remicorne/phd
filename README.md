# phd

Notebook is used to do the actual data processing
Module contains all necessary functions

getter (ie all functions strating with 'get') is used to get data. The getter will either generate the data or retrieve it from the cache if it exists

the cache

REMI: QUESTIONS FOR JASMINE JAS: i think this convention is ok - discuss

- I'm wondering what to do with figure naming because potentially many figure for the grouping will we be built. Im currently using an automatic naming convention
  {"histogram": 'f"{experiment}_for_{compound}_in_{region}"',
  "correlogram": 'f"{experiment}_{correlogram_type}_{buildCorrelogramFilenmae(to*correlate, columns)}"',
  "head_twitch_histogram": 'f"head_twitch_histogram*{experiment}_for_{to_plot}"',}
  but maybe the user would want to pick the name themself? probably better for them to remember what is what in the case of multiple stats choices

##### TO USE:

add csv with columns : mouse_id , group_id , COMPOUND_REGION... or BEHAVIOR_TIME (e.g. HT_20) to input folder

fill info in cell 1 of notebook (compound_ratio_mapping, ect)

perform outlier selection for experiment (including ratios chiosen in first cell)

generate quantitative histograms and aggregated stats table functions : REMI?

generate correlograms (use case for all three in functions) : REMI?
clasical_corellogram : getAndPlotSingleCorrelogram(filename, experiment='agonist_antagonist', correlogram_type='compound',  
 to_correlate='GLU', p_value_threshold=0.05, n_minimum=5, from_scratch= True)

    square_correlogram       :      getAndPlotSingleCorrelogram(filename, experiment='agonist_antagonist', correlogram_type='compound',
                                                                to_correlate='GLU-GABA', p_value_threshold=0.05, n_minimum=5, from_scratch= True)


    bar_corellogram
                                                    #see whatsapp image 3/5/23
        within BR       /       within compound
