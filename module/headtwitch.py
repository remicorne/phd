from matplotlib import pyplot as plt
import pandas as pd
from module.quantitative import processQuantitativeStats
from module.getters import getHeadTwitchDf, getExperimentalInfo
from module.histogram import buildHeadTwitchHistogramData, buildHistogram
from module.metadata import applyTreatmentMapping
from module.utils import get_or_add
import seaborn as sns


## TODO:REMI stats logic, outliers, prompts


@get_or_add("head_twitch_histogram")
def headTwitchHistogram(
    HT_filename,
    experiment=None,
    vairable=None,
    outlier_test=None,
    p_value_threshold=None,
    from_scratch=None,
):
    HT_data=getHeadTwitchDf(HT_filename)
    data, order, palette = buildHeadTwitchHistogramData(HT_filename, experiment, vairable)
    
    title = vairable.replace('_', ' at ') + ' min'
    ylabel = 'events / min'

    #REMI add stats calcultaion to input to histogram builder and outlier same as for quantativeHistograms()
    # the last quantitative test is coded to return the labels directly, thus the need for the bool
    (is_significant,
    significance_infos,
    ) = processQuantitativeStats(getExperimentalInfo(HT_filename)[experiment], data, p_value_threshold)

    fig = buildHistogram(title, 
                         ylabel, 
                         data, 
                         order, 
                         palette, 
                         hue=None, 
                         significance_infos=significance_infos if is_significant else None
                         )
    return fig

