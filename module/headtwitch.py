# from matplotlib import pyplot as plt
# import pandas as pd
# from module.quantitative import processQuantitativeStats
# from module.getters import getHeadTwitchDf, getExperimentalInfo, updateQuantitativeStats
# from module.histogram import buildHeadTwitchHistogramData, buildHistogram
# from module.metadata import applyTreatmentMapping
# from module.utils import figure_cache
# import seaborn as sns


# ## TODO:REMI stats logic, outliers, prompts


# @figure_cache("head_twitch_histogram")
# def headTwitchHistogram(
#     HT_filename,
#     experiment=None,
#     vairable=None,
#     outlier_test=None,
#     p_value_threshold=None,
#     from_scratch=None,
# ):
#     HT_data = getHeadTwitchDf(HT_filename)
#     data, order, palette = buildHeadTwitchHistogramData(
#         HT_filename, experiment, vairable
#     )

#     title = vairable.replace("_", " at ") + " min"
#     ylabel = "events / min"

#     # REMI add stats calcultaion to input to histogram builder and outlier same as for quantativeHistograms()
#     # the last quantitative test is coded to return the labels directly, thus the need for the bool
#     (is_significant, significance_infos, test_results) = processQuantitativeStats(
#         getExperimentalInfo(HT_filename)[experiment], data, p_value_threshold
#     )

#     updateQuantitativeStats(
#         HT_filename,
#         [
#             {
#                 "data_type": "HT",
#                 "experiment": experiment,
#                 "compound": None,
#                 "region": None,
#                 **test_result,
#             }
#             for test_result in test_results
#         ],
#     )

#     fig = buildHistogram(
#         title,
#         ylabel,
#         data,
#         order,
#         hue='treatment',
#         palette=palette,
#         significance_infos=significance_infos if is_significant else None,
#     )
#     return fig
