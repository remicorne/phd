# phd

Notebook is used to do the actual data processing
Module contains all necessary functions

getter (ie all functions strating with 'get') is used to get data. The getter will either generate the data or retrieve it from the cache if it exists

the cache

QUESTIONS FOR JASMINE

- Finish filling up experimental info for dose response and add statistical test/outlier test you might think about
- Also wouldnt the quantitative tests and outlier test be generic and not specific to experiments?
- add region subclassification
- I'm wondering what to do with figure naming because potentially many figure for the grouping will we be built. Im currently using an automatic naming convention
  {"histogram": 'f"{experiment}_for_{compound}_in_{region}"',
  "correlogram": 'f"{experiment}_{correlogram_type}_{buildCorrelogramFilenmae(to*correlate, columns)}"',
  "head_twitch_histogram": 'f"head_twitch_histogram*{experiment}_for_{to_plot}"',}
  but maybe the user would want to pick the name themself? probably better for them to remember what is what in the case of multiple stats choices
