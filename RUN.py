from module.quantitative import (
    justStats
)

from module.utils import subselectDf

from module.getters import (
    getQuantitativeStats,
)


filename = "TCB2_data_HPLC.csv"

justStats(filename, 
          experiments=['dose_response'], 
          compounds=['5HIAA/5HT'], 
          regions=["OF","PL","CC", "IC","M", "SJ","SL1", "SL6", "SR6", "SR1", "AC", "V",  
                "Am", "dH", "vH", "NAc", "VM", "DM","VL", "DL", "MD",  "VPL",  "VPR", 
                "DG", "Y",  "SC","SN", "VTA", "DR","MR", "CE"], 
        p_value_threshold=0.05)

subselectDf(getQuantitativeStats(filename), {'experiment':'dose_response', 
                                            #  'is_significant':True, 
                                             'compound':'5HIAA/5HT', 
                                             'test':'one_way_anova',
                                             'region':["VPR"]})