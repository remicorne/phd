from module.getters import getCompoundAndRatiosDf, getHeadTwitchDf
from module.utils import CACHE_DIR, askYesorNo, saveJSON, getJSON
import os


treatment_mapping = (
    {  # TO DO : change name to treatment_info and add columns in df REMI
        1: {
            "treatment": "vehicles",
            "color": "white",
            "experiments": ["dose_response", "agonist_antagonist"],
        },
        2: {
            "treatment": "0.3mg/kgTCB",
            "color": "firebrick",
            "experiments": ["dose_response"],
        },
        3: {
            "treatment": "3mg/kgTCB",
            "color": "red",
            "experiments": ["dose_response", "agonist_antagonist"],
        },
        4: {
            "treatment": "10mg/kgTCB",
            "color": "salmon",
            "experiments": ["dose_response"],
        },
        5: {
            "treatment": "0.2mg/kgMDL",
            "color": "black",
            "experiments": ["agonist_antagonist"],
        },
        6: {
            "treatment": "TCB+MDL",
            "color": "grey",
            "experiments": ["agonist_antagonist"],
        },
    }
)

experimental_info = {
    "dose_response": {
        "groups": [1, 2, 3, 4],
        "independant_vars": ["TCB2"],
        "paired": False,
        "parametric": True,  # idealy this would be True / False / Check : check would involve checking the data using the spearman test which should already be done then taking the majority of the data to be parametric or not chose that
        "outliers": ["grubbs"],
        "quantitative_statistics": {
            "twoway_anova": False,  # DELETE ONCE STATS LOGIC IN
            "oneway_anova": True,
            "tukey": True,
        },
    },
    "agonist_antagonist": {
        "groups": [
            1,
            3,
            5,
            6,
        ],
        "independant_vars": ["TCB2", "MDL"],
        "paired": False,
        "parametric": True,
        "outliers": ["grubbs"],
        "quantitative_statistics": {
            "twoway_anova": True,  # DELETE ONCE STATS LOGIC IN
            "oneway_anova": True,
            "tukey": True,
        },
    },
}

from module.Filesystem import Filesystem


class Experiment:
    configuration_names = [
        "treatment_mapping",
        "experimental_info",
        "compound_ratio_mapping",
        "region_subclassification",
        "compound_subclassification",
    ]
    raw_data_names = ["HPLC", "HT"]

    def __init__(self, name):
        self.name = name
        Filesystem.initiate_experiment(name)
        self.data = getCompoundAndRatiosDf(f"{name}_HPLC.csv")
        self.ht = getHeadTwitchDf(f"{name}_HT.csv")


Experiment("TCB2", treatment_mapping, experimental_info)
