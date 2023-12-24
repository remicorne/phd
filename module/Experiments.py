# This file contains the Experiment class, which is the main class of the module.

from module.Datasets import HPLC, HT
from module.Filesystem import Filesystem
from module.Configuration import Configuration


class Experiment:
    """
    Contain everything ralated to an experiment
    """

    configuration_names = [
        "treatment_mapping",
        "experimental_info",
        "compound_ratio_mapping",
        "region_subclassification",
        "compound_subclassification",
    ]
    raw_data_names = ["HPLC", "HT"]

    def __init__(self, name):
        self.initialise_experiment(name)
        self.name = name  # TCB2
        self.config = Configuration(name)
        self.HPLC = HPLC(name, self.config)
        self.HT = HT(name, self.config)

    def initialise_experiment(self, experiment_name):
        """Test

        Args:
            experiment_name (_type_): _description_
        """
        Filesystem.check_directory(f"{Filesystem.INPUT}/{experiment_name}")
        Filesystem.check_directory(f"{Filesystem.OUTPUT}/{experiment_name}")
