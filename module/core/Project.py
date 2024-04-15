from __future__ import annotations
from dataclasses import dataclass
from module.core.Experiment import Experiment
from module.core.Questions import Questions
from module.core.Metadata import ExperimentInformation, TreatmentInformation
from module.core.HPLC import HPLC
from ..constants import ROOT
from module.outliers import OUTLIER_TESTS
import os

PROJECTS = f"{ROOT}/PROJECTS"

@dataclass
class Project:

    name: str
    
    def __post_init__(self):
        self.location = f"{PROJECTS}/{self.name}"
        if not os.path.exists(self.location):
            if Questions.yes_or_no(f"INITIALIZE NEW PROJECT: '{self.name}' ?"):
                os.mkdir(self.location)
            else:
                print(f"UNKNOWN PROJECT: {self.name}")
                print(f"KNOW PROJECTS ARE: {Project.list()}")
                exit(1)
        self.hplc = HPLC(self.location)
        self.experiment_information = ExperimentInformation(self.location)
        self.treatment_information = TreatmentInformation(self.location)
        self.experiments = [Experiment(self, experiment_information) for experiment_information in self.experiment_information.iterrows()]
        self.outlier_test =  Questions.select_one(OUTLIER_TESTS.keys(), "Select outlier test")

    ##todo: reset cache deletes pkls but not excel files
    
    # @property
    # def df(self):
    #     return self.hplc.df.merge(self.information.full, on="group_id")

    @classmethod
    def list():
        os.listdir(path=PROJECTS)
        
### TODO, factorize with experiment and seprat directory from data

if not os.path.exists(PROJECTS):
    os.mkdir(PROJECTS)
