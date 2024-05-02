from __future__ import annotations
from dataclasses import dataclass
from module.core.Experiment import Experiment
from module.core.questions import yes_or_no
from module.core.Metadata import (
    ExperimentInformation,
    TreatmentInformation,
    ProjectInformation,
)
from module.core.HPLC import RawHPLC, HPLC
from module.core.Outliers import Outliers

from module.core.Statistics import Statistics
import os

ROOT = f"{os.getcwd()}/PROJECTS" 
if not os.path.exists(ROOT):
        os.mkdir(ROOT)
@dataclass
class Project:

    name: str
    
    def __post_init__(self):
        self.location = f"{ROOT}/{self.name}"
        if not os.path.exists(self.location):
            if yes_or_no(f"INITIALIZE NEW PROJECT: '{self.name}' ?"):
                os.mkdir(self.location)
            else:
                print(f"UNKNOWN PROJECT: {self.name}")
                print(f"KNOW PROJECTS ARE: {Project.list()}")
                exit(1)
        self.project_information = ProjectInformation(self.location) 
        self.experiment_information = ExperimentInformation(self.location)
        self.experiments = {
            experiment["experiment"]: Experiment(self, experiment)
            for experiment in self.experiment_information.list
        }
        self.treatment_information = TreatmentInformation(self.location)
        self.treatments = self.treatment_information.list

        self.raw_data = RawHPLC(self.location, self.project_information["raw_data_filename"]    )
        self.hplc = HPLC(self.location, self.raw_data)
        self.outliers = Outliers(
            self.location,
            self.project_information["outlier_test"],
            self.project_information["p_value_threshold"],
            self.hplc,
        )
        self.statistics = Statistics(
            self.location,
            self.treatment_information,
            self.experiments,
            self.project_information["p_value_threshold"],
            self.hplc,
            self.outliers,
        )

    @classmethod
    def list(self):
        os.listdir(path=ROOT)

