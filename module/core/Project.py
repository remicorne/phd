from __future__ import annotations
from dataclasses import dataclass
from module.core.Experiment import Experiment
from module.core.Question import Question
from module.core.Metadata import (
    ExperimentInformation,
    TreatmentInformation,
    OutlierInformation,
)
from module.core.HPLC import RawHPLC, HPLC
from module.core.Outliers import Outliers

from module.core.Statistics import Statistics
from ..constants import ROOT
import os
from distutils.util import strtobool


PROJECTS = f"{ROOT}/PROJECTS"


@dataclass
class Project:

    name: str

    def __post_init__(self):
        self.location = f"{PROJECTS}/{self.name}"
        if not os.path.exists(self.location):
            if Question.yes_or_no(f"INITIALIZE NEW PROJECT: '{self.name}' ?"):
                os.mkdir(self.location)
            else:
                print(f"UNKNOWN PROJECT: {self.name}")
                print(f"KNOW PROJECTS ARE: {Project.list()}")
                exit(1)
        self.statistics_information = OutlierInformation(self.location)
        self.outlier_test = self.statistics_information["outlier_test"]
        self.p_value_threshold = self.statistics_information["p_value_threshold"]

        self.raw_data = RawHPLC(self.location)
        self.hplc = HPLC(self.location, self.raw_data)
        self.outliers = Outliers(
            self.location,
            self.outlier_test,
            self.p_value_threshold,
            self.hplc,
        )

        self.experiment_information = ExperimentInformation(self.location)
        self.experiments = {
            experiment["experiment"]: Experiment(self, experiment)
            for experiment in self.experiment_information.list
        }
        self.treatment_information = TreatmentInformation(self.location)
        self.treatments = self.treatment_information.list

        self.statistics = Statistics(
            self.location,
            self.treatment_information,
            self.experiments,
            self.p_value_threshold,
            self.hplc,
            self.outliers,
        )

    # def correlogram(self, )

    @classmethod
    def list():
        os.listdir(path=PROJECTS)


### TODO, factorize with experiment and seprat directory from data

if not os.path.exists(PROJECTS):
    os.mkdir(PROJECTS)


@dataclass
class Experiment:

    project: Project
    experiment_information: dict

    def __post_init__(self):
        self.name = self.experiment_information["experiment"]
        self.groups = [
            int(group)
            for group in self.experiment_information["groups"]
            .replace(" ", "")
            .split(",")
        ]
        self.independant_variables = (
            self.experiment_information["independant_variables"]
            .replace(" ", "")
            .split(",")
        )
        self.paired = (
            self.experiment_information["paired"]
            if isinstance(self.experiment_information["paired"], bool)
            else strtobool(self.experiment_information["paired"])
        )
        self.parametric = (
            self.experiment_information["parametric"]
            if isinstance(self.experiment_information["parametric"], bool)
            else strtobool(self.experiment_information["parametric"])
        )
        self.location = f"{self.project.location}/{self.name}"
        del self.experiment_information
