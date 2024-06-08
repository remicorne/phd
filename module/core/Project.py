from __future__ import annotations
from dataclasses import dataclass
from module.core.Experiment import Experiment
from module.core.questions import yes_or_no
from module.core.Metadata import (
    ExperimentInformation,
    TreatmentInformation,
    ProjectInformation,
)
from module.core.FullHPLC import FullHPLC
from module.core.HPLC import RawHPLC, HPLC
from module.core.Outliers import Outliers

from module.core.Statistics import Statistics
import os
import types

ROOT = f"{os.getcwd()}/PROJECTS" 
if not os.path.exists(ROOT):
        os.mkdir(ROOT)
@dataclass
class Project:
    """
    A class to manage and organize different components of a research project,
    including handling of experimental setups, data collection, and analysis.

    Attributes:
        name (str): The name of the project.
        location (str): The file path to the project's directory.
        project_information (ProjectInformation): Contains metadata about the project.
        experiment_information (ExperimentInformation): Handles the collection of experiments.
        experiments (dict): Stores instantiated Experiment objects, keyed by experiment names.
        treatment_information (TreatmentInformation): Manages information related to treatments in experiments.
        treatments (list): A list of treatments derived from treatment_information.
        raw_data (RawHPLC): An object handling raw high-performance liquid chromatography data.
        hplc (HPLC): An object to manage processed HPLC data.
        outliers (Outliers): Manages outlier detection based on project settings.
        statistics (Statistics): Handles statistical analysis of experimental data.
        full_df: the df with all values + outliers+ treatment nfo

    Methods:
        __post_init__(): Initializes a new Project instance, sets up necessary directories,
                         loads or prompts for project information, and prepares all related
                         components like experiments, treatments, and statistical analyses.
    """

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
        self.project_information = ProjectInformation(self.name)
        self.experiment_information = ExperimentInformation(self.name)
        self.treatment_information = TreatmentInformation(self.name)
        self.raw_data = RawHPLC(self.name)
        self.hplc = HPLC(self.name)
        self.outliers = Outliers(self.name)
        self.statistics = Statistics(self.name)        
        
        
    
    @property
    def full_df(self):
        """Get df with full information (values, outliers, treatments)

        Returns:
            pd.datafFrame: A merge of hplc, outliers, and treatment_information, and df.select method
        """
        
        return FullHPLC(self.name).df


    @classmethod
    def list(self):
        os.listdir(path=ROOT)


    def __repr__(self) -> str:
        experiments_string = '\n'.join(str(experiment )for experiment in self.experiments.values())
        return f""""
    Project: {self.name}
    
    Parameters
    {self.project_information}
    
    Experiments
    {experiments_string}"""
    
    
    
from distutils.util import strtobool

@dataclass
class Experiment:

    project: str
    name: str

    def __post_init__(self):
        experiment_information = ExperimentInformation(self.project).get_experiment(self.name)
        treatment_information = TreatmentInformation(self.project)
        self.name = experiment_information["experiment"]
        self.groups = [
            int(group)
            for group in experiment_information["groups"]
            .replace(" ", "")
            .split(",")
        ]
        self.treatments = [treatment_information.dict[group_id]["treatment"] for group_id in self.groups]
        self.palette = {
            treatment: treatment_information.dict[treatment]["color"]
            for treatment in self.treatments
        }
        self.independant_variables = (
            experiment_information["independant_variables"]
            .replace(" ", "")
            .split(",")
        )
        self.paired = (
            experiment_information["paired"]
            if isinstance(experiment_information["paired"], bool)
            else strtobool(experiment_information["paired"])
        )
        self.parametric = (
            experiment_information["parametric"]
            if isinstance(experiment_information["parametric"], bool)
            else strtobool(experiment_information["parametric"])
        )

    @property
    def df(self):
        """Get df with full information (values, outliers, treatments)

        Returns:
            pd.datafFrame: A merge of hplc, outliers, and treatment_information, and df.select method
        """
        common_columns = self.hplc.df.columns.intersection(self.outliers.df.columns).tolist()
        full_df = self.hplc.df.merge(self.outliers.df, on=common_columns, how='left').merge(
            self.treatment_information.df, on="group_id", how='left'
        )
        full_df.loc[:, "experiment"] = self.name
        full_df[self.independant_variables] = list(
            full_df.independant_variables.apply(
                lambda group_independant_variables: [
                    experiment_independant_variable in group_independant_variables
                    for experiment_independant_variable in self.independant_variables
                ],
            )
        )
        return full_df

    def __repr__(self) -> str:
        return f"""
        groups: {self.groups}
        independant variables: {self.independant_variables}
        paired: {self.paired}
        parametric: {self.parametric}
        """
