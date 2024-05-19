from __future__ import annotations
from dataclasses import dataclass
from module.core.Experiment import Experiment
from module.core.questions import yes_or_no
from module.core.Metadata import (
    ExperimentInformation,
    TreatmentInformation,
    ProjectInformation,
)
from module.core.Dataset import SelectableDataFrame
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
        self.project_information = ProjectInformation(self.location) 
        self.experiment_information = ExperimentInformation(self.location)
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
        self.experiments = {
            experiment_info["experiment"]: Experiment(self, experiment_info, self.treatment_information)
            for experiment_info in self.experiment_information.list
        }
        self.statistics = Statistics(
            self.location,
            self.experiments,
            self.project_information["p_value_threshold"],
        )        
        
        
    
    @property
    def full_df(self):
        """Get df with full information (values, outliers, treatments)

        Returns:
            pd.datafFrame: A merge of hplc, outliers, and treatment_information, and df.select method
        """
        common_columns = self.hplc.df.columns.intersection(self.outliers.df.columns).tolist()
        full_df = self.hplc.df.merge(self.outliers.df, on=common_columns, how='left').merge(
            self.treatment_information.df, on="group_id", how='left'
        )
        return SelectableDataFrame(full_df)


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