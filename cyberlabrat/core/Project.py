from __future__ import annotations
from dataclasses import dataclass
from cyberlabrat.core.Dataset import SelectableDataFrame
from cyberlabrat.core.questions import yes_or_no
from cyberlabrat.core.Metadata import (
    ExperimentInformation,
    TreatmentInformation,
    ProjectInformation,
)
from cyberlabrat.core.HPLC import RawHPLC, HPLC, Outliers
from cyberlabrat.core.Statistics import Statistics
import os

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
    def data(self) -> SelectableDataFrame:
        return self.hplc.extend(self.treatment_information).extend(self.outliers)

    @classmethod
    def list(self):
        os.listdir(path=ROOT)


    def _repr_html_(self) -> str:
        return f"""
            <h3>Project: {self.name}</h3>
            <h4>Project Information</h4>
            {self.project_information._repr_html_()}
            <h4>Experiments</h4>
            {self.experiment_information._repr_html_()}
            <h4>Treatments</h4>
            {self.treatment_information._repr_html_()}
        """
    