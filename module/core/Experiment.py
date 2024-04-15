from __future__ import annotations
import os
from dataclasses import dataclass
from Project import Project, PROJECTS
from Questions import Questions


@dataclass
class Experiment:

    project: Project
    experiment_info: str

    def __post_init__(self):
        self.name = self.experiment_info.experiment
        self.groups = self.experiment_info.groups.replace(' ', '').split(",")
        self.independent_variables = self.experiment_info.independent_variables.replace(' ', '').split(",")
        self.paired = self.experiment_info.paired
        self.parametric = self.experiment_info.parametric
        self.location = f"{self.project.location}/{self.name}"
        if not os.path.exists(self.location):
            if Questions.yes_or_no(f"INITIALIZE NEW EXPERIMENT: '{self.name}' ?"):
                os.mkdir(self.location)
            else:
                print(f"UNKNOWN EXPERIMENT: {self.name}")
                print(f"KNOWN EXPERIMENTS ARE: {os.listdir(self.location)}")
                exit(1)

    @property
    def df(self):
        return self.project.hplc.df.merge(self.information.df, on="group_id")
    
    @classmethod
    def list():
        experiments = {}
        for project in os.listdir(path=PROJECTS):
            experiments.expand([(project.name, experiment.name) for experiment in os.listdir(project.path)])
        return experiments
    
    @classmethod
    def find(cls, experiment):
        possible_projects = [project for project, experiment in Experiment.list() if experiment == 'experiment']
        if possible_projects:
            project = possible_projects[0] if len(possible_projects == 1) else Questions.select_one(possible_projects)
            return cls(experiment, project)
        raise(f"Experiment {experiment} not found in any project")