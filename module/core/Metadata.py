import pandas as pd
from module.core.Dataset import Dataset
from dataclasses import dataclass
from typing import ClassVar
from typing import List


class _ProjectSettings(Dataset):

    def validate(self, data):
        if data.isnull().any():
            raise ValueError("NO EMPTY CELS ALLOWED")

    def initialize(self):
        super().initialize()
        self.edit_excel_and_save()        
        

@dataclass
class ExperimentInformation(_ProjectSettings):
    experiment: str
    _name: ClassVar[str] = "experiment_info"
    _template: ClassVar[list] = {
        "treatment": [],
        "group_id": [],
        "color": [],
    }

    def __post_init__(self):
        self._name = f"{self.experiment}_{self._name}"
        return super().__post_init__()
    
    def generate(self):
        project_info = ProjectInformation(self.project).get_experiment(self.experiment)
        for variable in project_info['independant_variables']:
            self._template[variable] = []
            
        for treatment in project_info['treatments'].split(","):
            self._template['treatment'].append(treatment)
            self._template['group_id'].append(f"{treatment}_group_number")
            self._template['color'].append(f"{treatment}_color")
            for variable in self.independant_variables:
                self._template[variable].append(f"{treatment}_{variable} present/absent (1/0)")
    
    @property
    def treatments(self):
        data = self.get()
        return data.to_dict('records')
        
    def get_treatment(self, name):
        return next(filter(lambda treatment: treatment['treatment'] == name, self.treatments))



class ProjectInformation(_ProjectSettings):
    _name: ClassVar[str] = "project_info"
    _template: ClassVar[list] = {
        "experiment": ["agonist antagonist", "dose response"],
        "number_of_treatments": [4, 4],
        "treatments": [
            "vehicles, TCB, TCB+MDL, MDL",
            "vehicles, 0.3mg/kg, 3mg/kg, 10mg/kg",
        ],
        "independant_variables": ["TCB2, MDL", "TCB2"],
        "paired": [False, False],
        "parametric": [True, True],
        "outliers": ["grubbs", "grubbs"],
    }

    def generate(self):
        return pd.DataFrame(self._template)
    
    def validate(self, data):
        super().validate()
        errors = data.filter(lambda row: len(row.treatments.split('')) != row.number_of_treatments)
        if not errors.empty:
            raise ValueError(', '.join([f'{row.treatments}: {row.number_of_treatments} EXPECTED' for row in errors]))
    
    def initialize(self):
        super().initialize()
        experiment_infos = []
        for experiment in self.experiments:
            try:
                experiment_infos.append(ExperimentInformation(experiment))
            except SystemExit:
                for experiment_info in experiment_infos:
                    experiment_info.delete()
                self.delete()
        
    @property
    def experiments(self):
        data = self.get()
        return data.to_dict('records')
        
    
    def get_experiment(self, name):
        return next(filter(lambda experiment: experiment['experiment'] == name, self.experiments))
    
