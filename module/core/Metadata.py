import pandas as pd
import os
from module.core.Dataset import Dataset
from module.core.Questions import Questions
from dataclasses import dataclass
from typing import ClassVar
from distutils.util import strtobool



class _ProjectSettings(Dataset):
    
    
    def load(self):
        return pd.read_excel(self.excel_path)
    
    def validate(self, data):
        if data.isnull().any():
            raise ValueError("NO EMPTY CELS ALLOWED")
        
    def edit_excel(self):
        self._open_excel()
        question = "Press any key and ENTER when done editing file"
        user_finished = False
        while not user_finished:
            Questions.input(question)
            try:
                pd.read_excel(self.excel_path)
                user_finished = True
            except SystemExit:
                self.delete()
                print("System interuption, deleting file")
            except:
                question = "Error reading file. Press any key and ENTER when done correcting file"
                user_finished = False

    def generate(self):
        super().generate().to_excel(self.excel_path)

    def initialize(self):
        self.generate()
        self.edit_excel()
        

@dataclass
class TreatmentInformation(_ProjectSettings):
    experiment: str
    _name: ClassVar[str] = "experiment_info"
    _template: ClassVar[list] = {
        "group_id": [1, 5],
        "treatment": ["vehicles", "TCB2+MDL"],
        "color": ["white", "red"],
        "independant_variables": ["", "TCB2, MDL"],
    }
    
    def __post_init__(self):
        if not os.path.exists(self.location):
            os.mkdir(self.location)
        super().__post_init__()        
    
    def validate(self, data):
        try:
            super().validate()
            for variable in self.project_info['independant_variables']:
                try:
                    self._template[variable] = [bool(strtobool(str(val))) for val in self._template[variable]]
                except ValueError as e:
                    raise ValueError(f"{variable} must be a boolean ('True' or 'False')")
        except Exception as e:
            raise Exception(f"Error creating experiment info : {e}")


    # def get_treatment(self, name):
    #     return next(filter(lambda treatment: treatment['treatment'] == name, self.df.treatments.tolist()))
    
    

class ExperimentInformation(_ProjectSettings):
    _name: ClassVar[str] = "project_info"
    _template: ClassVar[list] = {
        "experiment": ["agonist antagonist"],
        "groups": ["1, 2, 5, 8"],
        "independant_variables": ["TCB2, MDL"],
        "paired": [False],
        "parametric": [True],
    }


    # Outlier is no at project level and not in the df
    def validate(self, data):
        try:
            super().validate()
            for variable in ['paired', 'parametric']:
                try:
                    self._template[variable] = [bool(strtobool(str(val))) for val in self._template[variable]]
                except ValueError as e:
                    raise ValueError(f"{variable} must be a boolean ('True' or 'False')")
        except Exception as e:
            raise Exception(f"Error creating project info : {e}")
    
    @property
    def experiments(self):
        return self.df.experiment.tolist()

    def get_experiment(self, name):
        return next(filter(lambda experiment_information: experiment_information['experiment'] == name, self.experiments))
    
