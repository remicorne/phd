import pandas as pd
import numpy as np
import os
from module.core.Dataset import Dataset
from module.core.Question import Question
from module.core.JSON import JSONMapping
from dataclasses import dataclass
from typing import ClassVar
from distutils.util import strtobool
from module.core.Outliers import OUTLIER_TESTS


class _ProjectSettings(Dataset):
    
    def generate(self):
        return pd.DataFrame(self._template)
    
    def validate(self, data):
        if data.isnull().any():
            raise ValueError("NO EMPTY CELS ALLOWED")
        
    def edit_excel(self):
        self._open_excel()
        question = "Press any key and ENTER when done editing file"
        user_finished = False
        while not user_finished:
            Question.input(question)
            try:
                pd.read_excel(self.excel_path)
                user_finished = True
            except SystemExit:
                self.delete()
                print("System interuption, deleting file")
            except:
                question = "Error reading file. Press any key and ENTER when done correcting file"
                user_finished = False

    def initialize(self):
        super().initialize()
        self.edit_excel()
        
    def save(self, data):
        print(f"Saving {self._name} dataframe")
        data.to_excel(self.excel_path)
        print(f"Saved {self._name} datarame to {self.excel_path}")
    
    def load(self):
        return pd.read_excel(self.excel_path)
       
    @property
    def is_saved(self):
        return self.is_exceled
        
    @property
    def list(self):
        return self.df.to_dict(orient="index").values()



@dataclass
class TreatmentInformation(_ProjectSettings):

    _name: ClassVar[str] = "treatment_information"
    _template: ClassVar[list] = {
        "group_id": [1, 5, 3, 4],
        "treatment": ["vehicles", "MDL", "TCB2", "TCB2+MDL"],
        "color": ["white", "pink", "orange", "red"],
        "independant_variables": ["", "MDL", "TCB2", "TCB2, MDL"],
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


    def load(self):
        return super().load().replace(np.nan, "") # vehicle.independant_var == nan, problem for var in independant_var (nan not iterable)
    
    
@dataclass
class ExperimentInformation(_ProjectSettings):
    _name: ClassVar[str] = "experiment_information"
    _template: ClassVar[list] = {  
        "experiment": ["agonist antagonist"],
        "groups": ["1, 5, 3, 4"],
        "independant_variables": ["TCB2, MDL"],
        "paired": [False],
        "parametric": [True],
    }

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

    def get_experiment(self, name):
        return next(filter(lambda experiment_information: experiment_information['experiment'] == name, self.experiments))
    

@dataclass
class OutlierInformation(JSONMapping):
    
    name: str = "statistics_information"
    _template = {
        "outlier_test": None,
        "p_value_threshold": None,
    }
    
    def generate(self):
        data = self._template
        data["outlier_test"] = Question.select_one("Select outlier test", OUTLIER_TESTS.keys())
        data["p_value_threshold"] = float(Question.input("Enter p value threshold"))
        return data
    
    