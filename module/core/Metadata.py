import pandas as pd
import numpy as np
import os
from module.core.Dataset import ExcelDataset
from module.core.questions import select_one
from module.core.JSON import JSONMapping
from dataclasses import dataclass
from typing import ClassVar
from distutils.util import strtobool


@dataclass(repr=False)
class _ProjectSettings(ExcelDataset):

    project: str
    _template: ClassVar[dict] = None
    _template_types: ClassVar[dict] = None

    def generate(self):
        return pd.DataFrame(self._template)

    def validate(self, data):
        super().validate(data)

    def make_user_edit_excel(self):
        self.open()
        question = "Press any key and ENTER when done editing file"
        user_finished = False
        while not user_finished:
            input(question)
            try:
                self.validate(self.load())
                user_finished = True
            except SystemExit:
                self.delete()
                print("System interuption, deleting file")
            except:
                question = "Error reading file. Press any key and ENTER when done correcting file"
                user_finished = False

    def initialize(self):
        super().initialize()
        self.make_user_edit_excel()

    def load(self):
        # vehicle.independant_var == nan, problem for "var in independant_var" (nan not iterable)
        data = super().load().replace(np.nan, "")
        try:
            return self.convert_dtypes(data)
        except:
            raise ValueError(f"Wrong data types, please correct {self.filename}")

    def convert_dtypes(self, df):
        for col_name, col_info in self._template_types.items():
            if col_info["type"] == list:
                df[col_name] = df[col_name].apply(lambda val: [col_info["subtype"](subval) for subval in val.replace(" ", "").split(",")])
            elif col_info["type"] == bool:
                df[col_name] = df[col_name].apply(lambda val: bool(strtobool(str(val))))
            else:
                df[col_name] = df[col_name].apply(col_info["type"])
        return df
@dataclass(repr=False)
class TreatmentInformation(_ProjectSettings): # TODO: generalize to GroupInformation?

    filename: ClassVar[str] = "treatment_information"
    _template: ClassVar[dict] = {
        "group_id": [1, 5, 3, 4],
        "treatment": ["vehicles", "MDL", "TCB2", "TCB2+MDL"],
        "color": ["white", "pink", "orange", "red"],
        "independant_variables": ["", "MDL", "TCB2", "TCB2, MDL"],
    }
    _template_types: ClassVar[dict] = {
        "group_id": {"type": int},
        "treatment": {"type": str},
        "color": {"type": str},
        "independant_variables": {"type": list, "subtype": str},
    }
    
    def get_palette(self):
        return {group["treatment"]: group["color"] for group in self.list}

    # @property
    # def dict(self):
    #     """Generates a treatment: infos, group_id: infos mapping
    #     for easy access
    #     Returns:
    #         dict: { treatment_x: {infos_x}, goup_id_x: {infos_x}}
    #     """
    #     mapping = {}
    #     for treatment in self.list:
    #         mapping[treatment["treatment"]] = treatment
    #         mapping[treatment["group_id"]] = treatment
    #     return mapping


@dataclass(repr=False)
class ExperimentInformation(_ProjectSettings):

    filename: ClassVar[str] = "experiment_information"
    _template: ClassVar[dict] = {
        "experiment": ["agonist antagonist"],
        "groups": ["1, 5, 3, 4"],
        "independant_variables": ["TCB2, MDL"],
        "paired": [False],
        "parametric": [True],
    }
    _template_types: ClassVar[dict] = {
        "experiment": {"type": str},
        "groups": {"type": list, "subtype": int},
        "independant_variables": {"type": list, "subtype": str},
        "paired": {"type": bool},
        "parametric": {"type": bool},
    }
    
    @property
    def experiments(self):
        return self.df.experiment.tolist()

    def get_experiment(self, name):
        return next(
            filter(
                lambda experiment_information: experiment_information["experiment"]
                == name,
                self.list,
            )
        )
        
    def get_experiments(self, group_id):
        return list(self.df[self.df.groups.apply(lambda groups: group_id in groups)].experiment)
    
    
    # def label(self, df):
        

def is_valid_file(file_path):
    if not os.path.isfile(file_path):
        print("Not found", file_path)
        return False
    extension = os.path.splitext(file_path)[1].lower()
    if extension not in [".xlsx", ".csv"]:
        print("Invalid extension:", extension)
        return False

    return True


@dataclass(repr=False)
class ProjectInformation(JSONMapping):

    project: str
    filename: ClassVar[str] = "project_information"
    _template = {
        "outlier_test": None,
        "p_value_threshold": None,
        "raw_data_filename": None,
    }

    def generate(self):
        from module.core.Outliers import OUTLIER_TESTS

        data = self._template
        data["outlier_test"] = select_one("Select outlier test", OUTLIER_TESTS.keys())
        data["p_value_threshold"] = float(input("Enter p value threshold"))
        data["raw_data_filename"] = self.get_valid_filename()
        return data

    def get_valid_filename(self):
        # raw_data_filename = easygui.fileopenbox(title="Select raw HPLC file", filetypes=["*.xls", "*.xlsx"])
        raw_data_filename = input("Enter HPLC excel filename (must be in phd/)")
        file_path = f"{os.getcwd()}/{raw_data_filename}"

        while not is_valid_file(file_path):
            print(raw_data_filename, "NOT FOUND")
            raw_data_filename = input("Enter excel HPLC filename (must be in phd/)")
            file_path = f"{os.getcwd()}/{raw_data_filename}"

        return file_path

    @property
    def outlier_test(self):
        return self.dict["outlier_test"]

    @property
    def p_value_threshold(self):
        return self.dict["p_value_threshold"]

    @property
    def raw_data_filename(self):
        return self.dict["raw_data_filename"]
