import pandas as pd
import numpy as np
import os, platform, subprocess
from module.core.Dataset import Dataset
from module.core.questions import select_one
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

    def _open_excel(self):
        if self.is_exceled:
            if platform.system() == "Windows":
                os.startfile(self.excel_path)
            elif platform.system() == "Darwin":
                subprocess.call(("open", self.excel_path))
            elif platform.system() == "Linux":
                print("Can't handle Linux")
            else:
                raise OSError("Unknown operating system")
        else:
            raise FileNotFoundError(self.excel_path)

    def edit_excel(self):
        self._open_excel()
        question = "Press any key and ENTER when done editing file"
        user_finished = False
        while not user_finished:
            input(question)
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


@dataclass(repr=False)
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
            for variable in self.project_info["independant_variables"]:
                try:
                    self._template[variable] = [
                        bool(strtobool(str(val))) for val in self._template[variable]
                    ]
                except ValueError as e:
                    raise ValueError(
                        f"{variable} must be a boolean ('True' or 'False')"
                    )
        except Exception as e:
            raise Exception(f"Error creating experiment info : {e}")

    def load(self):
        return (
            super().load().replace(np.nan, "")
        )  # vehicle.independant_var == nan, problem for var in independant_var (nan not iterable)

    @property
    def dict(self):
        """Generates a treatment: infos, group_id: infos mapping
        for easy access
        Returns:
            dict: { treatment_x: {infos_x}, goup_id_x: {infos_x}}
        """
        mapping = {}
        for treatment in self.list:
            mapping[treatment["treatment"]] = treatment
            mapping[treatment["group_id"]] = treatment
        return mapping


@dataclass(repr=False)
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
            for variable in ["paired", "parametric"]:
                try:
                    self._template[variable] = [
                        bool(strtobool(str(val))) for val in self._template[variable]
                    ]
                except ValueError as e:
                    raise ValueError(
                        f"{variable} must be a boolean ('True' or 'False')"
                    )
        except Exception as e:
            raise Exception(f"Error creating project info : {e}")

    def get_experiment(self, name):
        return next(
            filter(
                lambda experiment_information: experiment_information["experiment"]
                == name,
                self.experiments,
            )
        )


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

    _name: ClassVar[str] = "project_information"
    _template = {
        "outlier_test": None,
        "p_value_threshold": None,
        "raw_data_filename": None,
    }

    def __post_init__(self):
        self.filepath = f"{self.location}/{self._name}.json"
        return super().__post_init__()

    def generate(self):
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
