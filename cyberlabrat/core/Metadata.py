import pandas as pd
import numpy as np
import os
from cyberlabrat.core.Dataset import ExcelDataset, SelectableDataFrame
from cyberlabrat.core.questions import select_one
from dataclasses import dataclass, field
from typing import ClassVar
from distutils.util import strtobool # Deprecated 3.12 https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python


@dataclass(repr=False)
class _ProjectSettings(ExcelDataset):
    """Base class for project settings.
    Handles loading, saving, and editing of project settings (excel files)

    Returns:
        ExcelDataset: Dataset with project settings
    """

    project: str = field(default=None)
    _template: ClassVar[dict] = None
    _template_types: ClassVar[dict] = None
    
    def __post_init__(self):
        """
        Load data from template and set columns as attributes
        """
        super().__post_init__()
        for key in self._template:
            self.__setattr__(key, self.df[key])

    def generate(self):
        """
        Generate template dataframe

        """
        return pd.DataFrame(self._template)

    def make_user_edit_excel(self):
        """
        Open excel file in edit mode to make user edit the config
        """
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
        """
        Same as chacheable initialize but also makes user edit as setting are user editable
        """
        super().initialize()
        self.make_user_edit_excel()

    def load(self):
        """
        Load data from template and handles type conversions for merges with raw data.
        This process makes sure the template is both human and programatic friendly.

        Raises:
            ValueError: If a cell is not of the correct type (ie user input unusable data)

        Returns:
            SelectableDataFrame: Contains the project settings
        """
        # vehicle.independant_var == nan, problem for "var in independant_var" (nan not iterable)
        data = super().load().replace(np.nan, "")
        try:
            return self.convert_dtypes(data)
        except:
            raise ValueError(f"Wrong data types, please correct {self.filename}")

    def convert_dtypes(self, df):
        for col_name, col_info in self._template_types.items():
            if col_info["type"] == list:
                df[col_name] = df[col_name].apply(
                    lambda val: [
                        col_info["subtype"](subval)
                        for subval in val.replace(" ", "").split(",")
                    ]
                )
            elif col_info["type"] == bool:
                df[col_name] = df[col_name].apply(lambda val: bool(strtobool(str(val))))
            else:
                df[col_name] = df[col_name].apply(col_info["type"])
        return df
    
    
    
    def __contains__(self, label):
        return label in self.df.label        


    def __getitem__(self, label) -> pd.Series:
        return self.df.select(**{"label": label})
    


@dataclass(repr=False)
class TreatmentInformation(_ProjectSettings):  # TODO: generalize to GroupInformation?

    filename: ClassVar[str] = "treatment_information"
    _template: ClassVar[dict] = {
        "group_id": [1, 5, 3, 4],
        "label": ["vehicles", "MDL", "TCB2", "TCB2+MDL"],
        "independant_variables": ["", "MDL", "TCB2", "TCB2, MDL"],
    }
    _template_types: ClassVar[dict] = {
        "group_id": {"type": int},
        "label": {"type": str},
        "independant_variables": {"type": list, "subtype": str},
    }
    control_group: ClassVar[list] = "vehicles"

    @property
    def palette(self):
        return {t.label: t.color for t in self}

    def get_hue_order(self):
        return self.df.sort_values(by="group_id").treatment.tolist()
    
    def load(self):
        data = super().load()
        data["treatment"] = data.label
        return data
    
    @property
    def treatments(self):
        return list(self.df.label)
    
    
@dataclass(repr=False)
class Palette(_ProjectSettings):  # TODO: generalize to GroupInformation?

    filename: ClassVar[str] = "palette"
    _template: ClassVar[dict] = {
        "treatment": ["vehicles", "MDL", "TCB2", "TCB2+MDL"],
        "color": ["white", "pink", "orange", "red"],
        "significance": ["*", "", "$", ""],
    }
    _template_types: ClassVar[dict] = {
        "treatment": {"type": str},
        "color": {"type": str},
        "significance": {"type": str},
    }

    @property
    def dict(self):
        return {t.treatment: t.color for t in self}
    
    
    def __contains__(self, value):
        return value in self.df.treatment
    
    def __getitem__(self, treatment) -> pd.Series:
        return self.df.select(**{"treatment": treatment})
    


@dataclass(repr=False)
class ExperimentInformation(_ProjectSettings):

    filename: ClassVar[str] = "experiment_information"
    _template: ClassVar[dict] = {
        "label": ["agonist antagonist"],
        "groups": ["1, 5, 3, 4"],
        "independant_variables": ["TCB2, MDL"],
        "paired": [False],
        "parametric": [True],        
        "control_group_id": [1]
    }
    _template_types: ClassVar[dict] = {
        "label": {"type": str},
        "groups": {"type": list, "subtype": int},
        "independant_variables": {"type": list, "subtype": str},
        "paired": {"type": bool},
        "parametric": {"type": bool},
        "control_group_id": {"type": int}
    }        

    def load(self):
        data = super().load()
        treatment_information = TreatmentInformation(self.project).df
        palette = Palette(self.project).dict
        full_experiment_info = []
        for _, experiment in data.iterrows():
            experiment["experiment"] = experiment.label
            experiment["treatments"] = treatment_information.select(group_id=experiment.groups).label.to_list()
            experiment["control_treatment"] = experiment.treatments[experiment.groups.index(experiment["control_group_id"])]
            experiment["palette"] = {t: palette[t] for t in experiment.treatments}
            full_experiment_info.append(experiment)
        return SelectableDataFrame(full_experiment_info)
    
    
    @property
    def experiments(self):
        return list(self.df.label)

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
class ProjectInformation(_ProjectSettings):

    filename: ClassVar[str] = "project_information"
    _template: ClassVar[dict] = {
        "label": ["TCB2"],
        "outlier_test": ["grubbs"],
        "p_value_threshold": [0.05],
        "raw_data_filename": ["raw_data.csv"],
    }
    _template_types: ClassVar[dict] = {
        "label": {"type": str},
        "outlier_test": {"type": str},
        "p_value_threshold": {"type": float},
        "raw_data_filename": {"type": str},
    }
    
    def generate(self):
        from cyberlabrat.core.HPLC import OUTLIER_TESTS
        data = super().generate()
        data["label"] = input("Enter project name")
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
    def df(self):
        data = super().df
        return data.iloc[0] if len(data) == 1 else data

    def _repr_html_(self) -> str:
        return f"""
            <b>Project Information:</b><br>
            <ul>
                <li>Project: {self.project}</li>
                <li>Raw Data Filename: {self.raw_data_filename}</li>
                <li>Outlier Test: {self.outlier_test}</li>
                <li>P-Value Threshold: {self.p_value_threshold}</li>
            </ul>
        """


@dataclass(repr=False)
class DatasetInformation(_ProjectSettings):

    filename: ClassVar[str] = "dataset_information"
    _template: ClassVar[dict] = {
        "label": ["hplc", "tissue_weight"],
        "grouping_variables": ["treatment, compound, region", "treatment, region"],
        "independant_variables": ["TCB2, MDL", ""],
        "raw_data_filename": ["raw_data.csv"],
    }
    _template_types: ClassVar[dict] = {
        "label": {"type": str},
        "outlier_test": {"type": str},
        "p_value_threshold": {"type": float},
        "raw_data_filename": {"type": str},
    }