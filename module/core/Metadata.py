import pandas as pd
import numpy as np
import os
from module.core.Dataset import ExcelDataset, SelectableDataFrame
from module.core.questions import select_one
from dataclasses import dataclass, field
from typing import ClassVar
from distutils.util import (
    strtobool,
)  # Deprecated 3.12 https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python


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
        except KeyError as e:
            raise ValueError(f"Missing column {e}")
        except ValueError as e:
            print(e)
            raise ValueError(f"Wrong data types, please correct {self.filename}")

    def convert_dtypes(self, df):
        for col_name, col_info in self._template_types.items():
            if col_info["type"] == list:
                df[col_name] = df[col_name].apply(
                    lambda val: [
                        col_info["subtype"](subval)
                        for subval in (val.replace(" ", "").split(",") if val else [])
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

    def select(self, **selector) -> SelectableDataFrame:
        df = super().select(**selector)
        return df.iloc[0] if len(df) == 1 else df


@dataclass(repr=False)
class GroupInformation(_ProjectSettings):  # TODO: generalize to GroupInformation?

    filename: ClassVar[str] = "group_information"
    _template: ClassVar[dict] = {
        "group_id": [1, 5, 3, 4],
        "label": ["vehicles", "MDL", "TCB2", "TCB2+MDL"],
        "independant_variables": ["", "MDL", "TCB2", "TCB2, MDL"],
        "mouse_id": [
            "2, 5, 7, 9, 11, 17, 20, 28, 32, 59, 67",
            "13, 14, 15, 26, 29, 34, 42, 48, 63, 65",
            "23, 24, 31, 36, 38, 40, 44, 50, 51, 57",
            "21, 25, 35, 41, 45, 49, 52, 58, 61, 66, 69",
        ],
    }
    _template_types: ClassVar[dict] = {
        "group_id": {"type": int},
        "label": {"type": str},
        "independant_variables": {"type": list, "subtype": str},
        "mouse_id": {"type": list, "subtype": int},
    }
    control_group: ClassVar[list] = "vehicles"

    @property
    def palette(self):
        return {t.label: t.color for t in self}

    @property
    def treatments(self):
        return list(self.df.label)

    def extend_dataset(self, dataset):
        return self.df.explode("mouse_id").extend(dataset)


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

    def get_significance_palette(self):
        return self._get_palette("significance")

    def get_color_palette(self):
        return self._get_palette("color")

    def _get_palette(self, palette_type):
        return {
            row.treatment: row[palette_type]
            for _, row in Palette(self.project).df.iterrows()
        }

    def __contains__(self, value):
        return value in self.df.treatment

    def __getitem__(self, treatment) -> pd.Series:
        return self.df.select(**{"treatment": treatment})


@dataclass(repr=False)
class ExperimentInformation(_ProjectSettings):

    filename: ClassVar[str] = "experiment_information"
    _template: ClassVar[dict] = {
        "label": ["agonist_antagonist"],
        "groups": ["1, 5, 3, 4"],
        "independant_variables": ["TCB2, MDL"],
        "paired": [False],
        "parametric": [True],
        "control_group_id": [1],
        # "data_source": ["hplc, behavior"],
    }
    _template_types: ClassVar[dict] = {
        "label": {"type": str},
        "groups": {"type": list, "subtype": int},
        "independant_variables": {"type": list, "subtype": str},
        "paired": {"type": bool},
        "parametric": {"type": bool},
        "control_group_id": {"type": int},
        # "data_source": {"type": list, "subtype": str},
    }

    def load(self):
        data = super().load()
        group_information = GroupInformation(self.project).df
        palette = Palette(self.project).get_color_palette()
        full_experiment_info = []
        for _, experiment in data.iterrows():
            experiment["experiment"] = experiment.label
            experiment["treatments"] = group_information.select(
                group_id=experiment.groups
            ).label.to_list()
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
        "subject_column": ["mouse_id"],
        "group_column": ["group_id"],
        "grouping_characteristic": ["treatment"],
    }
    _template_types: ClassVar[dict] = {
        "label": {"type": str},
        "outlier_test": {"type": str},
        "p_value_threshold": {"type": float},
        "raw_data_filename": {"type": str},
        "subject_column": {"type": str},
        "group_column": {"type": str},
        "grouping_characteristic": {"type": str},
    }

    def generate(self):
        from module.core.HPLC import OUTLIER_TESTS

        data = super().generate()
        data["label"] = self.project
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
        "label": ["hplc", "tissue_weight", "behavior"],
        "measurement_columns": ["compound, region", "region", "measure"],
        "unit": ["ng/mg", "mg", ""],
        "experiments": [
            "agonist_antagonist, dose_response",
            "agonist_antagonist, dose_response",
            "",
        ],
    }
    _template_types: ClassVar[dict] = {
        "label": {"type": str},
        "measurement_columns": {"type": list, "subtype": str},
        "unit": {"type": str},
        "experiments": {"type": list, "subtype": str},
    }
