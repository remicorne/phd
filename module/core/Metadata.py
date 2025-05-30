from dataclasses import dataclass, field
from distutils.util import (
    strtobool,
)  # Deprecated 3.12 https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python
from typing import ClassVar, Dict

import numpy as np
import pandas as pd

from module.core.Dataset import (
    DataframeWrapperMixin,
    ExcelCachedDataFrame,
    SelectableDataFrame,
)


class ValidationError(Exception):
    pass


@dataclass(repr=False)
class ProjectMetadata(ExcelCachedDataFrame):
    """Base class for project settings.
    Handles loading, saving, and editing of project settings (excel files)

    Returns:
        ExcelCachedDataFrame: Dataset with project settings
    """

    project: str = field(default=None)
    filename: ClassVar[str] = "metadata"

    def __post_init__(self):
        super().__post_init__()
        sheets = self.load()
        self.datasets = Datasets(sheets["datasets"])
        self.groups = Groups(sheets["groups"])
        self.experiments = Experiments(sheets["experiments"])
        self.palette = Palette(sheets["palette"])
        self.statistics = Statistics(sheets["statistics"])
        self.validate_consistency()

        self.subject_ids = self.groups.subject_ids
        self.experiments["group_names"] = self.experiments.group_ids.apply(
            lambda x: self.groups.group_name[self.groups.group_id.isin(x)].tolist()
        )
        self.p_value_threshold = self.statistics.p_value_threshold
        self.max_outliers = self.statistics.max_outliers

    def validate_consistency(self):
        experiment_group_ids = set()
        for experiment in self.experiments:
            experiment_group_ids.update(experiment.group_ids)

        groups_group_ids = set(self.groups.group_id)
        palette_group_ids = set(self.palette.group_id)

        unknown_experiment_groups = experiment_group_ids - groups_group_ids
        unknown_palette_groups = palette_group_ids - groups_group_ids

        if unknown_experiment_groups:
            raise ValidationError(
                f"Group IDs in experiments not found in groups: {unknown_experiment_groups}"
            )

        if unknown_palette_groups:
            raise ValidationError(
                f"Group IDs in palette not found in groups: {unknown_palette_groups}"
            )

    def generate(self):
        return {
            "datasets": Datasets.generate(),
            "groups": Groups.generate(),
            "experiments": Experiments.generate(),
            "palette": Palette.generate(),
            "statistics": Statistics.generate(),
        }

    def initialize(self):
        """
        Same as chacheable initialize but also makes user edit as setting are user editable
        """
        super().initialize()
        self.make_user_edit_excel()

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
                self.load()
                user_finished = True
            except SystemExit:
                self.delete()
                print("System interuption, deleting file")
            except Exception as e:
                question = "Error reading file. Press any key and ENTER when done correcting file"
                user_finished = False

    def save(self, content: Dict[str, pd.DataFrame]):
        """Save all sheet in content to multi sheet dataframe using keys as sheet names"""
        with pd.ExcelWriter(self.filepath) as writer:
            for sheet_name, df in content.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    def load(self):
        return pd.read_excel(self.filepath, sheet_name=None)


class SubSetting(DataframeWrapperMixin):
    """Base class for project settings subsets.
    Handles validation and data processing for individual sheets.
    """

    sheet_name: ClassVar[str] = None
    _default: ClassVar[dict] = None
    _types: ClassVar[dict] = None

    def __init__(self, df: pd.DataFrame):
        """
        Load data from template and handles type conversions for merges with raw data.
        This process makes sure the template is both human and programatic friendly.

        Raises:
            ValueError: If a cell is not of the correct type (ie user input unusable data)

        Returns:
            SelectableDataFrame: Contains the project settings
        """
        # vehicle.independant_var == nan, problem for "var in independant_var" (nan not iterable)
        df = SelectableDataFrame(df).replace(np.nan, "")
        self._df = self.convert_dtypes(df)

    @property
    def df(self):
        return self._df

    @classmethod
    def generate(cls):
        """Generate template dataframe"""
        return pd.DataFrame(cls._default)

    def convert_dtypes(self, df):
        errors = []
        for col_name, col_info in self._types.items():
            try:
                if col_info["type"] in [list, tuple]:
                    df[col_name] = df[col_name].apply(
                        lambda val: col_info["type"](
                            [
                                col_info["subtype"](subval)
                                for subval in (
                                    val.replace(" ", "").split(",") if val else []
                                )
                            ]
                        )
                    )
                elif col_info["type"] is bool:
                    df[col_name] = df[col_name].apply(
                        lambda val: bool(strtobool(str(val)))
                    )
                else:
                    df[col_name] = df[col_name].apply(col_info["type"])
            except KeyError as e:
                errors.append(f"Missing column {e}")
            except ValueError as e:
                errors.append(f"Wrong data types, please correct {self.filename}")
        if errors:
            raise ValidationError(errors)
        return df

    def __contains__(self, label):
        return label in self.df.label

    def __getitem__(self, label) -> pd.Series:
        return self.df.select(**{"label": label})

    def __iter__(self):
        return (row for _, row in self.df.iterrows())

    def select(self, **selector) -> SelectableDataFrame:
        df = super().select(**selector)
        if df.empty:
            raise ValueError(f"Empty selection for {self.filename}: {selector}")
        return df


class Datasets(SubSetting):
    sheet_name: ClassVar[str] = "datasets"
    _default: ClassVar[dict] = {
        "label": ["hplc", "tissue_weight", "behavior"],
        "measurement_columns": ["compound, region", "region", "measure"],
    }
    _types: ClassVar[dict] = {
        "label": {"type": str},
        "measurement_columns": {"type": list, "subtype": str},
    }


class Experiments(SubSetting):
    sheet_name: ClassVar[str] = "experiments"
    _default: ClassVar[dict] = {
        "label": ["agonist_antagonist"],
        "group_ids": ["1, 5, 3, 4"],
        "independant_variables": ["TCB2, MDL"],
        "paired": [False],
        "parametric": [True],
    }
    _types: ClassVar[dict] = {
        "label": {"type": str},
        "group_ids": {"type": list, "subtype": int},
        "independant_variables": {"type": list, "subtype": str},
        "paired": {"type": bool},
        "parametric": {"type": bool},
    }

    @property
    def experiments(self):
        return list(self.df.label)

    def select(self, **selector):
        if selector.get("label", selector.get("experiment")) in ["default", None]:
            return self.get_default_experiment()
        return super().select(**selector)

    def select_one(self, **selector):
        return self.select(**selector).iloc[0,:]
    
    def get_default_experiment(self):
        return pd.DataFrame(
            [
                dict(
                    independant_variables=["group_id"],
                    group_column="group_id",
                    group_ids=None,
                    paired=False,
                    parametric=True,
                    label="default",
                )
            ]
        )


class Groups(SubSetting):  # TODO: generalize to Groups?
    sheet_name: ClassVar[str] = "groups"
    _default: ClassVar[dict] = {
        "group_id": [1, 5, 3, 4],
        "group_name": ["vehicles", "MDL", "TCB2", "TCB2+MDL"],
        "independant_variables": ["", "MDL", "TCB2", "TCB2, MDL"],
        "subject_ids": [
            "2, 5, 7, 9, 11, 17, 20, 28, 32, 59, 67",
            "13, 14, 15, 26, 29, 34, 42, 48, 63, 65",
            "23, 24, 31, 36, 38, 40, 44, 50, 51, 57",
            "21, 25, 35, 41, 45, 49, 52, 58, 61, 66, 69",
        ],
    }
    _types: ClassVar[dict] = {
        "group_id": {"type": int},
        "group_name": {"type": str},
        "independant_variables": {
            "type": tuple,
            "subtype": str,
        },  # tuple because list unhashable in pandas
        "subject_ids": {"type": list, "subtype": int},
    }
    control_group: ClassVar[list] = "vehicles"

    @property
    def treatments(self):
        return list(self.df.group_name)

    @property
    def subject_ids(self):
        subject_ids_type = self._types["subject_ids"]["subtype"]
        return (
            self.df["subject_ids"]
            .explode("subject_ids")
            .astype(subject_ids_type)
            .to_list()
        )

    def extend_dataset(self, dataset):
        data = self.df.explode("subject_ids")
        data["subject_id"] = data.subject_ids.astype(int)
        return SelectableDataFrame(data.extend(dataset))


class Palette(SubSetting):
    sheet_name: ClassVar[str] = "palette"
    _default: ClassVar[dict] = {
        "group_id": [1, 5, 3, 4],
        "color": ["white", "pink", "orange", "red"],
        "significance_symbol": ["*", "", "$", ""],
    }
    _types: ClassVar[dict] = {
        "group_id": {"type": int},
        "color": {"type": str},
        "significance_symbol": {"type": str},
    }

    def get_significance_palette(self):
        return self._get_palette("significance")

    def get_color_palette(self):
        return self._get_palette("color")

    def _get_palette(self, palette_type):
        return {row.group_id: row[palette_type] for _, row in self.df.iterrows()}

    def __contains__(self, value):
        return value in self.df.group_id

    def __getitem__(self, group_id) -> pd.Series:
        return self.df.select(**{"group_id": group_id})


class Statistics(SubSetting):
    sheet_name: ClassVar[str] = "statistics"
    _default: ClassVar[dict] = {
        "p_value_threshold": [0.05],
        "max_outliers": [2],
    }
    _types: ClassVar[dict] = {
        "p_value_threshold": {"type": float},
        "max_outliers": {"type": int},
    }

    def __init__(self, df: pd.DataFrame):
        super().__init__(df)
        self.p_value_threshold = self.df.p_value_threshold.iloc[0]
        self.max_outliers = self.df.max_outliers.iloc[0]
