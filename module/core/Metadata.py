from dataclasses import dataclass, field
from module.core.utils import strtobool
from typing import ClassVar, Dict
from functools import partial
import numpy as np
import pandas as pd

from module.core.Dataset import (
    ExcelCachedDataFrame,
    CustomDataFrame,
)
from module.core.FileSystem import FileSystem
from module.core.questions import yes_or_no
from module.core.Dataset import SelectionError


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

    datasets: "Datasets" = field(init=False)
    groups: "Groups" = field(init=False)
    experiments: "Experiments" = field(init=False)
    palette: "Palette" = field(init=False)
    statistics: "Statistics" = field(init=False)

    def __post_init__(self):
        self.check_new_project()
        super().__post_init__()
        metadata = self.load()
        self.datasets: Datasets = metadata["datasets"]
        self.groups: Groups = metadata["groups"]
        self.experiments: Experiments = metadata["experiments"]
        self.palette: Palette = metadata["palette"]
        self.statistics: Statistics = metadata["statistics"]

        self.subject_ids = self.groups.subject_ids
        self.p_value_threshold = self.statistics.p_value_threshold
        self.max_outliers = self.statistics.max_outliers

    def check_new_project(self):
        if not FileSystem.project_exists(self.project):
            if not yes_or_no(
                f"Project '{self.project}' not found. Initialize new project {self.project}?"
            ):
                raise ValueError(f"Unknown project: {self.project}")

    def validate_consistency(self, metadata):
        experiment_group_ids = set()
        for experiment in metadata["experiments"]:
            experiment_group_ids.update(experiment.group_ids)

        groups_group_ids = set(metadata["groups"].df.group_id)
        palette_group_ids = set(metadata["palette"].df.group_id)

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

    def load(self) -> Dict[str, "SubSetting"]:
        metadata = {
            "datasets": Datasets(self.project),
            "groups": Groups(self.project),
            "experiments": Experiments(self.project),
            "palette": Palette(self.project),
            "statistics": Statistics(self.project),
        }
        self.validate_consistency(metadata)
        return metadata

    def initialize(self):
        """
        Same as cacheable initialize but also makes user edit as setting are user editable
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
                question = f"Correct errors in metadata.xlsx file: {e}. Press any key and ENTER when done correcting file"
                user_finished = False

    def save(self, content: Dict[str, pd.DataFrame]):
        """Save all sheet in content to multi sheet dataframe using keys as sheet names"""
        with pd.ExcelWriter(self.filepath) as writer:
            for sheet_name, df in content.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    def select(self, *, select_one=False, **kwarg):
        """
        Select a specific element from a subsetting based on its label

        Args:
            **kwarg: Keyword arguments with one key-value pair
                example: select(dataset="dataset1")

        Returns:
            pd.Series: The element from the subsetting
        """
        if len(kwarg) != 1:
            raise ValueError("Only one argument is allowed")
        subsetting, label = kwarg.popitem()
        if subsetting not in ["dataset", "group", "experiment"]:
            raise ValueError(
                "Method 'select' Only implemented for 'dataset', 'group' and 'experiment'"
            )
        subsetting: SubSetting = getattr(self, subsetting + "s")
        return subsetting.select(select_one=select_one, label=label)


def convert_iterable(values: list | tuple, iterable_type: type, value_type: type):
    values = values.replace(" ", "").split(",") if values else []
    return iterable_type([value_type(value) for value in values])


def convert_to_bool(value):
    return bool(strtobool(str(value)))


def get_converter(col_info):
    if col_info["type"] in [list, tuple]:
        converter = partial(
            convert_iterable,
            iterable_type=col_info["type"],
            value_type=col_info["subtype"],
        )
    elif col_info["type"] is bool:
        converter = convert_to_bool
    else:
        converter = col_info["type"]
    return converter


@dataclass
class SubSetting(ExcelCachedDataFrame):
    """Base class for project settings subsets.
    Handles validation and data processing for individual sheets.
    """

    project: str = field(default=None)
    filename: ClassVar[str] = ProjectMetadata.filename
    sheet_name: ClassVar[str] = None
    _default: ClassVar[dict] = None
    _types: ClassVar[dict] = None

    @classmethod
    def generate(cls):
        """Generate template dataframe"""
        return pd.DataFrame(cls._default)

    def convert_dtypes(self, df):
        errors = []
        if missing_columns := [col for col in self._types if col not in df]:
            errors.append(
                f"Sheet '{self.sheet_name}': missing columns: {missing_columns}"
            )

        for col_name in set(self._types.keys()) - set(missing_columns):
            col_types = self._types[col_name]
            values = []
            converter = get_converter(col_types)
            try:
                for value in df[col_name]:
                    value = converter(value)
                    values.append(value)
                df[col_name] = values
            except ValueError as _:
                type_hint = col_types["type"].__name__
                if "subtype" in col_types:
                    type_hint = f"{type_hint}[{col_types['subtype'].__name__}]"
                errors.append(
                    f"Sheet '{self.sheet_name}', column '{col_name}': '{value}' should be {type_hint}"
                )
        if errors:
            raise ValidationError(errors)
        return df

    def __contains__(self, label):
        return label in self.df.label

    def __getitem__(self, label) -> pd.Series:
        return self.df.select(**{"label": label})

    def __iter__(self):
        return (row for _, row in self.df.iterrows())

    def get(self, *, label) -> CustomDataFrame:
        try:
            return self.select(select_one=True, label=label)
        except SelectionError as e:
            if "'label'" in str(e):
                raise ValueError(
                    f"Unknown {self.sheet_name.rstrip('s')}: {label}, add to metadata"
                ) from e

    @property
    def df(self):
        return self.convert_dtypes(
            super(ExcelCachedDataFrame, self).df.replace(np.nan, "")
        )


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
        "label": ["dose_response", "agonist_antagonist"],
        "group_ids": ["1, 2, 3, 4", "1, 5, 3, 6"],
        "independant_variables": ["TCB2", "TCB2, MDL"],
        "paired": [False, False],
        "parametric": [True, True],
    }
    _types: ClassVar[dict] = {
        "label": {"type": str},
        "group_ids": {"type": list, "subtype": int},
        "independant_variables": {"type": tuple, "subtype": str},
        "paired": {"type": bool},
        "parametric": {"type": bool},
    }

    @property
    def experiments(self):
        return list(self.df.label)

    def _get_default_experiment(self):
        return pd.DataFrame(
            [
                dict(
                    independant_variables=["group_id"],
                    group_ids=Groups(self.project).df.group_id.values,
                    paired=False,
                    parametric=True,
                    label="default",
                )
            ]
        )

    @property
    def df(self):
        df = super().df
        if "default" in df.label.values:
            raise ValueError("Default experiment is a reserved keyword")
        return pd.concat([df, self._get_default_experiment()])

    def get_subjects(self, experiment):
        return self.select(label=experiment).group_ids


class Groups(SubSetting):  # TODO: generalize to Groups?
    sheet_name: ClassVar[str] = "groups"
    _default: ClassVar[dict] = {
        "group_id": [1, 2, 3, 4, 5, 6],
        "group_name": [
            "vehicles",
            "0,3mg/kg TCB",
            "3mg/kg TCB",
            "10mg/kg TCB",
            "0,2mg/kg MDL",
            "TCB2+MDL",
        ],
        "independant_variables": ["", "TCB2", "TCB2", "TCB2", "MDL", "TCB2, MDL"],
        "subject_ids": [
            "2, 5, 7, 9, 11, 17, 20, 28, 32, 59, 67",
            "13, 14, 15, 26, 29, 34, 42, 48, 63, 65",
            "16, 18, 19, 22, 27, 30, 37, 39, 46, 53, 55",
            "8, 10, 12, 47, 54, 56, 60, 62, 64, 68, 70",
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
        data.drop(columns=["subject_ids"], inplace=True)
        return CustomDataFrame(data.extend(dataset))


class Palette(SubSetting):
    sheet_name: ClassVar[str] = "palette"
    _default: ClassVar[dict] = {
        "group_id": [1, 2, 3, 4, 5, 6],
        "color": [
            "white",
            "lightgreen",
            "limegreen",
            "darkgreen",
            "lightgrey",
            "darkseagreen",
        ],
        "significance_symbol": ["*", "", "$", "", "", "#"],
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

    def __post_init__(self):
        super().__post_init__()
        self.p_value_threshold = self.df.p_value_threshold[0]
        self.max_outliers = self.df.max_outliers[0]
