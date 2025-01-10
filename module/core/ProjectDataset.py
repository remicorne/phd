import os, re, sys
from dataclasses import dataclass, field
from typing import ClassVar, Any
import pandas as pd
import numpy as np
import scipy
from itertools import chain
from tqdm import tqdm

from collections import namedtuple
from outliers import smirnov_grubbs as grubbs
from module.core.Dataset import PickleCachedDataFrame, SelectableDataFrame
from module.core.MeasuremenCharacteristics import MeasurementCharacteristics
from module.core.Constants import ConstantRegistry, ClassRegistry
from module.core.Metadata import (
    ProjectInformation,
    ExperimentInformation,
    GroupInformation,
    DatasetInformation,
    Palette,
)
from module.core.questions import input_escape
from module.core.utils import parallel_process
from module.core.Statistics import QuantitativeStatistic, QuantitativeStatisticBatch
from module.core.Figure import Histogram, SummaryHistogram
from module.core.FileSystem import FileSystem


class ProjectSelectableDataframe(SelectableDataFrame):
    """
    Shameful hack.
    TODO: Eliminate asap
    """

    def __init__(self, data=None, project=None, *args, **kwargs):
        super().__init__(data, *args, **kwargs)
        self.project = project

    @property
    def _constructor(self):
        return lambda *args, **kwargs: (
            ProjectSelectableDataframe(*args, project=self.project, **kwargs)
            if hasattr(self, "project")
            else SelectableDataFrame(*args, **kwargs)
        )

    def select(self, **selector) -> SelectableDataFrame:
        experiment = selector.pop("experiment", None)
        if experiment:
            experiment = ExperimentInformation(self.project).select(label=experiment)
            selector["group_id"] = experiment.groups
        data = SelectableDataFrame(self).select(**selector, value="notna")
        return (
            ProjectSelectableDataframe(data, project=self.project)
            if hasattr(self, "project")
            else data
        )


def label_group_outliers(df__test__p_value_threshold):
    # if standar variation is 0, we can't calculate outliers
    df, test, p_value_threshold = df__test__p_value_threshold
    only_values = df[df.value != 0].dropna()
    df["test"] = test
    if only_values.value.count() < 3:
        df["outlier_status"] = "Not enough data"
        return df
    outlier_test = OUTLIER_TESTS[test]
    normal_values = outlier_test(only_values.value.tolist(), p_value_threshold)
    df["is_outlier"] = df.value.apply(lambda value: value not in normal_values)
    df["outlier_status"] = df.is_outlier.apply(
        lambda is_outlier: "suspected" if is_outlier else "normal"
    )
    return df


def grubbs_test(values, p_value_threshold):
    """
    Takes a list of values on which to perform the test and returns normal values
    """
    return grubbs.test(values, alpha=float(p_value_threshold))


def iqr_test(values, p_value_threshold):
    k = {0.05: 1.5, 0.01: 2, 0.001: 3}[p_value_threshold]
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    lower_bound = q1 - k * iqr
    upper_bound = q3 + k * iqr
    normal_values = [x for x in values if lower_bound <= x <= upper_bound]
    return normal_values


OUTLIER_TESTS = {"grubbs": grubbs_test, "iqr": iqr_test}


@dataclass
class Dataset(
    PickleCachedDataFrame
):  # TODO seems to me that there is a confusion between a dataset and its linked onfo (outliers, stats..)

    project: str = field(kw_only=True)
    filename: str = field(kw_only=True)  # ClassVar[str] = "base"
    selector: dict = field(kw_only=True, default_factory=dict)  # ClassVar[str] = "base"
    with_validation: str = field(kw_only=True, default=True)

    def __post_init__(self):
        self.dataset_information = (
            DatasetInformation(self.project).select(label=self.filename).iloc[0]
        )
        self.measurement_columns = self.dataset_information.measurement_columns
        self.project_information = ProjectInformation(self.project).df
        self.experiment_information = ExperimentInformation(self.project).select(
            experiment=self.dataset_information.experiments
        )
        self.mandatory_columns = [
            self.project_information.subject_column,
            "value",
        ]
        super().__post_init__()
        self.columns = self.df.columns
        self.dataset_specific_columns = set(self.columns) - set(self.mandatory_columns)
        self.group_columns = set(self.columns) - {
            self.project_information.subject_column,
            "value",
        }

        self.data = self.get_full_df()
        self.statistics = []
        self.statistics_table = []
        self.is_generic = False

    def generate(self):
        filepath = input_escape(
            f"Enter {self.filename} filepath for {self.project} project"
        )
        if filepath.endswith(".xlsx"):
            df = pd.read_excel(filepath)
        elif filepath.endswith(".pkl"):
            df = pd.read_pickle(filepath)
        else:
            raise ValueError(f"Unsupported file type: {filepath}")
        return df

    # def save(self, data): #TODO restore validation
    #     if self.with_validation:
    #         self.validate(data)
    #     super().save(data)
    #     print(f"Data saved to {self.filepath}")

    def validate(self, df):
        """
        Validate that the dataframe has the required columns and that the values
        in these columns are valid according to the ConstantRegistry.
        """
        df_columns = df.columns
        if not all(col in df_columns for col in self.mandatory_columns):
            raise ValueError(
                f"{self.mandatory_columns} columns are mandatory, modify file and retry"
            )
        valid_mouse_ids = (
            Dataset(filename="groups", project=self.project)
            .df[self.project_information.subject_column]
            .unique()
        )
        df_mouse_ids = df[self.project_information.subject_column].unique()
        invalid_mouse_ids = set(df_mouse_ids) - set(valid_mouse_ids)
        if invalid_mouse_ids:
            raise ValueError(
                f"Invalid mouse ids: {invalid_mouse_ids}, modify file and retry"
            )
        for col_name in df_columns:
            try:
                registry = ConstantRegistry.get_registry(element_type=col_name)
                unique_values = df[col_name].unique()
                invalid_values = set(unique_values) - set(registry)
                if invalid_values:
                    correction_mapper = {
                        value: registry.choose_valid_value(value)
                        for value in unique_values
                    }
                    df[col_name].apply(correction_mapper.get)
            except FileNotFoundError:
                print(
                    f"No ConstantRegistry for element type '{col_name}', skipping validation"
                )
        return df

    def calculate_outliers(self):
        project_information = ProjectInformation(self.project)
        cases = []
        for _, subset_df in self.df[
            [
                self.project_information.subject_column,
                self.project_information.group_column,
                *self.measurement_columns,
                "value",
            ]
        ].groupby(
            [
                self.project_information.group_column,
                *self.measurement_columns,
            ]
        ):
            for test in OUTLIER_TESTS:
                cases.append(
                    (
                        subset_df,
                        test,
                        project_information.p_value_threshold,
                    )
                )
        results = parallel_process(
            cases, label_group_outliers, description="Calculating outliers"
        )  # TODO check what happens with nan
        return pd.concat(results).drop(
            columns=[self.project_information.group_column, "value"]
        )

    def calculate_group_statistics(self):
        result_ls = []
        group_columns = [
            self.project_information.group_column,
            *self.measurement_columns,
        ]
        for (group_column_values), groupby_df in tqdm(
            self.data.select(value="notna", is_outlier=False).groupby(group_columns),
            desc="Calculating group statistics",
        ):

            if len(groupby_df) >= 3:
                (
                    F,
                    p,
                ) = scipy.stats.shapiro(groupby_df["value"])
                is_parametric = p > 0.05
            else:
                F, p, is_parametric = np.nan, np.nan, np.nan

            mean, std, sem, values = [
                groupby_df.value.mean(),
                groupby_df.value.std(),
                groupby_df.value.sem(),
                groupby_df.value.values,
            ]
            result_ls.append(
                [
                    *group_column_values,
                    F,
                    p,
                    is_parametric,
                    mean,
                    std,
                    sem,
                    values,
                ]
            )
        return pd.DataFrame(
            result_ls,
            columns=[
                *group_columns,
                "shapiro_F",
                "shapiro_p",
                "is_parametric",
                "mean",
                "std",
                "sem",
                "values",
            ],
        )

    def calculate_quantitative_statistics(self, p_value_threshold=None):

        p_value_threshold = (
            p_value_threshold or self.project_information.p_value_threshold
        )
        stats_batch = QuantitativeStatisticBatch()
        measurement_cols = (
            ["measurement"] if self.is_generic else self.measurement_columns
        )
        for measurement, data in self.data.groupby(measurement_cols):
            metadata = dict(zip(measurement_cols, measurement))
            for experiment in self.experiment_information.itertuples():
                metadata["experiment"] = experiment.label
                stats_batch.add(
                    data,
                    "group_name",
                    experiment,
                    metadata,
                    p_value_threshold,
                )

        self.statistics, self.statistics_table = stats_batch.compute()

        return self

    def calculate_full_quantitative_statistics(self):
        self.calculate_quantitative_statistics()
        return self.statistics_table

    def get_linked_data(self, linked_data_type):
        filepath = self.filepath.replace(
            self.filename, f"{self.filename}_{linked_data_type}"
        )
        if not os.path.isfile(filepath):
            data = getattr(self, f"calculate_{linked_data_type}")()
            data.to_pickle(filepath)
        return SelectableDataFrame(pd.read_pickle(filepath))

    @property
    def group_statistics(self):
        return self.get_linked_data("group_statistics")

    @property
    def quantitative_statistics(self):
        return self.get_linked_data("experiment_statistics")

    @property
    def outliers(self):
        return self.get_linked_data("outliers")

    def get_full_df(self):
        return self.sort_values(self.df)

    @property
    def df(
        self,
    ):  # TODO clear up with full df, also, derived datasets should be able to be projectDatasets (network df)
        data = self.load()
        data = GroupInformation(self.project).extend_dataset(
            data
        )  # TODO should be groups.pkl here
        if self.dataset_information.unit:
            data["unit"] = self.dataset_information.unit
        # data = self.sort_values(data) #TODO check usefulness
        return self.sort_values(data)

    def select(self, **selector):
        self.selector = {**self.selector, **selector}
        if "experiment" in selector:
            self.experiment_information = ExperimentInformation(self.project).select(
                label=selector.pop("experiment")
            )
            selector["group_id"] = self.experiment_information.iloc[0].groups
        if "remove_outliers" in selector:
            test, remove_outliers = next(iter(selector["remove_outliers"].items()))
            self.data = self.data.extend(self.outliers.select(test=test))
            if remove_outliers == "eliminated":
                raise NotImplementedError
            elif remove_outliers == "calculated":
                selector.update(
                    is_outlier=lambda x: x is not True
                )  # nan considered not outlier
                del selector["remove_outliers"]
        for col in filter(lambda col: col in selector, self.measurement_columns):
            if ClassRegistry.exists(element_type=col):
                registry = ClassRegistry.get_registry(element_type=col)
                if selector[col] in registry:
                    selector[col] = registry[selector[col]]
        self.data = self.data.select(**selector)
        if self.data.empty:
            raise ValueError("No data left after selection")
        return self

    def sort_values(self, df):
        for col in self.measurement_columns:
            if ConstantRegistry.exists(element_type=col):
                registry = ConstantRegistry.get_registry(element_type=col)
                order = registry.keys()
                df[col] = pd.Categorical(df[col], categories=order, ordered=True)
            else:
                print(
                    f"No ConstantRegistry for element type '{col}', skipping validation"
                )
        return df.sort_values(
            by=[
                self.project_information.group_column,
                *self.measurement_columns,
            ]
        )

    def to_generic(self):
        if not self.is_generic:
            ordered_measurement = [
                # "dataset",
                *sorted(
                    self.measurement_columns,
                    key=lambda col: len(self.data[col].unique()),
                ),
            ]

            # self.data["measurement"] = self.data[self.measurement_columns].apply(
            #     MeasurementCharacteristics, axis=1
            # )

            self.data["measurement"] = self.data[self.measurement_columns].apply(
                tuple, axis=1
            )

            self.data["measurement"] = pd.Categorical(
                self.data["measurement"],
                categories=self.data["measurement"].unique(),
                ordered=True,
            )
            self.data["dataset"] = self.filename
            self.data.measurement = self.data.measurement.astype(object)
            self.data = self.data.drop(columns=ordered_measurement, axis=1)
            self.is_generic = True
        return self

    def get_palette(self, palette_type):
        palette = {}
        for (group_id, group_name), _ in self.data.groupby(["group_id", "group_name"]):
            palette[group_name] = (
                Palette(self.project).select(group_id=group_id).iloc[0][palette_type]
            )
        return palette

    def get_units(self):
        return self.data["unit"].unique()

    def get_selection_string(self):  # TODO: develop figure params classes
        return " in ".join(
            [
                (
                    (
                        self.selector[col]
                        if isinstance(self.selector[col], str)
                        else ", ".join(self.selector[col])
                    )
                    if col in self.selector
                    else f"all {col}s"
                )
                for col in self.measurement_columns
            ]
        )

    def __repr__(self):
        return self.data.__repr__()


# class GenericProjectDataset(ProjectDataset):
#     pass


@dataclass
class MergedDatasets:

    datasets: list[Dataset]

    def __post_init__(self):
        self.datasets = [dataset.to_generic() for dataset in self.datasets]
        self.selector = {}
        self.statistics = []
        self.statistics_table = []
        self.measurement_columns = ["dataset", "measurement"]

    def select(self, **selector):
        self.selector = {**self.selector, **selector}
        for dataset in self.datasets:
            selector = {
                col: val for col, val in selector.items() if col in dataset.columns
            }
            dataset.select(**selector)
        return self

    def calculate_quantitative_statistics(self):
        for dataset in self.datasets:
            dataset.calculate_quantitative_statistics()
        self.statistics_table = pd.concat(
            [dataset.statistics_table for dataset in self.datasets]
        )
        self.statistics = list(
            chain.from_iterable([dataset.statistics for dataset in self.datasets])
        )
        return self

    @property
    def data(self):
        data = pd.concat([dataset.data for dataset in self.datasets]).reset_index(
            drop=True
        )
        data.measurement = pd.Categorical(
            data.measurement,
            categories=data.measurement.unique(),
            ordered=True,
        )
        return data

    def get_palette(self, palette_type):
        palette = {}
        for dataset in self.datasets:
            palette.update(dataset.get_palette(palette_type))
        return palette

    def get_units(self):
        units = []
        for dataset in self.datasets:
            units.extend(dataset.get_units())
        return units

    def get_selection_string(self):
        return " and ".join(
            [dataset.get_selection_string() for dataset in self.datasets]
        )
