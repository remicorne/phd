from dataclasses import dataclass, field
from operator import attrgetter
from typing import Callable
import pandas as pd
import numpy as np
import scipy
from itertools import chain
from tqdm import tqdm

from outliers import smirnov_grubbs as grubbs
from module.core.Dataset import PickleCachedDataFrame, SelectableDataFrame
from module.core.Registry import Registry, ClassRegistry
from module.core.Metadata import (
    ProjectMetadata,
)
from module.core.questions import input_escape, yes_or_no
from module.core.utils import parallel_process, is_array_like
from module.core.Statistics import QuantitativeStatisticBatch
from module.core.Ratio import Ratio
from module.core.enums import DatasetColumn, SelectorColumns


def label_group_outliers(df__test__p_value_threshold__max_outliers):
    # if standar variation is 0, we can't calculate outliers
    df, test, p_value_threshold, max_outliers = (
        df__test__p_value_threshold__max_outliers
    )
    only_values = df[df.value != 0].dropna()
    df["test"] = test
    if only_values.value.count() < 3:
        df["outlier_status"] = "Not enough data"
        return df
    outlier_test: Callable = OUTLIER_TESTS[test]
    outliers = outlier_test(only_values.value.tolist(), p_value_threshold)
    if len(outliers) > max_outliers:
        outliers = outliers[:max_outliers]
    df["is_outlier"] = df.value.apply(lambda value: value in outliers)
    df["outlier_status"] = df.is_outlier.apply(
        lambda is_outlier: "suspected" if is_outlier else "normal"
    )
    return df


def grubbs_test(values, p_value_threshold):
    """
    Takes a list of values on which to perform the test and returns normal values
    """
    return grubbs.two_sided_test_outliers(values, alpha=p_value_threshold)


def iqr_test(values, p_value_threshold):
    k = {0.05: 1.5, 0.01: 2, 0.001: 3}[p_value_threshold]
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    lower_bound = q1 - k * iqr
    upper_bound = q3 + k * iqr
    outliers = [x for x in values if x < lower_bound or x > upper_bound]
    outliers = sorted(
        outliers, key=lambda x: abs(x - (q1 if x < lower_bound else q3)), reverse=True
    )
    return outliers


OUTLIER_TESTS = {"grubbs": grubbs_test, "iqr": iqr_test}


@dataclass
class ProjectDataset(
    PickleCachedDataFrame
):  # TODO seems to me that there is a confusion between a dataset and its linked onfo (outliers, stats..)
    project: str = field(kw_only=True)
    filename: str = field(kw_only=True)  # ClassVar[str] = "base"

    # Necessary to avoir infinite recursion with self.df -> self.data before data is init
    # Only problem in debug mode because debugger wants to display self.df
    # This mostly solves the issue but the debugger will still fail between self init and data init
    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls)
        self.data = None
        return self

    def __post_init__(self):
        self.metadata = ProjectMetadata(self.project)
        self.measurement_columns = self.metadata.datasets.select_one(
            label=self.filename
        ).measurement_columns
        self.subject_column = DatasetColumn.SUBJECT_ID
        self.group_column = DatasetColumn.GROUP_ID
        super().__post_init__()
        self.data: SelectableDataFrame = self.validate(self.load())
        self.columns = self.data.columns
        self.data = self.metadata.groups.extend_dataset(self.data)
        self.selector = {}
        self.selection = {}
        self.statistics = []
        self.statistics_table = []
        self.is_generic = False

    def generate(self):
        if not yes_or_no(
            f"Initialize new dataset '{self.filename}' for '{self.project}' project?"
        ):
            raise ValueError(f"Unknwon dataset {self.filename}")
        filepath = input_escape(
            f"Enter {self.filename} filepath for {self.project} project"
        )
        if filepath.endswith(".xlsx"):
            df = pd.read_excel(filepath, keep_default_na=False)
        if filepath.endswith(".csv"):
            df = pd.read_csv(filepath, keep_default_na=False)
        elif filepath.endswith(".pkl"):
            df = pd.read_pickle(filepath)
        else:
            raise ValueError(f"Unsupported file type: {filepath}")
        num_zero_values = len(df[df["value"] == 0])
        if num_zero_values and yes_or_no(
            f"{num_zero_values} values are equal to 0, replace with nan?"
        ):
            df[DatasetColumn.VALUE] = df[DatasetColumn.VALUE].replace(0, np.nan)
        df[DatasetColumn.VALUE] = df[DatasetColumn.VALUE].replace(
            {"NA": np.nan, "": np.nan}
        )
        return self.validate(df)

    def validate(self, df: pd.DataFrame) -> SelectableDataFrame:
        """
        Validate that the dataframe has the required columns and that the values
        in these columns are valid according to the Registry.
        First step of validation is the structure of the dataset. Expect columns are [subject_id, *measurement_colums, unit]
        measurement columns = columns characterising a measurement for a subject (compo+region, cell+signal, region, behavior, etc)
        These will then be the ones you use as parameters
        """
        df_columns = df.columns
        mandatory_columns = DatasetColumn.get_mandatory_columns()
        if missing_columns := [
            col for col in mandatory_columns if col not in df_columns
        ]:
            raise ValueError(
                f"{missing_columns} columns are missing, modify file and retry"
            )
        if invalid_subject_ids := set(df[self.subject_column].unique()) - set(
            self.metadata.subject_ids
        ):
            raise ValueError(
                f"Invalid mouse ids: {invalid_subject_ids}, modify file and retry"
            )

        for measurement_col in set(df_columns) - set(mandatory_columns):
            if Registry.exists(element_type=measurement_col):
                registry = Registry.get_registry(element_type=measurement_col)
                unique_values = df[measurement_col].unique()
                invalid_values = set(unique_values) - set(registry)
                if invalid_values:
                    correction_mapper = {
                        value: registry.choose_valid_value(value)
                        for value in unique_values
                    }
                    df[measurement_col] = df[measurement_col].apply(
                        correction_mapper.get
                    )
            else:
                print(
                    f"No Registry for element type '{measurement_col}', skipping validation"
                )
        df.value = df.value.astype(float)
        return df

    def calculate_outliers(self, test: str):
        cases = []
        for _, subset_df in self.data[
            [
                self.subject_column,
                self.group_column,
                *self.measurement_columns,
                DatasetColumn.VALUE,
            ]
        ].groupby(
            [
                self.group_column,
                *self.measurement_columns,
            ]
        ):
            cases.append(
                (
                    subset_df.copy(),
                    test,
                    self.metadata.p_value_threshold,
                    self.metadata.max_outliers,
                )
            )
        results = parallel_process(
            cases,
            label_group_outliers,
            description="Calculating outliers",
            optimize=True,
        )  # TODO check what happens with nan
        self.outliers = pd.concat(results).drop(
            columns=[self.group_column, DatasetColumn.VALUE]
        )

    def calculate_group_statistics(self):
        result_ls = []
        group_columns = [
            self.group_column,
            *self.measurement_columns,
        ]
        for (group_column_values), groupby_df in tqdm(
            self.data.select(value="notna").groupby(group_columns),
            desc="Calculating group statistics",
        ):
            if len(groupby_df) >= 3:
                (
                    F,
                    p,
                ) = scipy.stats.shapiro(groupby_df[DatasetColumn.VALUE])
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
        self.group_statistics = pd.DataFrame(
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
        return self

    def calculate_quantitative_statistics(self, p_value_threshold=None):
        p_value_threshold = p_value_threshold or self.metadata.p_value_threshold
        stats_batch = QuantitativeStatisticBatch()
        measurement_cols = (
            ["measurement"] if self.is_generic else self.measurement_columns
        )
        for measurement, data in self.data.groupby(measurement_cols, observed=True):
            metadata = dict(zip(measurement_cols, measurement))
            metadata["experiment"] = self.selected_experiment.label
            stats_batch.add(
                data,
                "group_name",
                self.selected_experiment,
                metadata,
                p_value_threshold,
            )

        self.statistics, self.statistics_table = stats_batch.compute()
        return self

    def select(self, **selector):
        self.selector = {**self.selector, **selector}
        selection = {**selector}
        # Pop it so it doesnt go through classic "select" filtering
        outlier_config = selection.pop("remove_outliers", None)
        if "experiment" in selector:
            if not isinstance(selector["experiment"], str):
                raise ValueError("Experiment must be a string")
            experiment = selection.pop("experiment")
            self.selected_experiment = self.metadata.experiments.select_one(
                label=experiment
            )
            if "group_name" not in selection:
                selection["group_name"] = [
                    self.metadata.groups.select_one(group_id=group_id).group_name
                    for group_id in self.selected_experiment.group_ids
                ]
        for col in set.intersection(set(self.measurement_columns), set(selection)):
            if callable(selection[col]):
                selection[col] = list(filter(selection[col], self.data[col].unique()))
            elif ClassRegistry.exists(element_type=col):
                registry = ClassRegistry.get_registry(element_type=col)
                if selection[col] in registry:
                    selection[col] = registry[selection[col]]
                if isinstance(selection[col], str):
                    selection[col] = [selection[col]]
        self.selection = {**self.selection, **selection}
        data_with_ratios = self.build_ratios(self.data, selection)
        print(data_with_ratios)
        data_selected = data_with_ratios.select(**selection)
        self.data = self.sort_values(data_selected, selection)
        if outlier_config:
            self.data = self.remove_outliers(outlier_config)
        if self.data.empty:
            raise ValueError(
                f"No data for selection: {selection}, try different selection"
            )
        return self

    def remove_outliers(self, remove_outliers: dict[str, str]):
        test, how = next(iter(remove_outliers.items()))
        self.calculate_outliers(test)
        self.data = self.data.extend(self.outliers)
        if how == "eliminated":
            raise NotImplementedError("Manual outlier selection not implemented")
        elif how == "calculated":
            return self.data.select(is_outlier=lambda x: x is not True)

    def build_ratios(
        self, data: SelectableDataFrame, selection: dict
    ) -> SelectableDataFrame:
        ratios = []
        for col in self.measurement_columns:
            for value in selection.get(col, []):
                if Ratio.is_ratio(value):
                    ratios.append((col, value))
        if ratios:
            ratios_df = pd.concat(
                Ratio(column, value).compute(data) for column, value in ratios
            )
            data = pd.concat([data, ratios_df])
        return data

    def sort_values(self, data, categoricals):
        categories_to_create = {
            col: values
            for col, values in categoricals.items()
            if col in data and is_array_like(values) and len(values) > 1
        }
        for col, values in categories_to_create.items():
            data[col] = pd.Categorical(
                data[col], categories=values, ordered=True
            ).remove_unused_categories()  # Necessary, .loc assignment doesnt work
        return data.sort_values(by=list(categories_to_create))

    def to_generic(self):
        if not self.is_generic:
            ordered_measurement = [
                # "dataset",
                *sorted(
                    self.measurement_columns,
                    key=lambda col: len(self.data[col].unique()),
                ),
            ]

            self.data["measurement"] = self.data[self.measurement_columns].apply(
                tuple, axis=1
            )

            self.data["measurement"] = pd.Categorical(
                self.data["measurement"],
                categories=self.data["measurement"].unique(),
                ordered=True,
            )
            self.data["dataset"] = self.filename
            self.data = self.data.drop(columns=ordered_measurement, axis=1)
            self.is_generic = True
        return self

    def get_palette(self, palette_type):
        palette = {}
        for (group_id, group_name), _ in self.data.groupby(["group_id", "group_name"]):
            palette[group_name] = self.metadata.palette.select(group_id=group_id).iloc[
                0
            ][palette_type]
        return palette

    def get_units(self):
        return self.data[DatasetColumn.UNIT].unique()

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

    @property
    def df(self):
        return self.data


@dataclass
class MergedDatasets:
    datasets: list[ProjectDataset]

    def __post_init__(self):
        self.datasets = [dataset.to_generic() for dataset in self.datasets]
        self.selector = {}
        self.selection = {}
        self.statistics = []
        self.statistics_table = []
        self.measurement_columns = ["dataset", "measurement"]

    def select(self, **selector):
        for dataset in self.datasets:
            selector = {
                col: val
                for col, val in selector.items()
                if col in dataset.columns.union(list(SelectorColumns))
            }
            dataset.select(**selector)
            self.selection = {**self.selection, **dataset.selection}
            self.selector = {**self.selector, **dataset.selector}
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

    @property
    def selected_experiment(self):
        selected_experiments = [
            getattr(dataset, "selected_experiment", None) for dataset in self.datasets
        ]
        if len(set(map(attrgetter("label"), selected_experiments))) > 1:
            raise ValueError(
                f"Multiple experiments selected: {selected_experiments}, please select only one"
            )
        return selected_experiments[0]
