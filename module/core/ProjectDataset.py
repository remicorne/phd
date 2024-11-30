import os, re, sys
from dataclasses import dataclass, field
from typing import ClassVar, Any
import pandas as pd
import numpy as np
import scipy
from module.core.Dataset import PickleDataset, SelectableDataFrame
from module.core.Constants import ConstantRegistry
from module.core.Metadata import (
    ProjectInformation,
    ExperimentInformation,
    GroupInformation,
    DatasetInformation,
)
from module.core.questions import yes_or_no, input_escape
from module.core.Constants import REGIONS, COMPOUNDS
from tqdm import tqdm
from outliers import smirnov_grubbs as grubbs
from module.core.utils import parallel_process
from module.core.Statistics import QuantitativeStatistic


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
        data = SelectableDataFrame(self).select(**selector)
        return (
            ProjectSelectableDataframe(data, project=self.project)
            if hasattr(self, "project")
            else data
        )


def label_group_outliers(df__test__p_value_threshold):
    # if standar variation is 0, we can't calculate outliers
    df, test, p_value_threshold = df__test__p_value_threshold
    only_values = df[df.value != 0].dropna()
    if only_values.value.count() < 3:
        df["outlier_status"] = False
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


OUTLIER_TESTS = {"grubbs": grubbs_test}


@dataclass
class ProjectDataset(PickleDataset):

    project: str = field(kw_only=True)
    filename: str = field(kw_only=True)  # ClassVar[str] = "base"
    with_validation: str = field(kw_only=True, default=True)

    def __post_init__(self):
        self.dataset_information = DatasetInformation(self.project).select(
            label=self.filename
        )
        self.project_information = ProjectInformation(self.project)
        self.experiment_information = ExperimentInformation(self.project)
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

    def save(self, data):
        if self.with_validation:
            self.validate(data)
        super().save(data)
        print(f"Data saved to {self.filepath}")

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
            ProjectDataset(filename="groups", project=self.project)
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
        cases = [
            (
                subset_df,
                project_information.outlier_test,
                project_information.p_value_threshold,
            )
            for _, subset_df in self.df.groupby(
                [
                    self.project_information.group_column,
                    *self.dataset_information.measurement_columns,
                ]
            )
        ]
        results = parallel_process(
            cases, label_group_outliers, description="Calculating outliers"
        )
        pd.concat(results).to_pickle(self.filepath.replace(".pkl", "_outliers.pkl"))

    def calculate_group_statistics(self):
        result_ls = []
        group_columns = [
            self.project_information.group_column,
            *self.dataset_information.measurement_columns,
        ]
        for (group_column_values), groupby_df in tqdm(
            self.full_df.select(value="notna", is_outlier=False).groupby(group_columns),
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
        pd.DataFrame(
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
        ).to_pickle(self.filepath.replace(".pkl", "_group_statistics.pkl"))

    def calculate_experiment_statistics(self):
        groupings = []
        data = self.full_df
        experiment_infos = (
            [
                self.experiment_information.select(experiment=experiment)
                for experiment in self.dataset_information.experiments
            ]
            if self.dataset_information.experiments
            else [
                pd.Series(
                    dict(
                        independant_variables=["group_id"],
                        group_column=self.project_information.group_column,
                        groups=None,
                        paired=False,
                        parametric=True,
                        label=self.dataset_information.label,
                    )
                )
            ]
        )

        for experiment_info in experiment_infos:
            group_selector = {
                self.project_information.group_column: experiment_info.groups
            }
            groupings.extend(
                [
                    QuantitativeStatistic(
                        data=group_data,
                        group_column=self.project_information.group_column,
                        independant_variables=experiment_info.independant_variables,
                        is_paired=experiment_info.paired,
                        is_parametric=experiment_info.parametric,
                        p_value_threshold=self.project_information.p_value_threshold,
                        delay_execution=True,
                        metadata={
                            "project": self.project,
                            "experiment": experiment_info.label,
                            **{
                                name: value
                                for name, value in zip(
                                    self.dataset_information.measurement_columns,
                                    measurement_columns,
                                )
                            },
                        },
                    )
                    for measurement_columns, group_data in tqdm(
                        data.select(**group_selector).groupby(
                            self.dataset_information.measurement_columns,
                        ),
                        desc=f"Preparing statistical groupings for {experiment_info.label}",
                    )
                ]
            )

            statistics = parallel_process(
                groupings, description="Calculating statistics"
            )

        results = []
        for statistic in statistics:
            result = statistic.results
            result["fully_significant"] = statistic.is_significant
            results.append(result)

        pd.concat(results).to_pickle(
            self.filepath.replace(".pkl", "_experiment_statistics.pkl")
        )

    def get_linked_data(self, linked_data_type):
        filepath = self.filepath.replace(
            self.filename, f"{self.filename}_{linked_data_type}"
        )
        if not os.path.isfile(filepath):
            getattr(self, f"calculate_{linked_data_type}")()
        return ProjectSelectableDataframe(pd.read_pickle(filepath))

    @property
    def group_statistics(self):
        return self.get_linked_data("group_statistics")

    @property
    def experiment_statistics(self):
        return self.get_linked_data("experiment_statistics")

    @property
    def outliers(self):
        return self.get_linked_data("outliers")

    @property
    def full_df(self):
        data = self.df.extend(self.outliers)
        if self.dataset_information.unit:
            data["unit"] = self.dataset_information.unit
        data = GroupInformation(self.project).extend_dataset(data)
        self.sort_values(data)
        return ProjectSelectableDataframe(
            data,
            self.project,
        )

    def sort_values(self, df):
        sort_columns = [
            self.project_information.group_column,
            *self.dataset_information.measurement_columns,
        ]
        for col in sort_columns:
            try:
                registry = ConstantRegistry.get_registry(element_type=col)
                order = registry.order(df[col].unique())
                df[col] = pd.Categorical(df[col], categories=order, ordered=True)
            except FileNotFoundError:
                print(
                    f"No ConstantRegistry for element type '{col}', skipping validation"
                )
        return df.sort_values(by=sort_columns)

    def to_generic_dataset(self, selector=dict()):
        df = self.full_df
        df = df.select(**selector)
        df["dataset"] = self.filename
        ordered_measurement = sorted(
            self.dataset_information.measurement_columns,
            key=lambda col: len(df[col].unique()),
        )
        df["measurement"] = df[ordered_measurement].apply(tuple, axis=1)
        df["measurement"] = pd.Categorical(
            df["measurement"], categories=df["measurement"].unique(), ordered=True
        )
        return df.drop(columns=ordered_measurement, axis=1)
