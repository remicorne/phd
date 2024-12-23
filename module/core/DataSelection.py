from dataclasses import dataclass, field
from typing import get_args
import pandas as pd
from module.core.FileSystem import FileSystem
from module.core.Metadata import (
    DatasetInformation,
    ProjectInformation,
    GroupInformation,
    ExperimentInformation,
)
from module.core.HPLC import HPLC, Outliers
from module.core.ProjectDataset import Dataset, DatasetMerger
from module.core.Constants import ConstantRegistry
from module.core.utils import is_array_like
import matplotlib.pyplot as plt
import seaborn as sns
from module.core.questions import input_list, yes_or_no


def convert_parameter_to_list(parameter):
    if parameter is None:
        raise ValueError("None should be handled specifically")
    if isinstance(parameter, str):
        return parameter.replace(", ", ",").replace(" ,", "").split(",")
    elif is_array_like(parameter):
        return list(parameter)
    return [parameter]


def whisker_plot(data, title):

    swarm_palette = {
        "normal": "green",
        "eliminated": "red",
        "kept": "blue",
    }

    EXTRA_COLORS = [
        "red",
        "orange",
        "yellow",
        "pink",
        "purple",
        "brown",
    ]

    for hue in data.palette_hue.unique():
        if "suspected" in hue:
            swarm_palette[hue] = EXTRA_COLORS.pop(0)

    ax = sns.boxplot(y="value", data=data, dodge=False)
    sns.swarmplot(
        ax=ax,
        y="value",
        data=data,
        hue="palette_hue",
        size=5,
        palette=swarm_palette,
        dodge=False,
        edgecolor="k",
        linewidth=1,
    )
    ax.set_title(title)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


def handle_outlier_selection(title, data):
    finished = False
    while not finished:
        data["palette_hue"] = data.apply(
            lambda row: (
                f"suspected (mouse_id={row.mouse_id})"
                if row.outlier_status == "suspected"
                else row.outlier_status
            ),
            axis=1,
        )
        whisker_plot(data, title)
        eliminated = input_list(
            f"Select outliers for {title}: input mouse_ids to eliminated or write 'none'"
        )

        def label_eliminated(row):
            if row.is_outlier:
                return "eliminated" if str(row.mouse_id) in eliminated else "kept"
            return row.outlier_status

        data.palette_hue = data.apply(label_eliminated, axis=1)
        whisker_plot(data, title)
        finished = yes_or_no("Confirm selection?")
    return data.palette_hue


@dataclass
class DataSelection:

    project: str | list | tuple = field(kw_only=True, default=None)
    dataset_name: str | bool = field(kw_only=True, default="hplc")
    experiment: str = field(kw_only=True, default=None)
    treatment: str = field(kw_only=True, default=None)
    compound: str | list | tuple = field(kw_only=True, default=None)
    region: str | list | tuple = field(kw_only=True, default=None)
    remove_outliers: str | bool = field(kw_only=True, default=None)
    p_value_threshold: float = field(kw_only=True, default=None)
    pool: str = field(kw_only=True, default=None)
    request: dict = field(kw_only=True, default=None)

    def __post_init__(self):
        self.project_information = ProjectInformation(self.project)
        self.group_information = GroupInformation(self.project)
        self.experiment_information = ExperimentInformation(self.project)
        self.p_value_threshold = (
            self.p_value_threshold or self.project_information.p_value_threshold
        )
        if self.request:
            self.handle_request()
        else:
            if self.project not in FileSystem.list_projects():
                raise ValueError(f"Unknown project {self.project}")

            self.dataset = Dataset(
                project=self.project, filename=self.dataset_name
            )  # TODO clarify variables (plural?)
            self.selector = {}

            self._region = self.region
            self._compound = self.compound
            self.compound = COMPOUND_CLASSES.get_many(self.compound, self.compound)
            self.region = REGION_CLASSES.get_many(self.region, self.region)

            self.experiment_options = self.experiment_information.experiments
            self.treatment_options = self.group_information.label.unique()
            self.compound_options = self.data.compound.unique()
            self.region_options = self.data.region.unique()
            self.remove_outliers_options = ["calculated", "eliminated", False]

            for name in ["experiment", "treatment", "compound", "region"]:
                options = getattr(self, name + "_options")
                if getattr(self, name) is not None:
                    processed_parameter = self.process_parameter(name, options)
                    setattr(self, name, processed_parameter)
                if is_array_like(getattr(self, name)) and len(getattr(self, name)) == 1:
                    setattr(self, name, getattr(self, name)[0])

            self.selector.update(compound=self.compound, region=self.region)

            if self.experiment:
                self.experiment_information = self.experiment_information.select(
                    label=self.experiment
                )

                self.selector.update(experiment=self.experiment)

            if self.treatment:
                self.selector.update(treatment=self.treatment)
                self.treatments = (
                    [self.treatment]
                    if isinstance(self.treatment, str)
                    else self.treatment
                )

            if self.remove_outliers == "eliminated":
                self.process_outliers()
                self.selector.update(outlier_status=["normal", "kept"])
            elif self.remove_outliers == "calculated":
                self.selector.update(
                    is_outlier=lambda x: x is not True
                )  # nan considered not outlier

            self.selector.update(value="notna")
            self.data = self.dataset.select(**self.selector)
            if self.pool:
                self.data[self.pool] = "all"
        if self.data.empty:
            raise ValueError("No data selected")

    def handle_request(self):
        self.datasets = self.request["datasets"]
        self.selector = self.request.get("selector", {})
        self.dataset_information = DatasetInformation(self.project).select(
            label=list(self.datasets.keys())
        )
        self.dataset = (
            DatasetMerger(
                [
                    Dataset(project=self.project, filename=dataset)
                    for dataset in self.datasets
                ]
            )
            if len(self.datasets) > 1
            else Dataset(project=self.project, filename=list(self.datasets.keys()))
        )
        self.dataset = self.dataset.select(**self.selector)

    def process_parameter(self, name, options):
        parameter = getattr(self, name)
        parameter_to_list = convert_parameter_to_list(parameter)
        unknown_params = set(parameter_to_list) - set(options)
        if unknown_params:
            print(f"Unknown parameter(s) for {name}: {unknown_params}")
        return parameter_to_list

    def process_outliers(self):
        for (treatment, compound, region), data in self.data.groupby(
            ["treatment", "compound", "region"]
        ):
            if "suspected" in data.outlier_status.values:
                title = f"{compound} in {region} for {treatment}"
                data.outlier_status = handle_outlier_selection(title, data)
                Outliers(self.project).update(data)


@dataclass()
class QuantitativeDataSelection(DataSelection):

    def __post_init__(self):
        if self.experiment == "weight":
            self.compound = "weight"

        super().__post_init__()
        if self.experiment is not None and self.pool != "treatment":
            self.dataset.calculate_quantitative_statistics()
            self.statistics = self.dataset.statistics
            self.statistics_table = self.dataset.statistics_table
            if len(self.dataset.statistics) == 1:
                self.statistic = self.dataset.statistics[0]
