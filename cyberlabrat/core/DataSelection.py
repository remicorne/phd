from dataclasses import dataclass, field
from typing import get_args
import pandas as pd
from cyberlabrat.core.FileSystem import FileSystem
from cyberlabrat.core.Metadata import (
    ProjectInformation,
    TreatmentInformation,
    ExperimentInformation,
)
from cyberlabrat.core.HPLC import HPLC, Outliers
from cyberlabrat.core.Statistics import QuantitativeStatistic
from cyberlabrat.core.utils import is_array_like
import matplotlib.pyplot as plt
import seaborn as sns
from cyberlabrat.core.questions import input_list, yes_or_no
from cyberlabrat.core.Constants import REGION_CLASSES, COMPOUND_CLASSES


def convert_parameter_to_list(parameter):
    if parameter is None:
        raise ValueError("None should be handled specifically")
    if isinstance(parameter, str):
        return parameter.replace(", ", ",").replace(" ,", "").split(",")
    elif is_array_like(parameter):
        return list(parameter)
    return [parameter]


def whisker_plot(data, grouping):
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
    ax.set_title(grouping)
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

    project: str | list | tuple = field(kw_only=True)
    experiment: str = field(kw_only=True, default=None)
    treatment: str = field(kw_only=True, default=None)
    compound: str | list | tuple = field(kw_only=True, default=None)
    region: str | list | tuple = field(kw_only=True, default=None)
    remove_outliers: str | bool = field(kw_only=True, default=None)
    data_source: str | bool = field(kw_only=True, default="hplc")
    p_value_threshold: float = field(kw_only=True, default=None)
    pool: str = field(kw_only=True, default=None)

    def __post_init__(self):
        if self.project not in FileSystem.list_projects():
            raise ValueError(f"Unknown project {self.project}")

        self.project_information = ProjectInformation(self.project)
        self.treatment_information = TreatmentInformation(self.project)
        self.experiment_information = ExperimentInformation(self.project)

        self.p_value_threshold = (
            self.p_value_threshold or self.project_information.p_value_threshold
        )

        self.data = HPLC(self.project).full_df
        self.compound = COMPOUND_CLASSES.get(self.compound, self.compound)
        self.region = REGION_CLASSES.get(self.region, self.region)

        self.experiment_options = self.experiment_information.experiments
        self.treatment_options = self.treatment_information.label.unique()
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

        self.data = self.data.select(compound=self.compound, region=self.region)

        if self.experiment:
            self.experiment_information = self.experiment_information.select(
                label=self.experiment
            )

            self.data = self.data.select(experiment=self.experiment)

        if self.treatment:
            self.data = self.data.select(treatment=self.treatment)
            self.treatments = (
                [self.treatment] if isinstance(self.treatment, str) else self.treatment
            )
        else:
            self.treatments = [
                label
                for label in self.treatment_information.label
                if label in self.data.treatment.unique()
            ]

        if self.remove_outliers == "eliminated":
            self.process_outliers()
            self.data = self.data.select(outlier_status=["normal", "kept"])
        elif self.remove_outliers == "calculated":
            self.data = self.data.select(
                is_outlier=lambda x: x != True
            )  # nan considered not outlier

        self.data = self.data.select(value="notna")
        if self.pool: 
            self.data[self.pool] = "all"
            

    def process_parameter(self, name, options):
        parameter = getattr(self, name)
        parameter_to_list = convert_parameter_to_list(parameter)
        unknown_params = set(parameter_to_list) - set(options)
        if unknown_params:
            raise ValueError(f"Invalid parameter(s) for {name}: {unknown_params}")
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

    statistics_pipeline: list = field(kw_only=True, default=None)

    def __post_init__(self):
        if self.experiment == "weights":
            self.compound = "weight"
            
        super().__post_init__()
        if self.experiment is not None and self.pool != "treatment":
            self.statistics, self.statistics_table = (
                QuantitativeStatistic.calculate_from_selection(
                    self.data,
                    (
                        [self.experiment_information]
                        if self.experiment == 1
                        else self.experiment_information
                    ),
                    self.p_value_threshold,
                )
            )
            if len(self.statistics) == 1:
                self.statistic = self.statistics[0]
        else:
            self.statistics, self.statistics_table, self.statistic = [], pd.DataFrame(), None
