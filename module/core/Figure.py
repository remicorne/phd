from dataclasses import dataclass, field
from module.core.Cacheable import Cacheable
from module.core.Dataset import SelectableDataFrame
from module.core.Statistics import Statistics, QuantitativeStatistic
from module.core.Outliers import Outliers
from module.core.Metadata import (
    ProjectInformation,
    ExperimentInformation,
    TreatmentInformation,
    Palette,
)
from module.Matrix import Matrices
from module.Network import Network
from module.core.Constants import COMPOUNDS_AND_REGIONS
from matplotlib import pyplot as plt
import seaborn as sns
from typing import ClassVar
import numpy as np
from threading import Thread
from IPython.display import Image, display
import pandas as pd
import os
import concurrent.futures
from module.core.utils import flatten, parallel_process
from module.core.HPLC import HPLC
from statannotations.Annotator import Annotator
from module.core.Statistics import Statistics
from module.core.questions import input_list, yes_or_no
from module.core.DataSelection import DataSelection

@dataclass
class Figure(Cacheable, DataSelection):
    """
    Base class for all figures. Handles loading, saving, and plotting of figures.

    Attributes:
        project (str): The name of the project.
        compound (str|list, optional): The name of the compound. Defaults to None.
        region (str|list, optional): The name of the region. Defaults to None.
        experiment (str, optional): The name of the experiment. Defaults to None.
        p_value_threshold (float, optional): The p-value threshold used for statistical analysis defaults to value set at project level.
        handle_outliers (float, optional): Whether to handle outlier selection if any present. Defaults to True.
        remove_outliers (str, optional): Whether to remove outliers. Defaults to "eliminated".
        custom_params (dict, optional): Custom parameters for the figure. Defaults to an empty dictionary.
        extension (ClassVar[str]): The file extension for the figure. Defaults to "png".
    """

    project: str
    compound: str | list = field(kw_only=True, default=None)
    region: str | list = field(kw_only=True, default=None)
    experiment: str = field(kw_only=True, default=None)
    p_value_threshold: float = field(kw_only=True, default=None)
    remove_outliers: str = field(kw_only=True, default="calculated")
    custom_params: dict = field(kw_only=True, default_factory=dict)
    extension: ClassVar[str] = "png"

    def __post_init__(self):
        if isinstance(self.compound, str):
            self.figure_type = "compound"
            self.to_plot = "region"
        elif isinstance(self.compound, str):
            self.figure_type = "region"
            self.to_plot = "compound"
        self.define_filename()
        DataSelection.__post_init__(self)
        Cacheable.__post_init__(self)

    def define_filename(self):
        self.filename = (
            f"{self.compound} in {self.region if self.region else 'all regions'}"
        )
        
    def save(self):
        self.fig.savefig(self.filepath)
        print(f"SAVED {self.filepath}")
        filepath_no_extension, _ = os.path.splitext(self.filepath)
        self.fig.savefig(f"{filepath_no_extension}.svg")
        print(f"SAVED {filepath_no_extension}.svg")

    def load(self):
        display(Image(filename=f"{self.filepath}"))

    def _repr_html_(self):
        if not hasattr(self, "fig"):
            if self.is_saved:
                self.load()
            else:
                raise Exception("Could not load figure")


@dataclass
class Histogram(Figure):
    """
    Generate a histogram of treatments. If only one compound or region is specified, a simple histogram is generated.
    If multiple compounds or regions are specified, a summary histogram is generated.
    """

    __doc__ += Figure.__doc__

    plot_swarm: bool = field(default=True)
    figure_type: ClassVar[str] = "histogram"

    def setup_plotter_parameters(self):
        self.is_summary = isinstance(self.compound, list) or isinstance(self.region, list)
        self.hue = self.swarm_hue = "treatment"
        self.palette = (
            self.experiment_information.palette
            if self.experiment
            else self.treatment_information.palette
        )
        self.significance_palette = {
            row.treatment: row.significance for row in Palette(self.project)
        }
        self.hue_order = self.treatments
        self.title = (
            f"{self.compound or 'all compounds'} in {self.region or 'all regions'}"
        )
        self.ylabel = "" if "/" in self.compound else "ng/mg of tissue"
        if self.is_summary:
            self.ylabel = (
                f"{self.__getattribute__(self.figure_type)} {self.ylabel} +/-98CI"
            )
            self.order = [
                item
                for item in COMPOUNDS_AND_REGIONS[self.to_plot]
                if item in self.data[self.to_plot].unique()
            ]
            self.title = ""

    def generate(self):
        self.setup_plotter_parameters()
        self.plot()
        if self.experiment:
            (
                self.label_summary_stats()
                if self.is_summary
                else self.label_histogram_stats()
            )

    def plot(self, custom_params=dict()):
        custom_params = {**self.custom_params, **custom_params}
        custom_params["palette"] = {
            **self.palette,
            **self.custom_params.get("palette", {}),
            **custom_params.get("palette", {}),
        }
        if self.is_summary:
            self.plot_summary(custom_params)
        else:
            self.plot_histogram(custom_params)

    def plot_histogram(self, custom_params):
        self.fig, self.ax = plt.subplots(figsize=(20, 10))
        ax = sns.barplot(
            x="treatment",
            y="value",
            data=custom_params.get("data", self.data),
            hue="treatment",
            palette=custom_params.get("palette", self.palette),
            errorbar=custom_params.get("errorbar", ("ci", 68)),
            edgecolor=custom_params.get("edgecolor", ".2"),
            errcolor=custom_params.get("errcolor", ".2"),
            capsize=custom_params.get("capsize", 0.1),
            alpha=custom_params.get("alpha", 0.8),
            order=custom_params.get("hue_order", self.hue_order),
            dodge=custom_params.get("dodge", False),
        )
        if custom_params.get("plot_swarm", self.plot_swarm):
            ax = sns.swarmplot(
                data=custom_params.get("data", self.data),
                x="treatment",
                y="value",
                hue=custom_params.get("swarm_hue", self.swarm_hue),
                size=custom_params.get("size", 5),
                palette=custom_params.get(
                    "swarm_palette", custom_params.get("palette", self.palette)
                ),
                legend=custom_params.get("legend", False),
                order=custom_params.get("hue_order", self.hue_order),
                edgecolor=custom_params.get("edgecolor", "k"),
                linewidth=custom_params.get("linewidth", 1),
                linestyle=custom_params.get("linestyle", "-"),
                dodge=custom_params.get("dodge", False),
            )

        ax.tick_params(labelsize=custom_params.get("labelsize", 24))
        ax.set_ylabel(
            custom_params.get("ylabel", self.ylabel),
            fontsize=custom_params.get("ylabel_fontsize", 24),
        )
        ax.set_xlabel(
            " ", fontsize=custom_params.get("xlabel_fontsize", 20)
        )  # treatments
        ax.set_title(
            custom_params.get("title", self.title),
            y=custom_params.get("y", 1.04),
            fontsize=custom_params.get("fontsize", 34),
        )  # '+/- 68%CI'
        sns.despine(left=False)

    def plot_summary(self, custom_params):
        fig_width = custom_params.get("fig_width", 1.48 + 2 * len(self.order))
        self.fig, self.ax = plt.subplots(figsize=(fig_width, 10))
        self.ax = sns.barplot(
            x=self.to_plot,
            y="value",
            data=self.data,
            hue=self.hue,
            palette=custom_params.get("palette", self.palette),
            errorbar=custom_params.get("errorbar", ("ci", 68)),
            edgecolor=custom_params.get("edgecolor", ".2"),
            errcolor=custom_params.get("errcolor", ".2"),
            capsize=custom_params.get("capsize", 0.1),
            alpha=custom_params.get("alpha", 0.8),
            order=custom_params.get("order", self.order),  # self.order,
            hue_order=custom_params.get("hue_order", self.hue_order),
            errwidth=custom_params.get("errwidth", 1),
            dodge=custom_params.get("dodge", True),
            width=custom_params.get("bar_width", 0.8)
        )
        self.ax.tick_params(labelsize=36)
        self.ax.set_ylabel(self.ylabel, fontsize=38, labelpad=custom_params.get("labelpad", 100))
        self.ax.yaxis.set_label_coords(custom_params.get("ylabel_x",  -0.459 / fig_width), 0.5)
        self.ax.set_xlabel(" ", fontsize=20)  # remove x title
        self.ax.set_title(self.title, y=1.04, fontsize=34)
        self.ax.legend(loc="upper right")  # , bbox_to_anchor=(0.1, 1))
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["right"].set_visible(False)
        plt.tight_layout()

    def label_histogram_stats(self):
        if not hasattr(self, "statistics"):
            self._statistics = QuantitativeStatistic(
                self.data,
                self.experiment_information.independant_variables,
                self.experiment_information.treatments,
                self.experiment_information.paired,
                self.experiment_information.parametric,
                self.p_value_threshold,
            )
            self.statistics = self._statistics.results
        if self._statistics.is_significant:
            pairs, p_values = self._statistics.significant_pairs
            annotator = Annotator(
                self.ax,
                pairs,
                data=self.data,
                x="treatment",
                y="value",
                order=self.hue_order,
            )
            annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
            annotator.set_pvalues_and_annotate(p_values)

    def label_summary_stats(self):
        if not hasattr(self, "statistics"):
            self._statistics = parallel_process(
                [
                    QuantitativeStatistic(
                        self.data.select(**{self.to_plot: x}),
                        self.experiment_information.independant_variables,
                        self.experiment_information.treatments,
                        self.experiment_information.paired,
                        self.experiment_information.parametric,
                        self.p_value_threshold,
                        delay_execution=True,
                        metadata={self.to_plot: x},
                    )
                    for x in self.order
                ]
            )
            self.statistics = SelectableDataFrame(
                pd.concat([s.results for s in self._statistics])
            )
        for statistics in self._statistics:
            if statistics.is_significant:
                # Font Scaling # HARDCODE JJB TODO - also add significance pairs!
                base_font_size = 36
                scaling_factor = 0.48
                dynamic_font_size = max(
                    base_font_size - (scaling_factor * len(self.order)), 6
                )
                for pair in statistics.significant_pairs[0]:
                    for i, (treatment, symbol) in enumerate(
                        self.significance_palette.items()
                    ):
                        if treatment in pair:
                            hue = pair[
                                pair.index(treatment) - 1
                            ]  # work because only two elements 0 -> -1, 1 -> 0
                            x_index = self.order.index(
                                statistics.metadata[self.to_plot]
                            )  # Statistics metadata stores grouping info
                            hue_index = self.hue_order.index(hue)
                            bar = self.ax.patches[hue_index * len(self.order) + x_index]
                            self.ax.text(
                                bar.get_x() + bar.get_width() / 2,
                                (bar.get_height() * (1.3 + i / 5)),
                                symbol,
                                ha="center",
                                va="bottom",
                                fontsize=40, #dynamic_font_size
                            )
                            break

    def set(self, **kwargs):
        # kwargs["palette"] = {**self.palette, **kwargs.get("palette", {})}
        self.custom_params = kwargs
        self.initialize()
        
        



@dataclass
class Table(Figure):

    def setup_plotter_parameters(self):
        pass

    def plot(self):
        pass

    def generate(self):
        super().generate()
        grouped = (
            self.data.groupby(["region", "compound", "treatment"])
            .agg(
                mean_value=("value", "mean"),
                sem_value=("value", lambda x: np.std(x, ddof=1) / np.sqrt(len(x))),
            )
            .reset_index()
        )

        # Combine mean and SEM into a single string
        grouped["mean ± SEM"] = grouped.apply(
            lambda row: f"{row['mean_value']:.2f} ± {row['sem_value']:.2f}", axis=1
        )

        # Pivot the DataFrame
        pivot_df = grouped.pivot_table(
            index="region",
            columns=["compound", "treatment"],
            values="mean ± SEM",
            aggfunc="first",
        )

        # Sort the multiindex columns
        return pivot_df.sort_index(axis=1)
