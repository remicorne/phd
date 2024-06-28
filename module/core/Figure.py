from dataclasses import dataclass, field
from module.core.Cacheable import Cacheable
from module.core.Dataset import SelectableDataFrame
from module.core.Statistics import Statistics, QuantitativeStatistic
from module.core.Outliers import Outliers
from module.core.Metadata import (
    ProjectInformation,
    ExperimentInformation,
    TreatmentInformation,
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


@dataclass
class Figure(Cacheable):
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
    compound: str|list = field(kw_only=True, default=None)
    region: str|list = field(kw_only=True, default=None)
    experiment: str = field(kw_only=True, default=None)
    p_value_threshold: float = field(kw_only=True, default=None)
    handle_outliers: float = field(kw_only=True, default=True)
    remove_outliers: str = field(kw_only=True, default="eliminated")
    custom_params: dict = field(kw_only=True, default_factory=dict)
    extension: ClassVar[str] = "png"

    def __post_init__(self):
        self.treatment_information = TreatmentInformation(self.project)
        self.handle_parameter_logic()
        self.data = self.get_data()
        self.p_value_threshold = (
                self.p_value_threshold
                or ProjectInformation(self.project).p_value_threshold
            )
        self.define_filename()
        super().__post_init__()
        
            

    def define_filename(self):
        self.filename = f"{self.compound} in {self.region if self.region else 'all regions'}"
        
        
    def handle_parameter_logic(self):
        if self.remove_outliers not in ["eliminated", "calculated", False]:
            raise ValueError("remove_outliers must be 'eliminated', 'calculated', or False")
        if self.experiment:
            self.experiment_information = ExperimentInformation(self.project).select(experiment=self.experiment)
            self.treatments = self.experiment_information.treatments
        else:
            self.treatments = self.treatment_information.df.treatment.to_list()
            
        if not (self.compound or self.region):
            raise ValueError(
                "Must specify either compound or region to generate histogram"
            )
        self.region = self.region[0] if isinstance(self.region, list) and len(self.region) == 1 else self.region
        self.compound = self.compound[0] if isinstance(self.compound, list) and len(self.compound) == 1 else self.compound
        for param in ["compound", "region"]:
            value = getattr(self, param)
            if isinstance(param, list) and len(param) == 1:
                setattr(self, param, value[0])
        
        if isinstance(self.compound, list) and isinstance(self.region, list):
            raise ValueError(
                "Cannot do summary for multiple compounds and regions at the same time"
            )
            
    def get_data(self):
        data = HPLC(self.project).extend(TreatmentInformation(self.project)).extend(Outliers(self.project))
        data = data.select(treatment=self.treatments)
        if self.compound:
            data = data.select(compound=self.compound)
        if self.region:
            data = data.select(region=self.region)
        return data.select(nan=False)

    def generate(self):
        if self.remove_outliers == "eliminated":
            self.handle_outlier_selection()
            self.data = self.get_data()
            self.data = self.data.select(outlier_status=["normal", "kept"])
        elif self.remove_outliers == "calculated":
            self.data = self.data.select(is_outlier=False)
        
        
    def handle_outlier_selection(self):
        raise NotImplementedError("Figure must handle outlier selection")
    
    def setup_plotter_parameters(self):
        raise NotImplementedError("Figure must handle setting up plotter parameters")

    def plot(self):
        raise NotImplementedError("Figure must handle plotting")

    def save(self):
        def target():
            self.fig.savefig(self.filepath)
            print(f"SAVED {self.filepath}")
            filepath_no_extension, _ = os.path.splitext(self.filepath)
            self.fig.savefig(f"{filepath_no_extension}.svg")
            print(f"SAVED {filepath_no_extension}.svg")

        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(target)

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

    figure_type: ClassVar[str] = "histogram"
    
    def handle_parameter_logic(self):
        super().handle_parameter_logic()             
        if self.region is None or isinstance(self.region, list):
            self.summary_type = "compound"
            self.to_plot = "region"
        elif self.compound is None or isinstance(self.compound, list):
            self.summary_type = "region"
            self.to_plot = "compound"
        self.is_summary = hasattr(self, "summary_type")
        
    def setup_plotter_parameters(self):
        self.hue = self.swarm_hue = "treatment"
        self.palette = self.experiment_information.palette if self.experiment else self.treatment_information.palette
        self.hue_order = self.treatments
        self.title = (
            f"{self.compound or 'all compounds'} in {self.region or 'all regions'}"
        )
        self.ylabel = "" if "/" in self.compound else "ng/mg of tissue"
        if self.is_summary:
            self.ylabel = f"{self.__getattribute__(self.summary_type)} {self.ylabel} +/-98CI"
            self.order = [
                item
                for item in COMPOUNDS_AND_REGIONS[self.to_plot]
                if item in self.data[self.to_plot].unique()
            ]
            self.title = ""

    def generate(self):
        super().generate()
        self.setup_plotter_parameters()
        self.plot()
        if self.experiment:
            self.label_summary_stats() if self.is_summary else self.label_histogram_stats()
    
    def handle_outlier_selection(self):
        if self.handle_outliers:
            for subselection in [subselection for _, subselection in self.data.groupby(by=["compound", "region"])] if self.is_summary else [self.data]:
                if not subselection.select(outlier_status="suspected").empty:
                    finished = False
                    while not finished:
                        title = subselection.compound.unique()[0] + ' in ' + subselection.region.unique()[0]
                        base_swarm_palette = {"normal": "green", "eliminated": "red", "kept": "blue"}
                        extra_colors = ["red", "orange", "yellow", "pink", "purple", "blue"]
                        subselection["updated_outlier_status"] = subselection.apply(lambda row: row.outlier_status if row.outlier_status != "suspected" else f"suspected (mouse_id={row.mouse_id})", axis=1)
                        swarm_palette = {**base_swarm_palette, **{hue: extra_colors.pop(0) for hue in subselection.updated_outlier_status.unique() if "suspected" in hue}}
                        self.plot_histogram(custom_params={"data": subselection, "title": title, "swarm_hue": "updated_outlier_status", "swarm_palette": swarm_palette, "legend": True, "size": 10})
                        plt.show()
                        eliminated = input_list(f"Select outliers for {title} ({len(subselection.select(is_outlier=True))}): input mouse_ids to eliminated or write 'none'")
                        def label_eliminated(row):
                            if row.is_outlier:
                                return "eliminated" if str(row.mouse_id) in eliminated else "kept"
                            return row.outlier_status
                        subselection["updated_outlier_status"] = subselection.apply(label_eliminated, axis=1)
                        self.plot_histogram(custom_params={"data": subselection, "title": title, "swarm_hue": "updated_outlier_status", "swarm_palette": base_swarm_palette, "legend": True, "size": 10})                    
                        plt.show()
                        finished = yes_or_no("Confirm selection?")
                    Outliers(self.project).update(subselection)


    def plot(self, custom_params=dict()):
        custom_params = {**self.custom_params, **custom_params}
        custom_params["palette"] = {**self.palette, **self.custom_params.get("palette", {}), **custom_params.get("palette", {})}
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
        ax = sns.swarmplot(
            data=custom_params.get("data", self.data),
            x="treatment",
            y="value",
            hue=custom_params.get("swarm_hue", self.swarm_hue),
            size=custom_params.get("size", 5),
            palette=custom_params.get("swarm_palette", custom_params.get("palette", self.palette)),
            legend=custom_params.get("legend", False),
            order=custom_params.get("hue_order", self.hue_order),
            edgecolor=custom_params.get("edgecolor", "k"),
            linewidth=custom_params.get("linewidth", 1),
            linestyle=custom_params.get("linestyle", "-"),
            dodge=custom_params.get("dodge", False),
        )

        ax.tick_params(labelsize=custom_params.get("labelsize", 24))
        ax.set_ylabel(custom_params.get("ylabel", self.ylabel), fontsize=custom_params.get("ylabel_fontsize", 24))
        ax.set_xlabel(" ", fontsize=custom_params.get("xlabel_fontsize", 20))  # treatments
        ax.set_title(
            custom_params.get("title", self.title), y=custom_params.get("y", 1.04), fontsize=custom_params.get("fontsize", 34)
        )  # '+/- 68%CI'
        sns.despine(left=False)

    def plot_summary(self, custom_params):
        self.fig, self.ax = plt.subplots(figsize=(20, 10))
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
        )
        self.ax.tick_params(labelsize=16)
        self.ax.set_ylabel(self.ylabel, fontsize=24)
        self.ax.yaxis.set_label_coords(-0.035, 0.5)
        self.ax.set_xlabel(" ", fontsize=20)  # remove x title
        self.ax.set_title(self.title, y=1.04, fontsize=34)
        self.ax.legend(loc="upper right")  # , bbox_to_anchor=(0.1, 1))
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
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
        if self._statistics.is_significant:
            pairs, p_values = self._statistics.significant_pairs
            annotator = Annotator(
                self.ax, pairs, data=self.data, x="treatment", y="value", order=self.hue_order
            )
            annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
            annotator.set_pvalues_and_annotate(p_values)
            self.statistics = self._statistics.results 
            

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
                            metadata=x,
                        )
                        for x in self.order
                    ]
                )
            self.statistics = SelectableDataFrame(pd.concat([s.results for s in self._statistics])) 
        for statistics in self._statistics:
            if statistics.is_significant:
                significant_vs_control = set(
                    flatten(
                        [
                            pair
                            for pair in statistics.significant_pairs[0]
                            if self.experiment_information.control_treatment in pair
                        ]
                    )
                )
                if significant_vs_control:
                    significant_vs_control.remove(
                        self.experiment_information.control_treatment
                    )
                    for hue in significant_vs_control:
                        x_index = self.order.index(
                            statistics.metadata
                        )  # Statistics metadata stores grouping info
                        hue_index = self.hue_order.index(hue)
                        bar = self.ax.patches[hue_index * len(self.order) + x_index]
                        self.ax.text(
                            bar.get_x() + bar.get_width() / 2,
                            bar.get_height() + 0.06,
                            "*",
                            ha="center",
                            va="bottom",
                            fontsize=14,
                        )
                
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
        grouped = self.data.groupby(['region', 'compound', 'treatment']).agg(
        mean_value=('value', 'mean'),
        sem_value=('value', lambda x: np.std(x, ddof=1) / np.sqrt(len(x)))
        ).reset_index()

        # Combine mean and SEM into a single string
        grouped['mean ± SEM'] = grouped.apply(
            lambda row: f"{row['mean_value']:.2f} ± {row['sem_value']:.2f}", axis=1
        )

        # Pivot the DataFrame
        pivot_df = grouped.pivot_table(
            index='region',
            columns=['compound', 'treatment'],
            values='mean ± SEM',
            aggfunc='first'
        )

        # Sort the multiindex columns
        return pivot_df.sort_index(axis=1)
# @dataclass()
# class MatricesFigure(Figure):

#     experiment: str
#     type: str
#     to_correlate: str
#     columns: list[str] = None
#     n_minimum: int = 5
#     method: str = "pearson"
#     pvalue_threshold: float = 0.05

#     def __post_init__(self):
#         if self.type not in ["region", "compound"]:
#             raise ValueError("Type must be compound or region")
#         self.experiment = Experiment(self.experiment, self.project)
#         self.variables = self.to_correlate.split(("-"))
#         self.accross = "compound" if self.type == "region" else "region"
#         self.data = self.get_data()
#         super().__post_init__()

#     def get_data(self):
#         sub_selector = {}
#         sub_selector[self.type] = self.variables
#         if self.columns:
#             sub_selector[self.accross] = self.columns
#         return self.experiment_information.df.select(sub_selector)

#     def generate(self):
#         self.matrices = self.get_matrices()
#         self.fig, self.axs = self.generate_figure()
#         self.plot_figure()

#     def get_matrices(self):
#         return Matrices(
#             data=self.get_data(),
#             group_by="experiment",
#             between=self.type,
#             variables=self.variables,
#             accross=self.accross,
#             columns=self.columns,
#             method=self.method,
#             pvalue_threshold=self.pvalue_threshold,
#             n_minimum=self.n_minimum,
#         ).matrices

#     def generate_figure(self):
#         """
#         Generic function to create subplots of correct dimentions
#         input: experimental data listed by treatment, plotter function that takes single treatment data
#         ~optional_experimental_info may be passed such that the plotter_cb may scal axis the same for instance
#         output: plotted and saved figure at experimental level
#         """
#         # determin number of treatments to corrispond to number of subplots
#         num_treatments = len(self.matrices)
#         num_cols = min(int(np.sqrt(num_treatments)), 2)  # max of 2 columns
#         num_rows = (
#             num_treatments + num_cols - 1
#         ) // num_cols  # Compute the number of rows

#         # define the base size and a scaling factor for the figure size
#         base_size = 11
#         scale_factor = 1

#         # create subplots
#         fig, axs = plt.subplots(
#             num_rows,
#             num_cols,
#             figsize=(
#                 num_cols * base_size * scale_factor,
#                 num_rows * base_size * scale_factor,
#             ),
#             constrained_layout=True,
#         )
#         # fig.tight_layout(pad=2)
#         # fig.subplots_adjust(hspace=0.4, wspace=0.4)
#         return fig, axs.flatten()

#     def plot_ax(self):
#         pass

#     def plot_figure(self):
#         [self.plot_ax(matrix, ax) for matrix, ax in zip(self.matrices, self.axs)]

#     def identifier(self):
#         return f"{self.group_by}_{self.between}_{self.variables}_{self.accross}"


# class Correlogram(MatricesFigure):

#     def plot_ax(self, matrix, ax):
#         """
#         Correlogram plotter for single correlation matrix ~ to be fed to plotExperiment()
#         input: a single element from matricies i.e. for one treatment
#         output:  ax with graph plotted
#         """

#         ax.set_title(
#             matrix.get_title(), fontsize=28, pad=20, y=1
#         )  # Adjust the y position of the title manually for square correlogram

#         sns.heatmap(
#             matrix.corr_masked,
#             vmin=-1,
#             vmax=1,
#             square=True,
#             annot=True,
#             cmap="coolwarm",
#             annot_kws={"size": 8},
#             ax=ax,
#             cbar_kws={"shrink": 0.7},  # adj color bar size
#         )
#         ax.set_xticklabels(
#             ax.get_xticklabels()
#         )  # rotation=45, horizontalalignment='right',

#         ax.set_ylabel(matrix.var1, fontsize=28)
#         ax.set_xlabel(matrix.var2, fontsize=28)


# class NetworkFigure(MatricesFigure):

#     def plot_ax(self, matrix, ax):
#         """
#         Correlogram plotter for single correlation matrix ~ to be fed to plotExperiment()
#         input: a single element from matricies i.e. for one treatment
#         output:  ax with graph plotted
#         """

#         Network(matrix).plot_ax(ax)
