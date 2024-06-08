from dataclasses import dataclass, field
from module.core.Cacheable import Cacheable
from module.core.Experiment import Experiment
from module.Matrix import Matrices
from module.Network import Network
from matplotlib import pyplot as plt
import seaborn as sns
from typing import ClassVar
import numpy as np
from threading import Thread
from IPython.display import Image, display
import pandas as pd
import os
import concurrent.futures


@dataclass
class Figure(Cacheable):
    
    figure_type: ClassVar[str] = "histogram"
    plotter: ClassVar = None
    
    def get_parameters(self):
        raise NotImplementedError 
        
    def generate(self):
        self.fig, self.ax = self.plotter(
            *self.get_parameters()
        ).generate()
        self.label_stats()
        return self.fig
        
    def label_stats(self):
        pass

    def save(self, fig):
        def target():
            fig.savefig(self.filepath)
            print(f"SAVED {self.filepath}.png")
            filepath_no_extension, _ = os.path.splitext(self.filepath)
            fig.savefig(f"{filepath_no_extension}.svg")
            print(f"SAVED {filepath_no_extension}.svg")

        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(target)

    def load(self):
        display(Image(filename=f"{self.filepath}"))
        
    def set(self, **kwargs):
        for attibute, value in kwargs.items():
            if not hasattr(self, attibute):
                raise ValueError(f"Unknown parameter {attibute}")
            self.__setattr__(attibute, value)
        self.initialize()

    def _repr_html_(self):
        if self.is_saved:
            self.load()
        else:
            return self.fig

@dataclass
class BarPlotter:

    data: pd.DataFrame
    x: str
    y: str
    hue: str
    palette: dict
    title: str
    ylabel: str
    order: list[str]

    def generate(self):
        fig, ax = plt.subplots(figsize=(20, 10))
        ax = self.draw_ax(ax)
        self.refine_ax(ax)
        return fig, ax

    def draw_ax(self, ax):
        return sns.barplot(**self.bar_parameters())
    
    def refine_ax(self, ax):
        ax.set_ylabel(self.ylabel, fontsize=24)
        ax.set_xlabel(" ", fontsize=20)
        ax.set_title(self.title, y=1.04, fontsize=34)
        sns.despine(left=False)
        

    def common_parameters(self):
        return {
            "data": self.data,
            "x": self.x,
            "y": self.y,
            "order": self.order,
            "hue": self.hue,
            "palette": self.palette
        }

    def bar_parameters(self):
        return {
            **self.common_parameters(),
            "errorbar": ("ci", 68),
            "capsize": 0.1,
            "alpha": 0.8,
            "errcolor": ".2",
            "edgecolor": ".2",
        }


@dataclass
class SwarmPlotter(BarPlotter):
    
    def draw_ax(self, ax):
        ax = super().draw_ax(ax)
        return sns.swarmplot(**self.swarm_parameters())
    
    def refine_ax(self, ax):
        ax.tick_params(labelsize=24)
    

    def common_parameters(self):
        return {
            **super().common_parameters(),
            "dodge": False,
        }

    def swarm_parameters(self):
        return {
            **self.common_parameters(),
            "edgecolor": "k",
            "linewidth": 1,
            "linestyle": "-",
            "legend": False,
        }


@dataclass
class OutlierPlotter(SwarmPlotter):

    swarm_hue: str
    swarm_palette: list[str]

    def swarm_parameters(self):
        swarm_parameters = super().swarm_parameters()
        swarm_parameters["hue"] = self.swarm_hue
        swarm_parameters["palette"] = self.swarm_palette
        swarm_parameters["legend"] = True
        return swarm_parameters


@dataclass
class SummaryPlotter(BarPlotter):
    
    hue_order: list[str]

    def refine_ax(self, ax):
        ax.tick_params(labelsize=16)
        ax.yaxis.set_label_coords(-0.035, 0.5)
        ax.legend(loc="upper right")
        plt.tight_layout()

    def bar_parameters(self):
        return {
            **super().common_parameters(),
            "errwidth": 1,
            "hue_order": self.hue_order,
        }
        
from module.core.FullHPLC import FullHPLC
from module.core.Metadata import TreatmentInformation

from statannotations.Annotator import Annotator
from module.core.Statistics import Statistics

@dataclass
class Histogram(Figure):
    
    project: str
    experiment: str
    compound: str
    region: str
    show_outliers: bool = field(default=False)
    plotter: ClassVar[SwarmPlotter] = SwarmPlotter
    figure_type: ClassVar[str] = "histogram"
    extension: ClassVar[str] = "png"
    
    def __post_init__(self):
        self.data = FullHPLC(self.project, self.experiment).select(compound=self.compound, region=self.region, nan=False)
        self.is_ratio = "/" in self.compound
        self.x = self.hue = "treatment"
        self.y = "value"
        self.order = self.data.sort_values(by="group_id").treatment.unique().tolist()
        self.palette = TreatmentInformation(self.project).get_palette()
        self.title = f"{self.compound} in {self.region}"
        self.ylabel = " " if self.is_ratio else "ng/mg of tissue"
        self.filename = self.title
        statistics = Statistics(self.project)
        self.statistics = statistics.select(experiment=self.experiment, compound=self.compound, region=self.region)
        self.is_significant = statistics.is_signigicant(self.experiment, self.compound, self.region)
        self.significant_pairs = statistics.get_significance_pairs(self.experiment, self.compound, self.region)
        super().__post_init__()
        
    def get_parameters(self):
        return self.data,self.x,self.y,self.hue,self.palette,self.title,self.ylabel, self.order
        
    def label_stats(self):
        if self.is_significant:
            pairs, p_values = self.significant_pairs
            annotator = Annotator(self.ax, pairs, data=self.data, x=self.x, y=self.y, order=self.order)
            annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
            annotator.set_pvalues_and_annotate(p_values)

    
class OutlierHistogram(Histogram):
    
    plotter: ClassVar[OutlierPlotter] = OutlierPlotter
    
    def __post_init__(self):
        super().__post_init__()
        self.swarm_palette = {**self.palette, "current": "black", "outlier": "grey"}
        self.swarm_hue = 'treatment'
        
    def get_parameters(self):
        return *super().get_parameters(), self.swarm_hue, self.swarm_palette
        
    
@dataclass
class TreatmentOutlierFigure:
    
    data: pd.DataFrame
    compound: str
    region: str
    treatment: str
    
    def __post_init__(self):
        self.x = 'treatment'
        self.y = 'value'
        self.hue = 'treatment'
        self.palette = None
        self.title = f"Outliers for {self.compound} {self.region} {self.treatment}"
        self.ylabel = " " if "/" in self.compound else "ng/mg of tissue"
    
    def generate(self):
        return BarPlotter(
            self.data,
            self.x,
            self.y,
            self.hue,
            self.palette,
            self.title,
            self.ylabel,
        )
    
    
    

@dataclass()
class MatricesFigure(Figure):

    experiment: str
    type: str
    to_correlate: str
    columns: list[str] = None
    n_minimum: int = 5
    method: str = "pearson"
    pvalue_threshold: float = 0.05

    def __post_init__(self):
        if self.type not in ["region", "compound"]:
            raise ValueError("Type must be compound or region")
        self.experiment = Experiment(self.experiment, self.project)
        self.variables = self.to_correlate.split(("-"))
        self.accross = "compound" if self.type == "region" else "region"
        self.data = self.get_data()
        super().__post_init__()

    def get_data(self):
        sub_selector = {}
        sub_selector[self.type] = self.variables
        if self.columns:
            sub_selector[self.accross] = self.columns
        return self.experiment.df.select(sub_selector)

    def generate(self):
        self.matrices = self.get_matrices()
        self.fig, self.axs = self.generate_figure()
        self.plot_figure()

    def get_matrices(self):
        return Matrices(
            data=self.get_data(),
            group_by="experiment",
            between=self.type,
            variables=self.variables,
            accross=self.accross,
            columns=self.columns,
            method=self.method,
            pvalue_threshold=self.pvalue_threshold,
            n_minimum=self.n_minimum,
        ).matrices

    def generate_figure(self):
        """
        Generic function to create subplots of correct dimentions
        input: experimental data listed by treatment, plotter function that takes single treatment data
        ~optional_experimental_info may be passed such that the plotter_cb may scal axis the same for instance
        output: plotted and saved figure at experimental level
        """
        # determin number of treatments to corrispond to number of subplots
        num_treatments = len(self.matrices)
        num_cols = min(int(np.sqrt(num_treatments)), 2)  # max of 2 columns
        num_rows = (
            num_treatments + num_cols - 1
        ) // num_cols  # Compute the number of rows

        # define the base size and a scaling factor for the figure size
        base_size = 11
        scale_factor = 1

        # create subplots
        fig, axs = plt.subplots(
            num_rows,
            num_cols,
            figsize=(
                num_cols * base_size * scale_factor,
                num_rows * base_size * scale_factor,
            ),
            constrained_layout=True,
        )
        # fig.tight_layout(pad=2)
        # fig.subplots_adjust(hspace=0.4, wspace=0.4)
        return fig, axs.flatten()

    def plot_ax(self):
        pass

    def plot_figure(self):
        [self.plot_ax(matrix, ax) for matrix, ax in zip(self.matrices, self.axs)]

    def identifier(self):
        return f"{self.group_by}_{self.between}_{self.variables}_{self.accross}"


class Correlogram(MatricesFigure):

    def plot_ax(self, matrix, ax):
        """
        Correlogram plotter for single correlation matrix ~ to be fed to plotExperiment()
        input: a single element from matricies i.e. for one treatment
        output:  ax with graph plotted
        """

        ax.set_title(
            matrix.get_title(), fontsize=28, pad=20, y=1
        )  # Adjust the y position of the title manually for square correlogram

        sns.heatmap(
            matrix.corr_masked,
            vmin=-1,
            vmax=1,
            square=True,
            annot=True,
            cmap="coolwarm",
            annot_kws={"size": 8},
            ax=ax,
            cbar_kws={"shrink": 0.7},  # adj color bar size
        )
        ax.set_xticklabels(
            ax.get_xticklabels()
        )  # rotation=45, horizontalalignment='right',

        ax.set_ylabel(matrix.var1, fontsize=28)
        ax.set_xlabel(matrix.var2, fontsize=28)


class NetworkFigure(MatricesFigure):

    def plot_ax(self, matrix, ax):
        """
        Correlogram plotter for single correlation matrix ~ to be fed to plotExperiment()
        input: a single element from matricies i.e. for one treatment
        output:  ax with graph plotted
        """

        Network(matrix).plot_ax(ax)

