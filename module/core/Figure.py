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


@dataclass
class Figure(Cacheable):

    def __post_init__(self):
        self.filepath = f"{self.location}/{self.identifier}"
        self.png_path = f"{self.filepath}.png"
        self.svg_path = f"{self.filepath}.svg"
        super().__post_init__()

    
    def save(self):
        def target():
            self.fig.savefig(f"{self.filepath}.svg")  # dpi also?
            print(f"SAVED {self.filepath}.svg")
            self.fig.savefig(f"{self.filepath}.png", dpi=self.fig.dpi)
            print(f"SAVED {self.filepath}.png")

        Thread(target=target).start()

    def load(self):
        display(Image(filename=f"{self.filepath}.png"))


@dataclass
class RegionCompoundParameters:
    compound: str | list
    region: str | list

    def __post_init__(self):
        if isinstance(self.compound, str):
            pass


@dataclass
class MatrixParameters:

    type: str
    var1: str
    columns: list[str]
    var2: str = None

    def __post_init__(self):
        self.var2 = self.var2 or self.var1


@dataclass
class FigureParameters:
    compound: str | list
    region: str | list
    experiment: str


@dataclass
class HistogramParameters(FigureParameters):
    compound: str
    region: str
    experiment: str


@dataclass
class SummaryHistogramParameters(FigureParameters):
    compound: list
    region: str
    experiment: str


@dataclass()
class Parameters:
    compound: list
    region: str
    experiment: str

    def generate(self):
        pass  # format parameters appropriately


@dataclass()
class DataSelection:

    full_hplc: pd.DataFrame
    parameters: Parameters

    def generate(self):
        self.hplc.select(self.parameters)


@dataclass()
class Figure:

    data_selection: DataSelection

    def generate(self):
        pass  # sns something with dataselection

    def indentifier(self):
        pass  # format params into indetifier string


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


@dataclass
class BarPlotter:

    data: pd.DataFrame
    x: str
    y: str
    order: list[str]
    hue: str
    palette: dict
    title: str
    ylabel: str

    def generate(self):
        fig, ax = plt.subplots(figsize=(20, 10))
        ax = self.draw_ax()
        self.refine_ax(ax)
        return fig, ax

    def draw_ax(self):
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
            "palette": self.palette,
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
        swarm_parameters = self.common_parameters()
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

@dataclass
class HistogramFigure(Cacheable):
    
    experiment: Experiment
    compound: str
    region: str
    
    def __post_init__(self):
        self.is_ratio = "/" in self.compound
        self.data = self.experiment.df.select(compound=self.compound, region=self.region),
        self.x, self.hue = "treatment"
        self.y = "value"
        self.order = self.experiment.groups,
        self.palete = self.experiment.palette,
        self.title = self.title,
        self.ylabel = " " if self.is_ratio else "ng/mg of tissue"
        
        
    def generate(self):
        self.figure = SwarmPlotter(
            self.data,
            self.x,
            self.y,
            self.order,
            self.hue,
            self.palete,
            self.title,
            self.ylabel,
            
        )
        
    @property
    def title(self):
        return f"{self.compound} in {self.region}"