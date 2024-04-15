from dataclasses import dataclass
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

    from_scratch = None
    
    def save(self):
        def target():
            self.fig.savefig(f"{self.filepath}.svg")  # dpi also?
            print(f"SAVED {self.filepath}.svg")
            self.fig.savefig(f"{self.filepath}.png", dpi=self.fig.dpi)
            print(f"SAVED {self.filepath}.png")

        Thread(target=target).start()

    def load(self):
        display(Image(filename=f"{self.filepath}.png"))
        
    def get(self):
        if self.from_scratch:
            return self.fig            
        return super().get()  

    @property
    def identifier(self):
        pass

    @property
    def filepath(self):
        return f"{self.location}/{self.identifier}"

    @property
    def png_filepath(self):
        return f"{self.filepath}.png"

    @property
    def svg_filepath(self):
        return f"{self.filepath}.svg"


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
        self.variables = self.to_correlate.split(('-'))
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

    def identify(self):
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





Correlogram(project, exp, to_correlate, columns)

Project(project).experiment(experiment).correlogram(to_correlate, columns)