from dataclasses import dataclass, field

import networkx as nx
from scipy.stats import norm
from scipy import stats

from module.core.Cacheable import Cacheable
from module.core.Dataset import SelectableDataFrame, ExcelDataset
from module.core.Statistics import QuantitativeStatistic
from module.core.Metadata import (
    Palette,
)
from module.core.Matrix import Matrix
from module.core.Matrix import Network as NetworkModel
from module.core.Constants import COMPOUNDS_AND_REGIONS, REGIONS
from matplotlib import pyplot as plt
import seaborn as sns
from typing import ClassVar
import numpy as np
from IPython.display import Image, display
import pandas as pd
import os
from module.core.utils import parallel_process
from statannotations.Annotator import Annotator
from module.core.DataSelection import DataSelection, QuantitativeDataSelection
from module.core.Constants import REGION_CLASSES, COMPOUND_CLASSES


class Figure:

    def __new__(self, data_handler):

        @dataclass
        class Figure(Cacheable, data_handler):
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

            custom_params: dict = field(kw_only=True, default_factory=dict)
            extension: ClassVar[str] = "png"

            def __post_init__(self):
                data_handler.__post_init__(self)
                self.setup()
                self.define_filename()
                Cacheable.__post_init__(self)

            def setup(self):
                self.compound_or_region = "compound" if self.is_compound() else "region"
                self.to_plot = "region" if self.is_compound() else "compound"
                self.order = COMPOUNDS_AND_REGIONS[self.to_plot].order(
                    self.data[self.to_plot].unique()
                )

            def is_compound(self):
                return isinstance(self.compound, str)

            def define_filename(self):
                if self.is_compound() == True:
                    if self.region in REGION_CLASSES.values():
                        self.filename = f"{self.compound} in {REGION_CLASSES.reversed.get(tuple(self.region), self.region)}"
                    else:
                        self.filename = f"{self.compound} in {self.region if self.region else 'all regions'}"
                else:
                    if self.compound in COMPOUND_CLASSES.values():
                        self.filename = f"{self.region} in {COMPOUND_CLASSES.reversed.get(tuple(self.compound), self.compound)}"
                    else:
                        self.filename = f"{self.region} in {self.compound if self.compound else 'all compounds'}"

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

        return Figure


@dataclass
class Histogram(Figure(QuantitativeDataSelection)):
    """
    Generate a histogram of treatments. If only one compound or region is specified, a simple histogram is generated.
    If multiple compounds or regions are specified, a summary histogram is generated.
    """

    plot_swarm: bool = field(default=True)
    figure_type: str = "histogram"

    def setup_plotter_parameters(self):
        multiple_regions = isinstance(self.region, list) or self.region is None
        multiple_compounds = isinstance(self.compound, list) or self.compound is None
        self.swarm_hue = self.hue = self.x = "treatment"

        self.palette = (
            self.experiment_information.palette
            if self.experiment
            else Palette(self.project).dict
        )
        self.significance_palette = {
            row.treatment: row.significance for row in Palette(self.project)
        }
        self.hue_order = self.treatments

        if self.pool == "treatment":
            self.x = self.to_plot
            self.hue = self.compound_or_region
            self.hue_order = self.order
            self.palette = self.significance_palette = None
            self.is_summary = multiple_compounds and multiple_regions
        else:
            self.is_summary = (
                isinstance(self.compound, list)
                or isinstance(self.region, list)
                or self.compound is None
                or self.region is None
            )
            
        self.title = (
            f"{self.compound or 'all compounds'} in {self.region or 'all regions'}"
        )
        self.ylabel = "" if "/" in self.compound else "ng/mg of tissue"
        if self.is_summary:
            self.ylabel = f"{self.compound_or_region} {self.ylabel} +/-98CI"
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
        custom_params["palette"] = (
            {
                **self.palette,
                **self.custom_params.get("palette", {}),
                **custom_params.get("palette", {}),
            }
            if not self.pool == "treatment"
            else None
        )
        if self.is_summary:
            self.plot_summary(custom_params)
        else:
            self.plot_histogram(custom_params)

    def plot_histogram(self, custom_params):
        self.fig, self.ax = plt.subplots(figsize=(20, 10))
        ax = sns.barplot(
            data=custom_params.get("data", self.data),
            x=self.x,
            y="value",
            hue=self.hue,
            palette=custom_params.get("palette", self.palette),
            errorbar=custom_params.get(
                "errorbar", "sd"
            ),  # ("ci", 68) remi :) would be nice to also have SEM intergrated as an option - it is not inbuild like ci or sd in seabourn
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
                x=self.x,
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
        )
        sns.despine(left=False)

    def plot_summary(self, custom_params):
        fig_width = custom_params.get("fig_width", 1.48 + 2 * len(self.order))
        self.fig, self.ax = plt.subplots(figsize=(fig_width, 10))
        self.ax = sns.barplot(
            data=self.data,
            x=self.to_plot,
            y="value",
            hue=self.hue,
            palette=custom_params.get("palette", self.palette),
            errorbar=custom_params.get("errorbar", "sd"),
            edgecolor=custom_params.get("edgecolor", ".2"),
            errcolor=custom_params.get("errcolor", ".2"),
            capsize=custom_params.get("capsize", 0.1),
            alpha=custom_params.get("alpha", 0.8),
            order=custom_params.get("order", self.order),  # self.order,
            hue_order=custom_params.get("hue_order", self.hue_order),
            errwidth=custom_params.get("errwidth", 1),
            dodge=custom_params.get("dodge", True),
            width=custom_params.get("bar_width", 0.8),
        )
        self.ax.tick_params(labelsize=40)
        self.ax.set_ylabel(
            self.ylabel, fontsize=44, labelpad=custom_params.get("labelpad", 100)
        )
        self.ax.yaxis.set_label_coords(
            custom_params.get("ylabel_x", -0.459 / fig_width), 0.5
        )
        self.ax.set_xlabel(" ", fontsize=20)  # remove x title
        self.ax.set_title(self.title, y=1.04, fontsize=34)
        self.ax.legend(loc="upper right")  # , bbox_to_anchor=(0.1, 1))
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["right"].set_visible(False)
        plt.tight_layout()

    def label_histogram_stats(self):
        if self.statistic and self.statistic.is_significant:
            pairs, p_values = self.statistic.significant_pairs
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
        for statistics in self.statistics:
            if statistics.is_significant:
                # Font Scaling # HARDCODE JJB TODO - also add significance pairs!
                base_font_size = 48
                scaling_factor = 0.2
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
                                fontsize=dynamic_font_size,
                            )
                            break

    def set(self, **kwargs):
        # kwargs["palette"] = {**self.palette, **kwargs.get("palette", {})}
        self.custom_params = kwargs
        self.initialize()


@dataclass
class MatricesFigure(Figure(DataSelection)):

    columns: list[str] = field(default=None)
    n_minimum: float = field(default=5)
    method: float = field(default="pearson")

    def __post_init__(self):
        if self.compound and "-" in self.compound:
            self.compound = self.compound.split("-")
        if self.region and "-" in self.region:
            self.region = self.region.split("-")
        super().__post_init__()

    def setup(self):
        super().setup()
        c_or_r = getattr(self, self.compound_or_region)
        self.var1 = c_or_r[0] if isinstance(c_or_r, list) else c_or_r
        self.var2 = c_or_r[-1] if isinstance(c_or_r, list) else c_or_r
        self.is_square = self.var1 != self.var2

    def is_compound(self):
        return super().is_compound() or len(self.compound) == 2

    def setup_plotter_parameters(self):
        self.build_matrices()
        self.homogenize_matrices()

    def build_matrices(self):
        cases = [
            Matrix(
                self.data.select(treatment=treatment),
                treatment,
                self.compound_or_region,
                self.var1,
                self.var2,
                self.to_plot,
                self.order,
                self.n_minimum,
                self.method,
                self.p_value_threshold,
            )
            for treatment in self.treatments
        ]  # Setup multiprocessing pool
        self.matrices = parallel_process(cases, description="Creating matrices")

    def homogenize_matrices(self):
        conserved_rows = set.intersection(
            *(set(matrix.corr_masked.index) for matrix in self.matrices)
        )
        conserved_cols = set.intersection(
            *(set(matrix.corr_masked.columns) for matrix in self.matrices)
        )

        for matrix in self.matrices:
            rows_to_drop = [
                row for row in matrix.corr_masked.index if row not in conserved_rows
            ]
            cols_to_drop = [
                col for col in matrix.corr_masked.columns if col not in conserved_cols
            ]
            matrix.corr_masked = matrix.corr_masked.drop(
                index=rows_to_drop, columns=cols_to_drop
            )

    def generate(self):
        self.setup_plotter_parameters()
        self.fig, self.axs = self.generate_figure()
        for i in range(len(self.axs)):
            self.plot_ax(i)

    def generate_figure(self):
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
        return fig, [axs] if num_rows == num_cols == 1 else axs.flatten()

    def plot_ax(self, i):
        raise NotImplementedError("Must be implemented in subclass")


@dataclass
class Correlogram(MatricesFigure):

    figure_type: str = "correlogram"

    def plot_ax(self, i):
        ax = self.axs[i]
        matrix = self.matrices[i]
        title = f"{'->'.join([self.var1, self.var2]) if self.is_square else self.var1} in {matrix.grouping}"

        ax.set_title(
            title, fontsize=28, pad=20, y=1
        )  # Adjust the y position of the title manually for square correlogram

        sns.heatmap(
            matrix.corr_masked,
            vmin=-1,
            vmax=1,
            square=True,
            # annot=True, #R value annotations
            cmap="coolwarm",
            annot_kws={"size": 8},
            ax=ax,
            cbar_kws={"shrink": 0.7},  # adj color bar size
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            ha="center",
            fontsize=12,
        )
        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=45,
            va="center",
            fontsize=12,
        )

        ax.set_ylabel(matrix.var1, fontsize=28)
        ax.set_xlabel(matrix.var2, fontsize=28)


@dataclass
class Network(MatricesFigure):

    figure_type: str = "network"
    node_positions: dict = None  # Optional parameter for manual node positions

    def define_filename(self):
        super().define_filename()
        self.filename = self.filename.replace("-", "->")

    def setup_plotter_parameters(self):
        super().setup_plotter_parameters()
        self.networks = parallel_process(
            [NetworkModel(matrix) for matrix in self.matrices],
            description="Creating networks",
        )

    def generate(self):
        super().generate()
        fig, axs = self.generate_figure()
        for i in range(len(self.axs)):
            ax = axs[i]
            network = self.networks[i]
            self.plot_degrees(ax, network)

    def plot_ax(self, i):
        ax = self.axs[i]
        network = self.networks[i]
        title = f"{'->'.join([self.var1, self.var2]) if self.is_square else self.var1} in {network.matrix.grouping}"

        positions = self.node_positions if self.node_positions else network.pos

        nx.draw_networkx_nodes(
            network.G,
            positions,
            node_size=2000,
            alpha=0.95,
            node_color="white",
            edgecolors="black",
            ax=ax,
        )
        nx.draw_networkx_edges(
            network.G,
            positions,
            width=list(nx.get_edge_attributes(network.G, "weight").values()),
            edge_color=list(nx.get_edge_attributes(network.G, "color").values()),
            ax=ax,
            node_size=2000,
            **({"arrowstyle": "->", "arrowsize": 20} if network.is_directed else {}),
        )
        # Add labels to nodes
        node_labels = {
            node: node for node in network.G.nodes()
        }  # Label nodes with their names
        nx.draw_networkx_labels(
            network.G, positions, labels=node_labels, font_size=22, ax=ax
        )
        # nx.draw_networkx_edge_labels(
        #     network.G, positions, edge_labels=network.edge_labels, font_size=18, ax=ax
        # )

        # Set the aspect ratio to 'equal' for the plot area
        ax.set_aspect("equal")
        ax.margins(0.1)

        # Set title for the graph
        ax.set_frame_on(False)
        ax.set_title(title, fontsize=28, pad=-10, y=1)

    def plot_degrees(self, ax, network):
        """
        Plots histogram of node degrees from network with a standard distribution overlay
        input:
            network object
            ax to plot
        returns:
            ax to plot
        """
        title = f"{'->'.join([self.var1, self.var2]) if self.is_square else self.var1} in {network.matrix.grouping}"
        G = network.G  # Access the graph from the Network object
        all_nodes = list(G.nodes())
        degree_sequence = [d for n, d in G.degree()]
        node_labels_with_degrees = [(n, d) for n, d in G.degree()]

        mean_degree = np.mean(degree_sequence)
        std_degree = np.std(degree_sequence)

        # Use the max_node_degree property from the Network class
        max_degree = network.max_degree
        mean_degree = network.average_degree

        x = np.linspace(0, max_degree, 100)
        y = norm.pdf(x, mean_degree, std_degree)
        ax.plot(x, y, "r-", lw=2, label=f"Standard Distribution std={std_degree:.2f}")

        # Create the histogram
        counts, bins, patches = ax.hist(
            degree_sequence,
            bins=np.arange(max_degree + 2) - 0.5,
            edgecolor="black",
            alpha=0.8,
        )
        # Check if the sum of counts matches the number of nodes
        total_nodes = len(all_nodes)
        total_counted = sum(counts)
        if total_nodes != total_counted:
            raise ValueError(
                f"Total nodes ({total_nodes}) does not match total counted ({total_counted})"
            )

        # Annotate each bar with the corresponding node labels
        for i, patch in enumerate(patches):
            bin_center = patch.get_x() + patch.get_width() / 2
            labels = [n for n, d in node_labels_with_degrees if d == i]
            if labels:
                ax.text(
                    bin_center,
                    0.06 * patch.get_height(),
                    ", ".join(labels),
                    ha="center",
                    va="bottom",
                    fontsize=28,
                    rotation=90,
                )

        ax.set_title(
            title, fontsize=28, pad=20, y=1
        )  # Use the get_title property from Network class
        ax.set_xlabel("Degree", fontsize=22)
        ax.set_ylabel("Frequency (n nodes)", fontsize=22)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=20)
        ax.tick_params(axis="x", labelsize=20)  # Adjust x-axis tick label size
        ax.tick_params(axis="y", labelsize=20)

        # HACKY PRINT
        title = network.get_title()
        print(title)
        print(
            f"edges = {network.total_edges}, pos = {network.pos_edges}, neg = {network.neg_edges}",
            f"density = {network.density}, max degree = {network.max_degree}, average degree = {network.average_degree}",
            f"unweighted clustering co = {network.avg_clust_coeff_unweighted}",
        )

        return ax


@dataclass
class Correlation(MatricesFigure):

    figure_type: str = "correlation"

    def define_filename(self):
        self.filename = f"{self.compound} in {self.region}"

    def generate(self):
        compound = self.compound if isinstance(self.compound, list) else [self.compound]
        region = self.region if isinstance(self.region, list) else [self.region]

        x_data = self.data.select(compound=compound[0], region=region[0])
        y_data = self.data.select(compound=compound[-1], region=region[-1])
        common_mouse_ids = list(set(x_data.mouse_id).intersection(set(y_data.mouse_id)))
        x_data = x_data.set_index("mouse_id").loc[common_mouse_ids].value.values
        y_data = y_data.set_index("mouse_id").loc[common_mouse_ids].value.values

        x_label = f"{compound[0]} in {region[0]} (ng/mg)"
        y_label = f"{compound[-1]} in {region[-1]} (ng/mg)"
        pearson_r, p_value = stats.pearsonr(x_data, y_data)
        color = "red" if pearson_r > 0 else "blue"

        # Create the plot
        self.fig, ax = plt.subplots()
        sns.scatterplot(x=x_data, y=y_data, ax=ax, marker="o", s=30, color="black")
        sns.regplot(
            x=x_data, y=y_data, ci=95, ax=ax, scatter=False, line_kws={"color": color}
        )

        ax.set_xlabel(x_label, fontsize=22)
        ax.set_ylabel(y_label, fontsize=22)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_title(self.treatment)

        # Add correlation values as text
        p_value_annotation = f"{p_value:.1e}" if p_value < 0.0001 else f"{p_value:.4f}"
        labels = f"Pearson R: {pearson_r:.2f}\np-value: {p_value_annotation}"
        ax.text(
            0.05,
            0.9,
            labels,
            transform=ax.transAxes,
            bbox=dict(facecolor="white", edgecolor="white", boxstyle="round"),
        )

        plt.tight_layout()
        plt.show()


@dataclass
class Table(ExcelDataset, Figure(QuantitativeDataSelection)):
    figure_type: str = field(default="table", init=False)

    def setup(self):
        self.order = REGIONS.order(self.data["region"].unique())

    def generate(self):
        grouped = (
            self.data.groupby(["region", "compound", "treatment"])
            .agg(
                mean_value=("value", "mean"),
                std_value=("value", lambda x: np.std(x, ddof=1)),
            )
            .reset_index()
        )

        # Combine mean and STD into a single string
        grouped["mean ± STD"] = grouped.apply(
            lambda row: f"{row['mean_value']:.3f} ± {row['std_value']:.3f}", axis=1
        )

        # Pivot the DataFrame
        pivot_df = grouped.pivot_table(
            index="region",
            columns=["compound", "treatment"],
            values="mean ± STD",
            aggfunc="first",
        )

        # Sort the multiindex columns
        return pivot_df.sort_index(axis=1).loc[self.order, self.compound]

    def load(self):
        return SelectableDataFrame(
            pd.read_excel(self.filepath, index_col=0, header=[0, 1])
        )


@dataclass
class StatisticsTable(Table):

    def define_filename(self):
        super().define_filename()
        self.filename += " STATS"

    def generate(self):

        
        if not self.statistics:
            return pd.DataFrame()
        
        results = []

        for statistic in self.statistics:
            data = statistic.results
            data = data[data["test"] == statistic.statistical_test][
                ["test", "region", "compound", "result_string"]
            ]
            results.append(data)
        results = pd.concat(results)
        
        results = results.pivot_table(
                index="region",
                columns=["test", "compound"],
                values="result_string",
                aggfunc="first",
            )
        
        results.index = pd.Categorical(results.index, categories=self.order, ordered=True)
        return results.sort_index()
        
        

    def load(self):
        return SelectableDataFrame(
            pd.read_excel(self.filepath, index_col=0, header=[0, 1])
        )
