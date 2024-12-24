from dataclasses import dataclass, field
import os
import networkx as nx
from scipy.stats import norm
from scipy import stats
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from module.core.Dataset import SelectableDataFrame, ExcelCachedDataFrame
from module.core.questions import input_escape
from IPython.display import display
import ipywidgets as widgets

from module.core.Matrix import Matrix, Network
from module.core.Matrix import Network as NetworkModel
from module.core.Constants import ConstantRegistry
from matplotlib import pyplot as plt
import seaborn as sns
from typing import ClassVar
import numpy as np
from IPython.display import Image, display
import pandas as pd
import os
from module.core.utils import parallel_process
from statannotations.Annotator import Annotator
from module.core.Statistics import QuantitativeStatistic
from typing import ClassVar


@dataclass
class Figure:

    title: str
    filepath: str
    custom_params: dict = field(kw_only=True, default_factory=dict)

    fig: plt.Figure = field(init=False)

    def __post_init__(
        self,
    ):
        self.generate_figure()
        self.plot()
        self.save()

    def generate_figure(self):
        raise NotImplementedError

    def plot(self):
        raise NotImplementedError

    def save(self):
        filepath, _ = os.path.splitext(self.filepath)
        os.makedirs(filepath, exist_ok=True)
        self.fig.savefig(f"{filepath}.svg")
        self.fig.savefig(f"{filepath}.png")
        print(f"SAVED {filepath}")


@dataclass
class Histogram(Figure):
    """
    Generate a histogram of treatments. If only one compound or region is specified, a simple histogram is generated.
    If multiple compounds or regions are specified, a summary histogram is generated.
    """

    figure_type: ClassVar[str] = "histogram"

    data: pd.DataFrame
    x: str
    hue: str
    statistic: QuantitativeStatistic = field(default=None)

    def generate_figure(self):
        self.fig, self.ax = plt.subplots(figsize=(20, 10))

    def plot(self):
        if self.custom_params.get("plot_bar", True):
            sns.barplot(
                data=self.data,
                x=self.x,
                y="value",
                hue=self.hue,
                palette=self.custom_params.get("palette"),
                errorbar=self.custom_params.get("errorbar", "sd"),
                edgecolor=self.custom_params.get("edgecolor", ".2"),
                errcolor=self.custom_params.get("errcolor", ".2"),
                capsize=self.custom_params.get("capsize", 0.1),
                alpha=self.custom_params.get("alpha", 0.8),
                order=self.custom_params.get("order"),
                dodge=self.custom_params.get("dodge", False),
            )
        if self.custom_params.get("plot_swarm", True):
            sns.swarmplot(
                data=self.data,
                x=self.x,
                y="value",
                hue=self.custom_params.get("swarm_hue", self.hue),
                size=self.custom_params.get("size", 5),
                palette=self.custom_params.get("palette"),
                # legend=False if self.custom_params.get("plot_bar") else "auto",
                order=self.custom_params.get("order"),
                edgecolor=self.custom_params.get("edgecolor", "k"),
                linewidth=self.custom_params.get("linewidth", 1),
                linestyle=self.custom_params.get("linestyle", "-"),
                dodge=self.custom_params.get("dodge", False),
            )

        self.ax.tick_params(labelsize=self.custom_params.get("labelsize", 24))
        self.ax.set_ylabel(
            self.custom_params.get("ylabel"),
            fontsize=self.custom_params.get("ylabel_fontsize", 24),
        )
        self.ax.set_xlabel(
            " ", fontsize=self.custom_params.get("xlabel_fontsize", 20)
        )  # treatments
        self.ax.set_title(
            self.title,
            y=self.custom_params.get("y", 1.04),
            fontsize=self.custom_params.get("fontsize", 34),
        )
        sns.despine(left=False)
        self.label_histogram_stats()

    def label_histogram_stats(self):
        if self.statistic and self.statistic.is_significant:
            pairs, p_values = self.statistic.significant_pairs
            annotator = Annotator(
                self.ax,
                pairs,
                data=self.data,
                x=self.x,
                y="value",
                order=self.custom_params.get("hue_order"),
            )
            annotator.configure(text_format="star", loc="inside", fontsize="xx-large")
            annotator.set_pvalues_and_annotate(p_values)


@dataclass
class SummaryHistogram(Figure):
    """
    Generate a histogram of treatments. If only one compound or region is specified, a simple histogram is generated.
    If multiple compounds or regions are specified, a summary histogram is generated.
    """

    figure_type: ClassVar[str] = "summary_histogram"

    data: pd.DataFrame
    x: str
    hue: str
    statistics: list[QuantitativeStatistic] = field(default_factory=list)

    def generate_figure(self):
        self.fig_width = self.custom_params.get(
            "fig_width", 1 + 4 * len(self.custom_params.get("order"))
        )
        self.fig, self.ax = plt.subplots(figsize=(self.fig_width, 10))

    def plot(self):
        if self.custom_params.get("plot_bar", True):
            sns.barplot(
                data=self.data,
                x=self.x,
                y="value",
                hue=self.hue,
                palette=self.custom_params.get("palette"),
                errorbar=self.custom_params.get("errorbar", "sd"),
                edgecolor=self.custom_params.get("edgecolor", ".2"),
                errcolor=self.custom_params.get("errcolor", ".2"),
                capsize=self.custom_params.get("capsize", 0.1),
                alpha=self.custom_params.get("alpha", 0.8),
                order=self.custom_params.get("order"),
                hue_order=self.custom_params.get("hue_order"),
                errwidth=self.custom_params.get("errwidth", 1),
                dodge=self.custom_params.get("dodge", True),
                width=self.custom_params.get("bar_width", 0.8),
            )
        if self.custom_params.get("plot_swarm", False):
            sns.swarmplot(
                data=self.data,
                x=self.x,
                y="value",
                hue=self.custom_params.get("swarm_hue", self.hue),
                hue_order=self.custom_params.get("hue_order"),
                palette=self.custom_params.get("palette"),
                errorbar=self.custom_params.get("errorbar", "sd"),
                errcolor=self.custom_params.get("errcolor", ".2"),
                capsize=self.custom_params.get("capsize", 0.1),
                alpha=self.custom_params.get("alpha", 0.8),
                order=self.custom_params.get("order"),
                errwidth=self.custom_params.get("errwidth", 1),
                width=self.custom_params.get("bar_width", 0.8),
                legend=False if self.custom_params.get("plot_bar") else "auto",
                edgecolor=self.custom_params.get("edgecolor", "k"),
                linewidth=self.custom_params.get("linewidth", 1),
                linestyle=self.custom_params.get("linestyle", "-"),
                dodge=self.custom_params.get("dodge", False),
            )

        self.ax.tick_params(labelsize=20)
        self.ax.set_ylabel(
            self.custom_params.get("ylabel"),
            fontsize=20,
            labelpad=self.custom_params.get("labelpad", 100),
        )
        self.ax.yaxis.set_label_coords(
            self.custom_params.get("ylabel_x", -0.5 / self.fig_width), 0.5
        )
        self.ax.set_xlabel(" ", fontsize=15)  # remove x title
        self.ax.set_title(self.custom_params.get("title"), y=1.04, fontsize=34)
        self.ax.legend(loc="upper right")  # , bbox_to_anchor=(0.1, 1))
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["right"].set_visible(False)
        plt.tight_layout()
        self.label_summary_stats()

    def label_summary_stats(self):
        for statistics in self.statistics:
            if statistics.is_significant:
                # Font Scaling # HARDCODE JJB TODO - also add significance pairs!
                base_font_size = 48
                scaling_factor = 0.2
                dynamic_font_size = max(
                    base_font_size
                    - (scaling_factor * len(self.custom_params.get("order"))),
                    6,
                )
                for pair in statistics.significant_pairs[0]:
                    for i, (treatment, symbol) in enumerate(
                        self.custom_params.get("significance_palette").items()
                    ):
                        if treatment in pair:
                            hue = pair[
                                pair.index(treatment) - 1
                            ]  # work because only two elements 0 -> -1, 1 -> 0
                            x_index = self.custom_params.get("order").index(
                                statistics.metadata[
                                    self.x
                                ]  # TODO too much knowledge of internal objects
                            )  # Statistics metadata stores grouping info
                            hue_index = self.custom_params.get("hue_order").index(hue)
                            bar = self.ax.patches[
                                hue_index * len(self.custom_params.get("order"))
                                + x_index
                            ]
                            self.ax.text(
                                bar.get_x() + bar.get_width() / 2,
                                (bar.get_height() * (1.3 + i / 5)),
                                symbol,
                                ha="center",
                                va="bottom",
                                fontsize=dynamic_font_size,
                            )
                            break


class MultiAxFigure(Figure):

    @property
    def num_axs(self):
        raise NotImplementedError("Must be implemented in subclass")

    def generate_figure(self):
        # determin number of treatments to corrispond to number of subplots
        num_cols = min(int(np.sqrt(self.num_axs)), 2)  # max of 2 columns
        num_rows = (
            self.num_axs + num_cols - 1
        ) // num_cols  # Compute the number of rows

        # define the base size and a scaling factor for the figure size
        square_side = 11

        # create subplots
        self.fig, axs = plt.subplots(
            num_rows,
            num_cols,
            figsize=(
                self.custom_params.get("width", num_cols * square_side),
                self.custom_params.get("height", num_rows * square_side),
            ),
            constrained_layout=True,
        )
        self.axs = [axs] if num_rows == num_cols == 1 else axs.flatten()
        self.fig.suptitle(self.title, fontsize=28)

    def plot(self):
        for i in range(self.num_axs):
            self.plot_ax(i)

    def plot_ax(self, i):
        raise NotImplementedError("Must be implemented in subclass")


@dataclass
class Correlogram(MultiAxFigure):

    figure_type: ClassVar[str] = "correlogram"

    matrices: list[Matrix]

    @property
    def num_axs(self):
        return len(self.matrices)

    def plot_ax(self, i):

        ax = self.axs[i]
        matrix = self.matrices[i]

        colormap = self.custom_params.get("colormap", "coolwarm")

        ax.set_title(
            self.title, fontsize=28, pad=20, y=1
        )  # Adjust the y position of the title manually for square correlogram

        sns.heatmap(
            matrix.corr_masked,
            vmin=-1,
            vmax=1,
            square=True,
            # annot=True, #R value annotations
            cmap=colormap,
            annot_kws={"size": 8},
            ax=ax,
            cbar_kws={"shrink": 0.7},  # adj color bar size
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=90,
            ha="center",
            fontsize=16,
        )
        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=0,
            va="center",
            fontsize=16,
        )

        ax.set_ylabel(matrix.var1, fontsize=28)
        ax.set_xlabel(matrix.var2, fontsize=28)


@dataclass
class NetworkFigure(MultiAxFigure):

    figure_type: ClassVar[str] = "network"

    networks: list[Network]
    positions: dict = field(kw_only=True, default=None)

    @property
    def num_axs(self):
        return len(self.networks)

    def plot_ax(self, i):
        show_edge_labels = self.custom_params.get("show_edge_labels", False)
        edge_thickness = self.custom_params.get(
            "edge_thickness", 3
        )  # 'weight' for thickness weighting
        colormap = self.custom_params.get("colormap", "coolwarm")

        ax = self.axs[i]
        network = self.networks[i]
        if not self.positions:
            self.positions = self.get_default_positions(network.matrix)
            ax.set_xlim(0, 27)
            ax.set_ylim(0, 15)

        nx.draw_networkx_nodes(
            network.G,
            self.positions,
            node_size=2000,
            alpha=0.95,
            node_color="white",
            edgecolors="black",
            ax=ax,
        )

        edge_weights = list(nx.get_edge_attributes(network.G, "weight").values())

        if edge_thickness == "weight":  # display weight by line thickness
            edge_colors = list(nx.get_edge_attributes(network.G, "color").values())
            weight_scaler = 3  # this should be log #TODO
            edge_weight_to_plot = [weight * weight_scaler for weight in edge_weights]

        else:  # display weight by colormap
            normalize = Normalize(vmin=-1, vmax=1)
            cmap = plt.get_cmap(colormap)
            edge_colors = [cmap(normalize(weight)) for weight in edge_weights]
            edge_weight_to_plot = edge_thickness
            sm = ScalarMappable(cmap=cmap, norm=normalize)
            plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.04)  # label='Edge Weight',

        nx.draw_networkx_edges(
            network.G,
            self.positions,
            width=edge_weight_to_plot,
            edge_color=edge_colors,
            ax=ax,
            node_size=2000,
            **({"arrowstyle": "->", "arrowsize": 20} if network.is_directed else {}),
        )

        # Label nodes
        node_labels = {node: node for node in network.G.nodes()}
        nx.draw_networkx_labels(
            network.G, self.positions, labels=node_labels, font_size=24, ax=ax
        )

        if show_edge_labels is True:
            rounded_edge_labels = {
                edge: f"{weight:.1f}"
                for edge, weight in nx.get_edge_attributes(network.G, "weight").items()
            }

            if edge_thickness == "weight":
                color_labels = {"red": [], "blue": []}
                for edge, weight in rounded_edge_labels.items():
                    color = nx.get_edge_attributes(network.G, "color")[edge]
                    color_labels[color].append((edge, weight))

                for color, edges in color_labels.items():
                    edge_labels = {edge: label for edge, label in edges}
                    nx.draw_networkx_edge_labels(
                        network.G,
                        self.positions,
                        edge_labels=edge_labels,
                        font_size=14,
                        ax=ax,
                        font_color=color,
                        label_pos=0.5,
                        bbox=dict(
                            facecolor="white",
                            edgecolor="none",
                            boxstyle="round,pad=0.01",
                            alpha=0.7,
                        ),  # White background box
                    )

            else:
                edge_colors = [
                    cmap(weight)
                    for weight in nx.get_edge_attributes(network.G, "weight").values()
                ]
                edge_label_colors = {}
                for (edge, weight), color in zip(
                    rounded_edge_labels.items(), edge_colors
                ):
                    edge_label_colors[edge] = color

                for edge, color in edge_label_colors.items():
                    label = rounded_edge_labels[edge]
                    nx.draw_networkx_edge_labels(
                        network.G,
                        self.positions,
                        edge_labels={
                            edge: label
                        },  # Pass a dict with the edge and its label
                        font_size=14,
                        ax=ax,
                        font_color=color,  # Use the unique color for each label
                        label_pos=0.5,
                        bbox=dict(
                            facecolor="white",
                            edgecolor="none",
                            boxstyle="round,pad=0.01",
                            alpha=0.9,
                        ),
                    )

        ax.set_aspect("equal")
        ax.margins(0.1)
        ax.set_frame_on(False)
        ax.set_title(network.title, fontsize=28, pad=-10, y=1)

    def get_default_positions(self, matrix):
        angles = np.linspace(
            0, 2 * np.pi, len(matrix.corr_masked.columns), endpoint=False
        )
        return {
            col: (np.cos(angles[i]), np.sin(angles[i]))
            for i, col in enumerate(matrix.corr_masked.columns)
        }


@dataclass
class NetworkDegreesFigure(MultiAxFigure):

    figure_type: ClassVar[str] = "network_degrees"

    networks: list[Network]

    @property
    def num_axs(self):
        return len(self.networks)

    def plot_ax(self, i):
        """
        Plots histogram of node degrees from network with a standard distribution overlay
        input:
            network object
            ax to plot
        returns:
            ax to plot
        """
        network = self.networks[i]
        ax = self.axs[i]
        # Access the graph from the Network object
        all_nodes = list(network.G.nodes())
        degree_sequence, node_labels_with_degrees = [], []
        for n, d in network.G.degree():
            degree_sequence.append(d)
            node_labels_with_degrees.append((n, d))

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

        ax.set_title(network.title, fontsize=28, pad=20, y=1)
        ax.set_xlabel("Degree", fontsize=22)
        ax.set_ylabel("Frequency (n nodes)", fontsize=22)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=20)
        ax.tick_params(axis="x", labelsize=20)  # Adjust x-axis tick label size
        ax.tick_params(axis="y", labelsize=20)

        # HACKY PRINT
        print(network.grouping)
        print(
            f"edges = {network.total_edges}, pos = {network.pos_edges}, neg = {network.neg_edges}",
            f"density = {network.density}, max degree = {network.max_degree}, average degree = {network.average_degree}",
            f"unweighted clustering co = {network.avg_clust_coeff_unweighted}",
        )


@dataclass
class Correlation(Figure):

    figure_type: ClassVar[str] = "correlation"

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
class Table(ExcelCachedDataFrame, Figure):
    figure_type: str = field(default="table", init=False)

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
        return pivot_df.sort_index(axis=1).loc[self.order]

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

        results.index = pd.Categorical(
            results.index, categories=self.order, ordered=True
        )
        return results.sort_index()

    def load(self):
        return SelectableDataFrame(
            pd.read_excel(self.filepath, index_col=0, header=[0, 1])
        )
