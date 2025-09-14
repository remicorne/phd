from itertools import chain
from dataclasses import dataclass, field
import networkx as nx
import scipy
import pandas as pd
import numpy as np
from module.core.utils import parallel_process
from module.core.Dataset import SelectableDataFrame
from statsmodels.stats.multitest import fdrcorrection
from module.core.enums import DatasetColumn


def calculate_correlation(method, x, y):
    """Calculate the correlation and p-value based on the specified method."""
    if method == "pearson":
        result = scipy.stats.pearsonr(x, y)
        return result.statistic, result.pvalue
    elif method == "spearman":
        result = scipy.stats.spearmanr(x, y)
        return result.correlation, result.pvalue
    elif method == "kendall":
        result = scipy.stats.kendalltau(x, y)
        return result.correlation, result.pvalue
    else:
        raise ValueError(f"Unknown method: {method}")


def get_correlation_callback(method, return_type):
    """Return a correlation function based on the specified method, p-value threshold, and return type."""

    def executor(x, y):
        correlation, pvalue = calculate_correlation(method, x, y)
        if return_type == "pvalues":
            return pvalue
        elif return_type == "correlations":
            return correlation
        else:
            raise ValueError(f"Unknown return type: {return_type}")

    return executor


@dataclass
class Matrix:
    """
    Creates a reusable matrix class for a given eperimnet.
    Args:
        data (pd.DataFrame):    The original dataset.
        group (str):        Identifier for the group (logging purposes).
        variables (str):         'var1-var2' to correlate from type 'between'. If only one: self correlation.
        pivot_columns (list[str]): Columns (orederd) used to pivot the data).
        n_minimum (int): Minumum occurnces of overlapping var1 and var2 to be correlated. Default = 5.
        method (str): Correlation method ('pearson', 'spearman', 'kendall'). Default = "pearson".
        pvalue_threshold (float): Threshold for significance in correlation. Defult = 0.05

    Returns:
        filtered_data (pd.DataFrame): Subselected data filtered based on n_minimum between vairables.
        pivot (pd.DataFrame): Pivot table of the filtered data.
        corr_masked (pd.DataFrame): Masked correlation matrix based on p-value threshold.
        correlations (pd.DataFrame): Full correlation matrix.
        pvalues (pd.DataFrame): Matrix of p-values for the correlations.
        missing_values (list): List of variables with missing data (< n_minimum).
        missing_overlap (list): List of variable pairs with insufficient data overlap.

    """

    data: pd.DataFrame
    grouping: str
    pivot_columns: list[str]
    # order: list[str] = None # TODO: use pdcategorical
    between: dict
    n_minimum: int = field(kw_only=True, default=5)
    method: str = field(kw_only=True, default="pearson")
    pvalue_threshold: float = field(kw_only=True, default=0.05)
    fdr_correction: bool = field(kw_only=True, default=False)
    density_thresholding: bool = field(kw_only=True, default=False)
    delay_execution: bool = field(default=True, kw_only=True)

    filtered_data: pd.DataFrame = field(init=False)
    pivot: pd.DataFrame = field(init=False)
    corr_masked: pd.DataFrame = field(init=False)
    correlations: pd.DataFrame = field(init=False)
    pvalues: pd.DataFrame = field(init=False)
    missing_values: list = field(init=False)
    missing_overlap: list = field(init=False)

    def __call__(self):
        self.__post_init__()
        return self

    def __post_init__(self):
        if self.delay_execution:
            self.delay_execution = False
        else:
            if self.between:
                between, (self.var1, self.var2) = next(iter(self.between.items()))
                self.pivot_columns.remove(between)
                self.pivot_columns.insert(0, between)
            else:
                self.var1 = self.var2 = self.data[self.pivot_columns[0]].unique()[0]
            self.is_square = self.var1 != self.var2
            self.filter_missing_values()
            self.pivot_data()
            self.order_columns()
            self.correlate()
            self.find_missing_overlap()
            self.process_triangle_correlogram()
            self.title = self.get_title()

    def filter_missing_values(self):
        """
        Filters out variables with occurrences less than n_minimum and updates missing_values list.
        """
        self.missing_values = []
        data = []
        for measurement_characteristics, df in self.data.groupby(by=self.pivot_columns):
            if df.value.notna().sum() < self.n_minimum:
                self.missing_values.append(measurement_characteristics)
            else:
                data.append(df)
        if self.missing_values:
            print(
                f"{self.grouping}, {self.between} missing data for {self.missing_values}, deleted from analysis"
            )
        self.filtered_data = pd.concat(data)

    def pivot_data(self):
        """
        Creates a pivot table from the filtered data.
        """
        self.pivot = self.filtered_data.pivot_table(
            values=DatasetColumn.VALUE,
            index="subject_id",
            columns=self.pivot_columns,
        )

    def order_columns(self):
        """
        Orders the columns of the pivot table based on the provided column list.
        """
        # columns = (
        #     sorted(
        #         self.pivot.columns,
        #         key=lambda x: (
        #             self.order.index(x[1]) if x[1] in self.order else float("inf")
        #         ),
        #     )
        #     if self.order
        #     else self.pivot.columns
        # )
        # self.pivot = self.pivot[columns]

    def correlate(self):
        """
        Creates self.correlations and self.pvalues self.corr_masked (specific to normalisation).
        """
        self.pvalues = self.create_corr_matrix("pvalues")
        self.correlations = self.create_corr_matrix("correlations")

        if (
            self.density_thresholding is not None
        ):  # Density-based masking (ignore p-values)
            if self.between:
                corr_flat = self.correlations.abs().stack()
            else:
                corr_values = self.correlations.where(
                    ~np.eye(self.correlations.shape[0], dtype=bool)
                )
                corr_flat = corr_values.abs().stack()

            cutoff = corr_flat.quantile(1 - self.density_thresholding)
            mask = self.correlations.abs() >= cutoff

        else:  #  P-value based masking
            if self.fdr_correction:
                self.uncorrected_pvalues = self.pvalues
                self.pvalues = self.apply_fdr_correction()

            mask = self.pvalues < self.pvalue_threshold
        self.corr_masked = self.correlations.where(mask, other=np.nan)

    # def correlate(self): # REMI CHECK I have replaed the bellow with the function above
    #     """
    #     Calculates and stores correlation and p-value matrices.
    #     """
    #     self.pvalues = self.create_corr_matrix("pvalues")
    #     self.correlations = self.create_corr_matrix("correlations")
    #     # self.corr_masked = self.correlations[self.pvalues < self.pvalue_threshold]

    #     if self.fdr_correction:
    #         self.uncorrected_pvalues = self.pvalues
    #         self.pvalues = self.apply_fdr_correction()

    #     mask = self.pvalues < self.pvalue_threshold
    #     self.corr_masked = self.correlations.where(mask, other=np.nan)

    def apply_fdr_correction(self):
        """
        Applies Benjamini-Hochberg FDR correction to p-values.
        """
        pvalues_corrected_matrix = np.full(
            self.pvalues.shape, np.nan
        )  # Start with NaN matrix

        if self.is_square:  # apply FDR to the entire matrix
            p_flat = self.pvalues.values.flatten()
            _, p_corrected = fdrcorrection(p_flat, method="indep")
            pvalues_corrected_matrix = p_corrected.reshape(
                self.pvalues.shape
            )  # Reshape back
        else:
            triu_indices = np.triu_indices_from(self.pvalues, k=1)
            p_flat = self.pvalues.values[triu_indices]
            _, p_corrected = fdrcorrection(p_flat, method="indep")

            pvalues_corrected_matrix[triu_indices] = p_corrected
            pvalues_corrected_matrix[triu_indices[::-1]] = (
                p_corrected  # Copy symmetrically
            )
            np.fill_diagonal(pvalues_corrected_matrix, self.pvalues.values.diagonal())

        print(f"Matrix for: {self.grouping}")
        print(
            f"Significant correlations before FDR correction: {np.sum(self.pvalues.values < self.pvalue_threshold)}"
        )
        print(
            f"Significant correlations after FDR correction: {np.sum(pvalues_corrected_matrix < self.pvalue_threshold)}"
        )
        return pd.DataFrame(
            pvalues_corrected_matrix,
            index=self.pvalues.index,
            columns=self.pvalues.columns,
        )

    def create_corr_matrix(self, result_type):
        """
        Creates a correlation matrix for either correlation values or p-values.

        Args:
            result_type (str): Type of result to return ('pvalues' or 'correlations').

        Returns:
            pd.DataFrame: A DataFrame containing the requested correlation matrix.
        """
        method = get_correlation_callback(self.method, result_type)
        matrix = self.pivot.corr(method=method, min_periods=self.n_minimum)
        return matrix.loc[
            self.var1, self.var2
        ]  # IMPROVE: use .corrwith to only caluclate necessary correlation for square corr

    def find_missing_overlap(self):
        """
        Identifies and reports variable pairs with insufficient data overlap.
        """
        stack = self.correlations.stack(dropna=False)
        stack.index = stack.index.rename(self.pivot_columns)
        stack = pd.DataFrame(stack.reset_index())
        stack.columns = list(stack.columns[:-1]) + [DatasetColumn.VALUE]
        stack = stack[stack.value.isna()]
        self.missing_overlap = stack[stack.value.isna()][self.pivot_columns].values
        if len(self.missing_overlap):
            print(
                f"{self.grouping} {self.between} insuficient overlapp for {self.missing_overlap} pairs"
            )
            print("Inspect with self.corr to adjust {columns} and redo analysis")

    def process_triangle_correlogram(self):
        """
        Masks the upper triangle of the correlogram if the matrix is not square.
        """
        if not self.is_square:
            mask = np.triu(np.ones(self.corr_masked.shape, dtype=bool), k=1)
            self.corr_masked[mask] = np.nan
            np.fill_diagonal(self.corr_masked.values, 1)

    @property
    def significant_correlations(self):
        return self.corr_masked.stack().items()

    def get_title(self):
        return f"{self.var1 if self.var1 == self.var2 else '->'.join([self.var1, self.var2])} in {self.grouping}"


@dataclass
class Network:
    """
    A class to represent a network/graph constructed from a correlation matrix.

    Attributes:
    matrix (Matrix): An instance of the Matrix class containing the data and correlation matrix.

    Methods:
    max_node_degree(): Returns the maximum degree of the nodes in the graph.
    is_directed(): Property that checks if the network is directed (based on the matrix being square).
    plot_ax(ax): Plots the network graph on the given matplotlib axis.
    """

    matrix: Matrix
    delay_execution: bool = field(default=True, kw_only=True)

    def __call__(self):
        self.__post_init__()
        return self

    def __post_init__(self):
        """
        Initializes the network/graph from the given matrix. Constructs a directed or undirected graph
        based on the matrix's properties and fills it with nodes and edges based on the correlation data.
        """
        if self.delay_execution:
            self.delay_execution = False
        else:
            self.title = self.matrix.title
            self.grouping = self.matrix.grouping
            self.between = self.matrix.between
            self.is_directed = self.matrix.is_square
            self.G = nx.MultiDiGraph() if self.matrix.is_square else nx.Graph()
            # directed edge -  to_correlate[0] --> to_correlate[1]
            self.G.clear()

            self.G.add_nodes_from(
                set(
                    self.matrix.corr_masked.columns.tolist()
                    + self.matrix.corr_masked.index.tolist()
                )
            )
            self.edge_labels = {}
            for (row, col), correlation in self.matrix.significant_correlations:
                if not (row == col and not self.is_directed):
                    self.G.add_edge(
                        row,
                        col,
                        weight=correlation,
                        color="red" if correlation > 0 else "blue",
                    )
                    self.edge_labels[(row, col)] = f"{correlation:.2f}"

            self.density = self.calculate_graph_density()
            (
                self.total_edges,
                self.pos_edges,
                self.neg_edges,
                self.neg_pos_edge_ratio,
                self.neg_edge_density,
            ) = self.calculate_edge_count()

            self.max_degree, self.average_degree, self.min_degree = (
                self.calculate_node_degree()
            )
            self.SD_node_degree, self.SD_node_strength = (
                self.calculate_SD_node_degree_strength()
            )

            self.clust_coeff_unweighted, self.clust_coeff_weighted = (
                self.calculate_clustering_coefficient()
            )
            self.global_efficiency_weighted, self.global_efficiency_unweighted = (
                self.calculate_global_efficiency()
            )  # unweighted
            self.local_efficiency_weighted, self.local_efficiency_unweighted = (
                self.calculate_local_efficiency()
            )  # unweighted

    def calculate_SD_node_degree_strength(self):
        """
        standard deviation of the node degrees (float)
        """
        degrees = dict(self.G.degree())
        strengths = dict(self.G.degree(weight="weight"))

        return np.std(list(degrees.values())), np.std(list(strengths.values()))

    def calculate_edge_count(self):
        total_edges = self.G.number_of_edges()
        total_nodes = self.G.number_of_nodes()

        pos_edges = 0
        neg_edges = 0
        for u, v, data in self.G.edges(data=True):
            color = data.get("color")
            if color == "red":
                pos_edges += 1
            elif color == "blue":
                neg_edges += 1

        if pos_edges > 0 and neg_edges > 0:
            neg_pos_edge_ratio = neg_edges / total_edges
        else:
            neg_pos_edge_ratio = 0

        neg_edge_density = neg_edges / (total_nodes * (total_nodes - 1))

        return total_edges, pos_edges, neg_edges, neg_pos_edge_ratio, neg_edge_density

    def calculate_node_degree(self):
        degrees = dict(self.G.degree())
        max_degree = max(degrees.values())
        min_degree = min(degrees.values())
        average_degree = np.mean(list(degrees.values()))
        return max_degree, average_degree, min_degree

    def calculate_graph_density(self):
        """
        graph_density (float) ie edges/all_possible_edges
        """
        num_edges = self.G.number_of_edges()
        num_nodes = self.G.number_of_nodes()
        if self.is_directed:  # directed graph have doubble possible edges
            max_edges = num_nodes * (num_nodes - 1)
        else:
            max_edges = num_nodes * (num_nodes - 1) / 2
        return num_edges / max_edges if max_edges > 0 else 0

    def calculate_clustering_coefficient(self):
        """
        Calculates the average clustering coefficient for both directed and undirected graphs.
        Handles both weighted and unweighted cases.

        Returns:
            avg_clust_coeff_unweighted (float): Average unweighted clustering coefficient for the graph.
            avg_clust_coeff_weighted (float): Average weighted clustering coefficient for the graph.
        """
        # For undirected graphs (Graph), we can directly use NetworkX's clustering function
        if not self.is_directed:
            clust_coeff_unweighted = nx.clustering(
                self.G
            )  # Unweighted clustering coefficient
            clust_coeff_weighted = nx.clustering(
                self.G, weight="weight"
            )  # Weighted clustering coefficient
        else:
            # For directed graphs (MultiDiGraph), clustering calculation requires special handling
            if isinstance(self.G, nx.MultiDiGraph):
                # Convert MultiDiGraph to DiGraph, as NetworkX doesn't support clustering on MultiDiGraph directly
                # Here we take the first edge between each pair of nodes (if multiple edges exist)
                simple_directed_G = (
                    nx.DiGraph()
                )  # Create a simple DiGraph from the MultiDiGraph
                for u, v, data in self.G.edges(data=True):
                    if not simple_directed_G.has_edge(
                        u, v
                    ):  # Only add the first edge between nodes
                        simple_directed_G.add_edge(u, v, weight=data["weight"])

                clust_coeff_unweighted = nx.clustering(
                    simple_directed_G
                )  # Unweighted clustering coefficient
                clust_coeff_weighted = nx.clustering(
                    simple_directed_G, weight="weight"
                )  # Weighted clustering coefficient
            else:
                # For standard directed graphs (DiGraph), just use the regular directed clustering method
                clust_coeff_unweighted = nx.clustering(
                    self.G
                )  # Unweighted clustering coefficient
                clust_coeff_weighted = nx.clustering(
                    self.G, weight="weight"
                )  # Weighted clustering coefficient

        # Calculate average clustering coefficient
        avg_clust_coeff_unweighted = (
            sum(clust_coeff_unweighted.values()) / len(clust_coeff_unweighted)
            if clust_coeff_unweighted
            else 0
        )
        avg_clust_coeff_weighted = (
            sum(clust_coeff_weighted.values()) / len(clust_coeff_weighted)
            if clust_coeff_weighted
            else 0
        )

        return avg_clust_coeff_unweighted, avg_clust_coeff_weighted

    def calculate_global_efficiency(self):
        """
        Calculates the global efficiency of the graph (unweighted only).

        Returns:
            global_efficiency_weighted (float): Global efficiency using edge weights.
            global_efficiency_unweighted (float): Global efficiency for the unweighted graph.
        """
        n = self.G.number_of_nodes()
        if n <= 1:
            return 0, 0

        # Weighted: convert weights to distances (larger weight = shorter distance).
        G_weighted = self.G.copy()
        for u, v, d in G_weighted.edges(data=True):
            d["distance"] = 1 / max(d.get("weight", 1e-6), 1e-6)

        # if directed shortest path respects directionality
        weighted_lengths = dict(
            nx.all_pairs_dijkstra_path_length(G_weighted, weight="distance")
        )
        unweighted_lengths = dict(nx.all_pairs_shortest_path_length(self.G))

        # Unreachable pairs are skipped
        weighted_vals = [
            1 / l
            for src in weighted_lengths
            for tgt, l in weighted_lengths[src].items()
            if src != tgt and l > 0
        ]
        unweighted_vals = [
            1 / l
            for src in unweighted_lengths
            for tgt, l in unweighted_lengths[src].items()
            if src != tgt and l > 0
        ]

        return np.mean(weighted_vals) if weighted_vals else 0, np.mean(
            unweighted_vals
        ) if unweighted_vals else 0

    def calculate_local_efficiency(self):
        """
        Calculates the local efficiency of the graph.

        Local efficiency is the average efficiency of each node’s neighbors.
        Returns:
            local_efficiency_weighted (float): Local efficiency using edge weights.
            local_efficiency_unweighted (float): Local efficiency ignoring weights.
        """
        local_efficiency_unweighted = []
        local_efficiency_weighted = []

        for node in self.G.nodes():
            neighbors = list(self.G.neighbors(node))
            if len(neighbors) < 2:
                continue  # Need at least two neighbors to compute efficiency

            # Subgraph of neighbors
            subgraph = self.G.subgraph(neighbors)

            # --- Unweighted ---
            unweighted_lengths = dict(nx.all_pairs_shortest_path_length(subgraph))
            unweighted_vals = [
                1 / l
                for src in unweighted_lengths
                for tgt, l in unweighted_lengths[src].items()
                if src != tgt and l > 0
            ]
            if unweighted_vals:
                local_efficiency_unweighted.append(np.mean(unweighted_vals))

            # --- Weighted ---
            subgraph_w = subgraph.copy()
            for u, v, d in subgraph_w.edges(data=True):
                d["distance"] = 1 / max(d.get("weight", 1e-6), 1e-6)

            weighted_lengths = dict(
                nx.all_pairs_dijkstra_path_length(subgraph_w, weight="distance")
            )
            weighted_vals = [
                1 / l
                for src in weighted_lengths
                for tgt, l in weighted_lengths[src].items()
                if src != tgt and l > 0
            ]
            if weighted_vals:
                local_efficiency_weighted.append(np.mean(weighted_vals))

        avg_local_unweighted = (
            np.mean(local_efficiency_unweighted) if local_efficiency_unweighted else 0
        )
        avg_local_weighted = (
            np.mean(local_efficiency_weighted) if local_efficiency_weighted else 0
        )

        return avg_local_weighted, avg_local_unweighted


@dataclass
class MatrixGroup:
    """
    Creates a collection of Matrix objects from a larger dataset. This class
    helps in processing and analyzing data by grouping, selecting relevant variables,
    and building individual matrices for further analysis.

    Attributes:
        data (pd.DataFrame): The original dataframe containing the data.
        group_by (str): The column name in 'data' to group by (generally 'experiment').
        between (str): The first variable to correlate (e.g., compound or region).
        variables (str): The specific variables to correlate from 'between'.
        accross (str): The second variable to correlate against 'between'.
        sub_selector (str): Additional filtering criteria for sub-selecting the data.
        columns (list[str]): Columns to select from the 'accross' column. Defaults to None.
        n_minimum (int): Minimum number of occurrences for a valid correlation. Defaults to 5.
        method (str): Correlation method, one of 'pearson', 'kendall', 'spearman'. Defaults to "pearson".
        pvalue_threshold (float): P-value threshold for significance. Defaults to 0.05.

    Returns:
        matrices (list[Matrix]): A list of Matrix objects created from the grouped data.
        var1 (str): The first variable derived from 'variables'.
        var2 (str): The second variable derived from 'variables'.
    """

    data: pd.DataFrame
    group_by: str
    pivot_columns: list[str]
    between: dict
    # order: list[str] = None
    n_minimum: int = field(kw_only=True, default=5)
    method: str = field(kw_only=True, default="pearson")
    pvalue_threshold: float = field(kw_only=True, default=0.05)
    fdr_correction: float = field(kw_only=True, default=None)
    density_thresholding: bool = field(kw_only=True, default=False)

    matrices: list[Matrix] = field(init=False)

    def __post_init__(self):
        self.build_matrices()
        self.homogenize_datasets()

    def build_matrices(self):
        batch = []
        col, cases = next(iter(self.between.items()))
        for col, cases in self.between.items():
            for group, group_df in self.data.groupby(by=self.group_by, sort=False):
                for between in cases:
                    batch.append(
                        Matrix(
                            group_df.select(**{col: between}),
                            group,
                            self.pivot_columns,
                            between={col: tuple(between)},
                            n_minimum=self.n_minimum,
                            method=self.method,
                            pvalue_threshold=self.pvalue_threshold,
                            fdr_correction=self.fdr_correction,
                            density_thresholding=self.density_thresholding,
                        )
                    )
        self.matrices = parallel_process(batch, description="Building matrices")

    def homogenize_datasets(self):
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

    def __iter__(self):
        for matrix in self.matrices:
            yield matrix


@dataclass
class NetworkGroup:
    matrix_group: MatrixGroup

    def __post_init__(self):
        self.networks = parallel_process(
            [Network(matrix) for matrix in self.matrix_group],
            description="Creating networks",
        )
        self.nodes = set(chain(*[network.G.nodes for network in self.networks]))

    def get_summary_df(self):
        data = []
        for network in self.networks:
            row = {self.matrix_group.group_by: network.grouping, **network.between}
            for variable in [
                "density",
                "neg_edge_density",
                "total_edges",
                "pos_edges",
                "neg_edges",
                "neg_pos_edge_ratio",
                "max_degree",
                "average_degree",
                "min_degree",
                "SD_node_degree",
                "SD_node_strength",
                "clust_coeff_unweighted",
                "clust_coeff_weighted",
                "global_efficiency_weighted",
                "global_efficiency_unweighted",
                "local_efficiency_weighted",
                "local_efficiency_unweighted",
            ]:
                row = {
                    **row,
                    "measurement": variable,
                    DatasetColumn.VALUE: getattr(network, variable),
                }
                data.append(row)
        return SelectableDataFrame(data)

    def __len__(self):
        return len(self.networks)

    def __getitem__(self, index):
        return self.networks[index]
