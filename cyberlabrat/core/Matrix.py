from dataclasses import dataclass, field
import networkx as nx
import scipy
import pandas as pd
import numpy as np
from cyberlabrat.core.utils import parallel_process




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


def correlate(method, return_type):
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
        treatment (str):        Identifier for the treatment group.
        between (str):          Variable type for correlation (e.g., 'compound' or 'region').
        variables (str):         'var1-var2' to correlate from type 'between'. If only one: self correlation.
        accross (str): The column that will constitute the rows/cols of the matrix.
        columns (list[str]): Columns to include in the analysis. If None, all columns are included.
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
    between: str
    var1: str
    var2: str
    accross: str
    order: list[str] = None
    n_minimum: int = 5
    method: str = "pearson"
    pvalue_threshold: float = 0.05
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
            self.is_square = self.var1 != self.var2
            self.filter_missing_values()
            self.pivot_data()
            self.order_columns()
            self.correlate()
            self.find_missing_overlap()
            self.process_triangle_correlogram()
            self.get_title

    def get_title(self):
        """
        Generates a title for the correlogram based on the matrix configuration.

        Returns:
            str: A title string.
        """
        if self.is_square:
            return f"{'-'.join([self.var1, self.var2])} in {self.grouping}"
        return f"{self.var1} in {self.grouping}"  
    
    def filter_missing_values(self):
        """
        Filters out variables with occurrences less than n_minimum and updates missing_values list.
        """
        self.missing_values = []
        missing_indices = []
        for col, df in self.data.groupby(by=[self.between, self.accross]):
            if df.value.notna().sum() < self.n_minimum:
                self.missing_values.append(col)
                missing_indices.extend(df.index)
        self.filtered_data = self.data.drop(missing_indices)
        if self.missing_values:
            print(
                f"{self.grouping} missing data for {self.missing_values}, deleted from analysis"
            )

    def pivot_data(self):
        """
        Creates a pivot table from the filtered data.
        """
        self.pivot = self.filtered_data.pivot_table(
            values="value",
            index=self.filtered_data["mouse_id"],
            columns=[self.between, self.accross],
        )

    def order_columns(self):
        """
        Orders the columns of the pivot table based on the provided column list.
        """
        columns = (
            sorted(
                self.pivot.columns,
                key=lambda x: self.order.index(x[1])
                if x[1] in self.order
                else float("inf"),
            )
            if self.order
            else self.pivot.columns
        )
        self.pivot = self.pivot[columns]
        
    def correlate(self):
        """
        Calculates and stores correlation and p-value matrices.
        """
        self.pvalues = self.create_corr_matrix('pvalues')
        self.correlations = self.create_corr_matrix('correlations')
        self.corr_masked = self.correlations[self.pvalues < self.pvalue_threshold]
        

    def create_corr_matrix(self, result_type):
        """
        Creates a correlation matrix for either correlation values or p-values.

        Args:
            result_type (str): Type of result to return ('pvalues' or 'correlations').

        Returns:
            pd.DataFrame: A DataFrame containing the requested correlation matrix.
        """
        method = correlate(self.method, result_type)
        return self.pivot.corr(method=method, min_periods=self.n_minimum).loc[
            tuple([self.var1, self.var2])
        ]

    def find_missing_overlap(self):
        """
        Identifies and reports variable pairs with insufficient data overlap.
        """
        self.missing_overlap = [
            ((self.var1, index), (self.var2, column))
            for (index, column), is_na in self.correlations.T.isna().stack().items()
            if is_na
        ]
        if self.missing_overlap:
            print(
                f"{self.grouping} insuficient overlapp for {self.missing_overlap} pairs"
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

@dataclass
class Network:
    """
    A class to represent a network/graph constructed from a correlation matrix.

    Attributes:
    matrix (Matrix): An instance of the Matrix class containing the data and correlation matrix.

    Methods:
    max_node_degree(): Returns the maximum degree of the nodes in the graph.
    is_directed(): Property that checks if the network is directed (based on the matrix being square).
    get_title(): Generates a title for the network graph.
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
            self.is_directed = self.matrix.is_square
            self.G = nx.MultiDiGraph() if self.matrix.is_square else nx.Graph()
            # directed edge -  to_correlate[0] --> to_correlate[1]
            self.G.clear()

            self.G.add_nodes_from(
                self.matrix.corr_masked.columns.tolist()
            )  # adds every BR as a node
            self.edge_labels = {}
            for (row, col), correlation in self.matrix.corr_masked.stack().dropna().items():
                # Add edge to the graph with edge weight and color
                # Avoid self sorrelation
                if not(row == col and not self.is_directed):
                    self.G.add_edge(
                        row,
                        col,
                        weight=abs(correlation),
                        color="red" if correlation > 0 else "blue",
                    )
                    
                    self.edge_labels[(row, col)] = f"{correlation:.2f}"

            angles = np.linspace(
                0, 2 * np.pi, len(self.matrix.corr_masked.columns), endpoint=False
            )
            self.pos = {
                col: (np.cos(angles[i]), np.sin(angles[i]))
                for i, col in enumerate(self.matrix.corr_masked.columns)
            }

            self.total_edges, self.pos_edges, self.neg_edges = self.edge_count()
            self.density = self.calculate_graph_density()
            self.max_degree, self.average_degree = self.calculate_node_degree()
            self.avg_clust_coeff_unweighted, self.avg_clust_coeff_weighted = self.calculate_clustering_coefficient()

            # self.local_efficiency = self.calculate_local_efficiency()
            # self.global_efficiency = self.calculate_global_efficiency()
            # self.characteristic_path_length = self.calculate_characteristic_path_length()

    def get_title(self):
        """Generates a formatted title for the network graph."""
        title = self.matrix.get_title()
        return title.replace('-', '->') if self.is_directed else title
    
    def edge_count(self):
        """
        Returns the total number of edges, positive edges, and negative edges in the graph.

        Returns:
            total_edges (int): The total number of edges in the graph.
            pos_edges (int): The number of edges with positive weights.
            neg_edges (int): The number of edges with negative weights.
        """
        total_edges = self.G.number_of_edges()
        pos_edges = 0
        neg_edges = 0

        for u, v, data in self.G.edges(data=True):
            color = data.get('color')  
            if color == 'red':
                pos_edges += 1
            elif color == 'blue':
                neg_edges += 1
        return total_edges, pos_edges, neg_edges
    
    def calculate_node_degree(self):
        """
        Returns:
          max_degree(int): the maximum degree of the nodes in the graph.
        """
        # Calculate degrees for all nodes and find the maximum and mean
        degrees = dict(self.G.degree())
        max_degree = max(degrees.values())
        average_degree = np.mean(list(degrees.values()))
        return max_degree, average_degree


    def calculate_graph_density(self):
        """
        Returns:
           graph_density (float): The density of the graph; edges/all_possible_edges.
        """
        num_edges = self.G.number_of_edges()
        num_nodes = self.G.number_of_nodes()

        if self.is_directed: #directed graph have doubble possible edges
            max_edges = num_nodes * (num_nodes - 1)
        else:
            max_edges = num_nodes * (num_nodes - 1) / 2
        return num_edges / max_edges if max_edges > 0 else 0
    

    
    def calculate_clustering_coefficient(self):
        """
        Calculates the average unweighted and weighted clustering coefficients for the graph.

        Returns:
            avg_clust_coeff_unweighted (float): The average unweighted clustering coefficient of the graph.
            avg_clust_coeff_weighted (float): The average weighted clustering coefficient of the graph.
        """
        if self.is_directed:
            # For directed graphs, use nx.clustering with 'directed' and 'weight' parameters
            clust_coeff_unweighted = None #nx.clustering(self.G.to_undirected(), weight=None)  # Unweighted
            clust_coeff_weighted = None #nx.clustering(self.G.to_undirected(), weight='weight')  # Weighted
        else:
            # For undirected graphs, use nx.clustering without 'directed' parameter
            clust_coeff_unweighted = nx.clustering(self.G)  # Unweighted
            clust_coeff_weighted = nx.clustering(self.G, weight='weight')  # Weighted

        avg_clust_coeff_unweighted = sum(clust_coeff_unweighted.values()) / len(clust_coeff_unweighted) if clust_coeff_unweighted else 0
        avg_clust_coeff_weighted = sum(clust_coeff_weighted.values()) / len(clust_coeff_weighted) if clust_coeff_weighted else 0

        return avg_clust_coeff_unweighted, avg_clust_coeff_weighted
    
    # def calculate_local_efficiency(self):
    #     """
    #     Returns:
    #         avg_local_eff_unweighted (float):  average unweighted local efficiency for i nodes.
    #         avg_local_eff_weighted (float):  average weighted local efficiency for i nodes.
    #     """
    #     local_eff_unweighted = 0
    #     local_eff_weighted = 0
    #     total_nodes = len(self.G.nodes())

    #     for node in self.G.nodes():
    #         # Determine neighbors for directed and undirected graphs
    #         if self.is_directed:
    #             neighbors = list(set(nx.predecessors(self.G, node)) | set(nx.successors(self.G, node)))
    #         else:
    #             neighbors = list(nx.neighbors(self.G, node))

    #         if len(neighbors) > 1:
    #             subgraph = self.G.subgraph(neighbors)
    #             # Calculate unweighted local efficiency
    #             local_eff_unweighted += nx.global_efficiency(subgraph)

    #             # Calculate weighted local efficiency
    #             # Ensure weights are considered in the subgraph efficiency calculation
    #             local_eff_weighted += nx.global_efficiency(subgraph, weight='weight')

    #     # Calculate average efficiencies
    #     avg_local_eff_unweighted = local_eff_unweighted / total_nodes if total_nodes > 0 else 0
    #     avg_local_eff_weighted = local_eff_weighted / total_nodes if total_nodes > 0 else 0
    #     return avg_local_eff_unweighted, avg_local_eff_weighted

    # def calculate_global_efficiency(self):
    #     """
    #     Calculates the average unweighted and weighted global efficiency of the graph.

    #     Returns:
    #         avg_global_eff_unweighted (float): The average unweighted global efficiency of the graph.
    #         avg_global_eff_weighted (float): The average weighted global efficiency of the graph.
    #     """
    #     # Unweighted Global Efficiency
    #     avg_global_eff_unweighted = nx.global_efficiency(self.G)

    #     # Weighted Global Efficiency
    #     # Inverting weights for efficiency calculation as smaller weights imply stronger connections
    #     G_copy = self.G.copy()
    #     inverted_weights = {(u, v): 1 / data['weight'] for u, v, data in G_copy.edges(data=True)}
    #     nx.set_edge_attributes(G_copy, inverted_weights, 'inverted_weight')
    #     avg_global_eff_weighted = nx.global_efficiency(nx.stochastic_graph(G_copy, weight='inverted_weight'))


    #     return avg_global_eff_unweighted, avg_global_eff_weighted

    

    # def calculate_characteristic_path_length(self):
    #     """
    #     Calculates the average unweighted and weighted characteristic path length of the graph.

    #     Returns:
    #         avg_path_length_unweighted (float): The average unweighted characteristic path length of the graph.
    #         avg_path_length_weighted (float): The average weighted characteristic path length of the graph.
    #     """
    #     if nx.is_connected(self.G):
    #         avg_path_length_unweighted = nx.average_shortest_path_length(self.G)
            
    #         # Inverting weights for path length calculation as smaller weights imply stronger connections
    #         G_copy = self.G.copy()
    #         inverted_weights = {(u, v): 1 / data['weight'] for u, v, data in G_copy.edges(data=True)}
    #         nx.set_edge_attributes(G_copy, inverted_weights, 'inverted_weight')
    #         avg_path_length_weighted = nx.average_shortest_path_length(G_copy, weight='inverted_weight')
    #     else:
    #         avg_path_length_unweighted = None
    #         avg_path_length_weighted = None

    #     return avg_path_length_unweighted, avg_path_length_weighted
    