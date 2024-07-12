import numpy as np
import networkx as nx
from module.Matrix import Matrix
from dataclasses import dataclass


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

    def __post_init__(self):
        """
        Initializes the network/graph from the given matrix. Constructs a directed or undirected graph
        based on the matrix's properties and fills it with nodes and edges based on the correlation data.
        """
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

    @property
    def is_directed(self):
        """
        Checks if the network is directed, based on whether the matrix is square. Directed graph are use for double correlations i.e. 'GLU-DA'. 
        """
        return self.matrix.is_square
    
    @property  
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
            weight = data.get('weight', 1)  # Default weight to 1 if not specified
            if weight > 0:
                pos_edges += 1
            elif weight < 0:
                neg_edges += 1

        return total_edges, pos_edges, neg_edges

    @property
    def node_degree(self):
        """
        Returns:
          max_degree(int): the maximum degree of the nodes in the graph.
          average_degree(int): average node degree of graph.
        """
        # Calculate degrees for all nodes and find the maximum and mean
        degrees = dict(self.G.degree())
        max_degree = max(degrees.values())
        average_degree = np.mean(list(degrees.values()))
        return max_degree, average_degree

    @property
    def graph_density(self):
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
    
    @property
    def local_efficiency(self):
        """
        Returns:
            avg_local_eff_unweighted (float):  average unweighted local efficiency for i nodes.
            avg_local_eff_weighted (float):  average weighted local efficiency for i nodes.
        """
        local_eff_unweighted = 0
        local_eff_weighted = 0
        total_nodes = len(self.G.nodes())

        for node in self.G.nodes():
            # Determine neighbors for directed and undirected graphs
            if self.is_directed:
                neighbors = list(set(nx.predecessors(self.G, node)) | set(nx.successors(self.G, node)))
            else:
                neighbors = list(nx.neighbors(self.G, node))

            if len(neighbors) > 1:
                subgraph = self.G.subgraph(neighbors)
                # Calculate unweighted local efficiency
                local_eff_unweighted += nx.global_efficiency(subgraph)

                # Calculate weighted local efficiency
                # Ensure weights are considered in the subgraph efficiency calculation
                local_eff_weighted += nx.global_efficiency(subgraph, weight='weight')

        # Calculate average efficiencies
        avg_local_eff_unweighted = local_eff_unweighted / total_nodes if total_nodes > 0 else 0
        avg_local_eff_weighted = local_eff_weighted / total_nodes if total_nodes > 0 else 0
        return avg_local_eff_unweighted, avg_local_eff_weighted
    
    @property 
    def global_efficiency(self):
        """
        Calculates the average unweighted and weighted global efficiency of the graph.

        Returns:
            avg_global_eff_unweighted (float): The average unweighted global efficiency of the graph.
            avg_global_eff_weighted (float): The average weighted global efficiency of the graph.
        """
        # Unweighted Global Efficiency
        avg_global_eff_unweighted = nx.global_efficiency(self.G)

        # Weighted Global Efficiency
        # Inverting weights for efficiency calculation as smaller weights imply stronger connections
        G_copy = self.G.copy()
        inverted_weights = {(u, v): 1 / data['weight'] for u, v, data in G_copy.edges(data=True)}
        nx.set_edge_attributes(G_copy, inverted_weights, 'inverted_weight')
        avg_global_eff_weighted = nx.global_efficiency(nx.stochastic_graph(G_copy, weight='inverted_weight'))


        return avg_global_eff_unweighted, avg_global_eff_weighted
    
    @property  
    def clustering_coefficient(self):
        """
        Calculates the average unweighted and weighted clustering coefficients for the graph.

        Returns:
            avg_clust_coeff_unweighted (float): The average unweighted clustering coefficient of the graph.
            avg_clust_coeff_weighted (float): The average weighted clustering coefficient of the graph.
        """
        if self.is_directed:
            # For directed graphs, use nx.clustering with 'directed' and 'weight' parameters
            clust_coeff_unweighted = nx.clustering(self.G, weight=None)  # Unweighted
            clust_coeff_weighted = nx.clustering(self.G, weight='weight')  # Weighted
        else:
            # For undirected graphs, use nx.clustering without 'directed' parameter
            clust_coeff_unweighted = nx.clustering(self.G)  # Unweighted
            clust_coeff_weighted = nx.clustering(self.G, weight='weight')  # Weighted

        avg_clust_coeff_unweighted = sum(clust_coeff_unweighted.values()) / len(clust_coeff_unweighted)
        avg_clust_coeff_weighted = sum(clust_coeff_weighted.values()) / len(clust_coeff_weighted)

        return avg_clust_coeff_unweighted, avg_clust_coeff_weighted
    
    @property
    def characteristic_path_length(self):
        """
        Calculates the average unweighted and weighted characteristic path length of the graph.

        Returns:
            avg_path_length_unweighted (float): The average unweighted characteristic path length of the graph.
            avg_path_length_weighted (float): The average weighted characteristic path length of the graph.
        """
        if nx.is_connected(self.G):
            avg_path_length_unweighted = nx.average_shortest_path_length(self.G)
            
            # Inverting weights for path length calculation as smaller weights imply stronger connections
            G_copy = self.G.copy()
            inverted_weights = {(u, v): 1 / data['weight'] for u, v, data in G_copy.edges(data=True)}
            nx.set_edge_attributes(G_copy, inverted_weights, 'inverted_weight')
            avg_path_length_weighted = nx.average_shortest_path_length(G_copy, weight='inverted_weight')
        else:
            avg_path_length_unweighted = None
            avg_path_length_weighted = None

        return avg_path_length_unweighted, avg_path_length_weighted
    
    @property
    def get_title(self):
        """Generates a formatted title for the network graph."""
        title = self.matrix.get_title()
        return title.replace('-', '->') if self.is_directed else title
    
    @property
    def network_id(self):
        group = self.matrix.grouping

    def plot_ax(self, ax):
        """
        Plots the network graph on the provided matplotlib axis.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the graph.
        """
        nx.draw_networkx_nodes(
            self.G,
            self.pos,
            node_size=1100,
            alpha=0.95,
            node_color="white",
            edgecolors="black",
            ax=ax,
        )
        nx.draw_networkx_edges(
            self.G,
            self.pos,
            width=list(nx.get_edge_attributes(self.G, "weight").values()),
            edge_color=list(nx.get_edge_attributes(self.G, "color").values()),
            ax=ax,
            node_size=1100,
            **(
                {"arrowstyle": "->", "arrowsize": 20}
                if self.is_directed
                else {}
            ),
        )
        # Add labels to nodes
        node_labels = {
            node: node for node in self.G.nodes()
        } 
        nx.draw_networkx_labels(
            self.G, self.pos, labels=node_labels, font_size=18, ax=ax
        )    
        # nx.draw_networkx_edge_labels(
        #     self.G, self.pos, edge_labels=self.edge_labels, font_size=18, ax=ax
        # )

        # Set title for the graph
        ax.set_frame_on(False)
        ax.set_title(self.get_title(), fontsize=28, pad=-10, y=1)
        self.matrix.is_square
