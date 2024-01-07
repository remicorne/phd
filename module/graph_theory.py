### AllenSDK
# #modulenotfounderror?
# from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi


import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from module.utils import (
    flatten,
    figure_cache,
    generate_figure,
    parallel_process,
)
from module.correlogram import corrSelector, buildExperimentCorrmatrices
from module.Network import Network

########## GENERATE NETWORKS
# network for SINGLE correlation matrix (all weighted --- symetric: undirected -- unsymetric: multiedge and directed)  #DONE
# METTA network made of MULTIPLE correlation matrices (multiedge -- directed -- weighted)

########## COMPARE NETWORKS
# graph metrics : average edge weight, average degree, clustering coefficient, degree distribution, average shortest path length (stats where possible)
# network alignment : GRAAL / NetAlign / Jaccard Index


# USER FUNCTIONS
@figure_cache("network")
def network(
    filename,
    experiment=None,
    correlogram_type=None,
    to_correlate=None,
    p_value_threshold=0.05,
    n_minimum=5,
    columns=None,
    from_scratch=None,
    corr_method="pearson",  # HARD CODE FIXME
    # hierarchical_clustering=None,
):
    """
    This is the function that is called by the user, it will call buildExperimentalCorrelogram that builds a signle correlogram
    The function can be called with parameters, user input will be required if not
    """
    (
        filename,
        experiment,
        correlogram_type,
        to_correlate,
        p_value_threshold,
        n_minimum,
        columns,
        corr_method,
        from_scratch,
    ) = corrSelector(  # generic prompter for selecting corr matrices, probably need to add pearson/spearman
        filename,
        experiment=experiment,
        correlogram_type=correlogram_type,
        to_correlate=to_correlate,
        p_value_threshold=p_value_threshold,
        n_minimum=n_minimum,
        columns=columns,
        corr_method=corr_method,
        from_scratch=from_scratch,
    )
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate

    matrices = buildExperimentCorrmatrices(
        filename,
        experiment,
        correlogram_type,  # compound / ratio
        to_correlate,  # whihc compound/s / ratio/s
        p_value_threshold,
        n_minimum,
        columns,  # ordered to plot
        corr_method,  # 'pearson' 'spearman' 'kendall'
        from_scratch,  # Used in decorator
    )

    fig, axs = generate_figure(matrices)
    networks = parallel_process(Network, matrices)
    
    [
        plotNetwork(matrix, ax, network_object)
        for matrix, ax, network_object in zip(matrices, axs, networks)
    ]

    return fig


# PLOTTER FUNCTIONS
    

def plotNetwork(matrix, ax, network):
    """
    Network/Graph plotter for single correlation matrix ~ to be fed to plotExperiment()
    input: a single element from matricies i.e. for one treatment
    output:  ax with graph plotted
    """

    # Draw nodes and edges
    nx.draw_networkx_nodes(
        network.G,
        network.pos,
        node_size=1100,
        alpha=0.95,
        node_color="white",
        edgecolors="black",
        ax=ax,
    )
    nx.draw_networkx_edges(
        network.G,
        network.pos,
        edgelist=network.edges,
        width=network.weights,
        edge_color=network.edge_colors,
        ax=ax,
        node_size=1100,
        arrowstyle="->",
        arrowsize=20,
    )
    # Add labels to nodes
    node_labels = {node: node for node in network.G.nodes()}  # Label nodes with their names
    nx.draw_networkx_labels(network.G, network.pos, labels=node_labels, font_size=18, ax=ax)

    # Set title for the graph
    ax.set_frame_on(False)
    ax.set_title(
        f"{'-'.join(matrix.variables)} in {matrix.treatment}", fontsize=28, pad=-10, y=1
    )
