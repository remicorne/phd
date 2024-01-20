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
from module.correlogram import corrSelector
from module.Matrix import Matrices
from module.Network import Network
from module.getters import getCompoundAndRatiosDf

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
    compound_and_ratios_df = getCompoundAndRatiosDf(filename)
    matrices = Matrices(
        data=compound_and_ratios_df,
        group_by = 'treatment',
        between =correlogram_type, 
        variables = to_correlate.split('-'),
        accross="compound" if correlogram_type == "region" else "region",
        sub_selector ={"experiment":experiment},
        columns=columns,  # ordered to plot
        method=corr_method,  # 'pearson' 'spearman' 'kendall'
        pvalue_threshold=p_value_threshold,
        n_minimum=n_minimum,
    ).matrices

    fig, axs = generate_figure(matrices)
    networks = parallel_process(Network, [(matrix,) for matrix in matrices])
    
    [
        network_object.plot_ax(ax)
        for ax, network_object in zip(axs, networks)
    ]

    return fig


# PLOTTER FUNCTIONS
    
