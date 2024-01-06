### AllenSDK
# #modulenotfounderror?
# from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi


import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from module.utils import flatten, get_or_add, plotExperiment
from module.correlogram import corrSelector, buildExperimentCorrmatrices


########## GENERATE NETWORKS
# network for SINGLE correlation matrix (all weighted --- symetric: undirected -- unsymetric: multiedge and directed)  #DONE
# METTA network made of MULTIPLE correlation matrices (multiedge -- directed -- weighted)

########## COMPARE NETWORKS
#graph metrics : average edge weight, average degree, clustering coefficient, degree distribution, average shortest path length (stats where possible)
#network alignment : GRAAL / NetAlign / Jaccard Index


# USER FUNCTIONS
def network(
    filename,
    experiment=None,
    correlogram_type=None,
    to_correlate=None,
    p_value_threshold=0.05,
    n_minimum=5,
    columns=None,
    from_scratch=None,
    corr_method = 'pearson' #HARD CODE FIXME
    # hierarchical_clustering=None,
):
    """
    This is the function that is called by the user, it will call buildExperimentalCorrelogram that builds a signle correlogram
    The function can be called with parameters, user input will be required if not
    """
    filename,experiment,correlogram_type,to_correlate,p_value_threshold,n_minimum,columns,corr_method,from_scratch=corrSelector( # generic prompter for selecting corr matrices, probably need to add pearson/spearman
                                                                                                                                filename,
                                                                                                                                experiment=experiment,
                                                                                                                                correlogram_type=correlogram_type,
                                                                                                                                to_correlate=to_correlate,
                                                                                                                                p_value_threshold=p_value_threshold,
                                                                                                                                n_minimum=n_minimum,
                                                                                                                                columns=columns,
                                                                                                                                corr_method = corr_method,
                                                                                                                                from_scratch=from_scratch )
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate 

    plotExperimentalNetworks(
        filename,
        experiment=experiment,
        correlogram_type=correlogram_type,
        to_correlate=to_correlate,
        p_value_threshold=p_value_threshold,
        n_minimum=n_minimum,
        columns=columns,
        corr_method=corr_method,
        from_scratch=from_scratch,
        # hierarchical_clustering=hierarchical_clustering,
    )




# PLOTTER FUNCTIONS

@get_or_add("network")
def plotExperimentalNetworks( 
    filename,
    experiment,
    correlogram_type,
    to_correlate,
    p_value_threshold,
    n_minimum,
    columns,
    from_scratch,  # Used in decorator
    corr_method,
):
    matrices=buildExperimentCorrmatrices(
                                            filename,
                                            experiment,
                                            correlogram_type, #compound / ratio
                                            to_correlate, # whihc compound/s / ratio/s
                                            p_value_threshold,
                                            n_minimum,
                                            columns, #ordered to plot
                                            corr_method, # 'pearson' 'spearman' 'kendall'
                                            from_scratch,  # Used in decorator
                                            )

    fig = plotExperiment(matrices, plotNetwork) 

    return fig


def plotNetwork(matrix, ax): 
    '''
    Network/Graph plotter for single correlation matrix ~ to be fed to plotExperiment()
    input: a single element from matricies i.e. for one treatment
    output:  ax with graph plotted
     '''

    G=buildNetwork(matrix)
   
    ##### Draw the graph
    # pos = nx.spring_layout(G, seed=42)  # using a seed for consistency need allensdk working 
    edges = G.edges()
    weights = list(nx.get_edge_attributes(G, 'weight').values())
    edge_colors = list(nx.get_edge_attributes(G, 'color').values())

    # Create a custom circular layout based on column order #FIXME
    column_order = list(matrix.corr.columns)
    num_nodes = len(column_order)
    angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
    pos = {col: (np.cos(angles[i]), np.sin(angles[i])) for i, col in enumerate(column_order)}

    # Draw nodes and edges
    nx.draw_networkx_nodes(G, pos, node_size=1100, alpha=0.95, node_color='white', edgecolors='black', ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights, edge_color=edge_colors, ax=ax, node_size=1100, arrowstyle='->', arrowsize=20)
    # Add labels to nodes
    node_labels = {node: node for node in G.nodes()}  # Label nodes with their names
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=18, ax=ax)

    # Set title for the graph
    ax.set_frame_on(False)  
    ax.set_title(f"{'-'.join(matrix.variables)} in {matrix.treatment}", fontsize=28, pad=-10, y=1)
    
    return ax

# BUILDER FUNCTIONS 
def buildExperimentalNetworks(matrices):
    '''
    input: matrices i.e. df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate for each treatment in an experiment
    output: networks i.e. G, treatment, to_correlate for each treatment in an experiment
    '''
    networks=[]
    for (df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate) in matrices:
        G=buildNetwork(df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate)
        networks.append([G, treatment, to_correlate])
    return networks


def buildNetwork(matrix):
    '''
    input:  matrices[n] = df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate 
    --- updates graph_stats df --- #TODO
    output: network/graph for a single treatment G, directed or not based off length of to_correlate 
    '''
    if matrix.is_square:
        #directed edge -  to_correlate[0] --> to_correlate[1] 
        G=nx.MultiDiGraph()
    else:
        G = nx.Graph()
    G.clear()
        
    G.add_nodes_from(matrix.corr.columns.tolist()) #adds every BR as a node
    row_indices, col_indices = np.where(matrix.corr != np.nan)
    for row_i, col_i in zip(row_indices, col_indices):
        correlation = matrix.corr.iloc[row_i, col_i]
        row, col = matrix.corr.index[row_i], matrix.corr.columns[col_i]
        edge_color = 'black'
        if matrix.is_square:
            edge_color = 'red' if correlation > 0 else 'blue'
        # Add edge to the graph with edge weight and color
        G.add_edge(row, col, weight=abs(correlation), color=edge_color, label=f"{correlation:.2f}") 

    return G





