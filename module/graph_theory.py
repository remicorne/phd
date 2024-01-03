### AllenSDK
# #modulenotfounderror?
# from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi


import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from module.utils import flatten, get_or_add
from module.correlogram import corrSelector, buildExperimentCorrMatricies

#https://networkx.org/documentation/stable/tutorial.html


########## GENERATE GRAPHS

# graph for correlations between the same compound 


#FUNCTIONS
def graph(
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
    filename,experiment,correlogram_type,to_correlate,p_value_threshold,n_minimum,columns,corr_method,from_scratch=corrSelector( # generic prompter for selecting corr matricies, probably need to add pearson/spearman
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

    buildExperimentalGraph(
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


@get_or_add("graph")
def buildExperimentalGraph( #TODO change the name as it should be experimental not single stupid remi
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
    matricies=buildExperimentCorrMatricies(
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

    fig = plotGraphs(matricies) 

    return fig


#builds figure for subplots loops on 
def plotGraphs(matricies):
    #JJB might be nice to have one color bar for all figures
    fig, axs = plt.subplots(2, 2, figsize=(22, 22))
    axs = flatten(axs)  # Put all the axes into a list at the same level

    # for matrix in matricies:
    for (df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate), ax in zip(
        matricies, axs
    ):
        print (f" treatment {treatment} correlating {to_correlate} ")

        plotGraph(df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate, ax)

    fig.tight_layout()
    return fig

#generates graph from matricies[0] - updates graph_stats df - plots graphs
def plotGraph(df_to_corr, correlation_matrix, T_F_mask_matrix, treatment, to_correlate, ax):
    
    ##### Build graph
    if len(to_correlate)>1:
        G=nx.MultiDiGraph()
    else:
        G = nx.Graph()
    G.clear()
        
    G.add_nodes_from(correlation_matrix.columns.tolist()) #adds every BR as a node
    
    # Iterate through the correlation matrix and significance matrix
    for i, col in enumerate(correlation_matrix.columns):
        for j, row in enumerate(correlation_matrix.index):
            if i < j:
                correlation = correlation_matrix.iloc[j, i]  # As correlation_matrix is symmetric
                significance = T_F_mask_matrix.iloc[j, i]

                # create edge where mask == False (significant correlations)
                if significance == False:
                    edge_color = 'red' if correlation > 0 else 'blue'
                    weight = abs(correlation)

                    # Add edge to the graph with edge weight and color
                    if len(to_correlate)>1:
                        G.add_edge(col, row, weight=abs(correlation), color=edge_color, label=f"{correlation:.2f}") #add only nodes for BRs with edge 
                    else:
                        G.add_edge(col, row, weight=abs(correlation), color=edge_color)

        ##### Draw the graph
        # pos = nx.spring_layout(G, seed=42)  # using a seed for consistency need allensdk working 
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        edge_colors = [G[u][v]['color'] for u, v in edges]

        # Create a custom circular layout based on column order #FIXME
        column_order = list(correlation_matrix.columns)
        num_nodes = len(column_order)
        angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
        pos = {col: (np.cos(angles[i]), np.sin(angles[i])) for i, col in enumerate(column_order)}

        # Draw nodes and edges
        nx.draw_networkx_nodes(G, pos, node_size=1100, alpha=0.95, node_color='white', edgecolors='black', ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights, edge_color=edge_colors, ax=ax)
        # Add labels to nodes
        node_labels = {node: node for node in G.nodes()}  # Label nodes with their names
        nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=18, ax=ax)

        # Set title for the graph
        ax.set_frame_on(False)  
        ax.set_title(f"{'-'.join(to_correlate)} in {treatment}", fontsize=28, pad=-10, y=1)
        
        return ax



#single graph multiple edge types corrisponding to each compound correlation DA-DA, DA-GLU

#G=nx.MultiGraph()

########## COMPARE GRAPHS

#graph metrics : average edge weight, average degree, clustering coefficient, degree distribution, average shortest path length (stats where possible)

#network alignment : GRAAL or NetAlign




