import numpy as np
import networkx as nx
from module.Matrix import Matrix
from dataclasses import dataclass


@dataclass
class Network:
    matrix: Matrix

    def __post_init__(self):
        
        if self.matrix.is_square:
            # directed edge -  to_correlate[0] --> to_correlate[1]
            self.G = nx.MultiDiGraph()
        else:
            self.G = nx.Graph()
        self.G.clear()

        self.G.add_nodes_from(
            self.matrix.corr_masked.columns.tolist()
        )  # adds every BR as a node
        row_indices, col_indices = np.where(self.matrix.corr_masked.notna())
        for row_i, col_i in zip(row_indices, col_indices):
            correlation = self.matrix.corr_masked.iloc[row_i, col_i]
            row, col = self.matrix.corr_masked.index[row_i], self.matrix.corr_masked.columns[col_i]
            edge_color = "black"
            if self.matrix.is_square:
                edge_color = "red" if correlation > 0 else "blue"
            # Add edge to the graph with edge weight and color
            self.G.add_edge(
                row,
                col,
                weight=abs(correlation),
                color=edge_color,
                label=f"{correlation:.2f}",
            )

            ##### Draw the graph
        # pos = nx.spring_layout(G, seed=42)  # using a seed for consistency need allensdk working
        self.edges = self.G.edges()
        self.weights = list(nx.get_edge_attributes(self.G, "weight").values())
        self.edge_colors = list(nx.get_edge_attributes(self.G, "color").values())

        # Create a custom circular layout based on column order #FIXME
        column_order = list(self.matrix.corr_masked.columns)
        num_nodes = len(column_order)
        angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
        self.pos = {
            col: (np.cos(angles[i]), np.sin(angles[i]))
            for i, col in enumerate(column_order)
        }
