import graph_tool.all as gt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import torch

def plot_directed_unweighted_without_y(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None
    ):
    """
    Visualize a directed and unweighted graph.

    Parameters:
        graph (torch_geometric.data.Data): Input graph in PyTorch Geometric format.
        node_names (list, optional): List of node names. Default is None.
        output_path (str, optional): File path to save the plot. Default is None (no saving).
        layout_options (dict, optional): Dictionary of layout options for the layout algorithm. Default is None.
        draw_options (dict, optional): Dictionary of drawing options for the graph drawing. Default is None.

    Returns:
        None
    """
    # 1. Prepare layout and drawing options as dictionaries
    layout_options = layout_options or {}
    draw_options = draw_options or {}

    # 2. Create a graph_tool graph -- (directed,unweighted) version 
    gt_graph = gt.Graph(directed=True)
    links = graph.edge_index
    num_nodes = graph.num_nodes

    # 3. Set node names if provided
    v_text_prop = None
    if node_names:
        v_text_prop = gt_graph.new_vertex_property("string")
        for v, name in enumerate(node_names):
            v_text_prop[v] = name

    # 5. Add nodes and edges to the graph -- unweighted version
    gt_graph.add_edge_list(links.t().tolist())

    # 6. Set draw_options
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
        
    # 7. Perform graph layout using sfdf_layout and draw the graph using graph_draw
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, vertex_text=v_text_prop, **draw_options)