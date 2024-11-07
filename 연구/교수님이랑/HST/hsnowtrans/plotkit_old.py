import graph_tool.all as gt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import torch

def plot_undirected_unweighted_without_y(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None
    ):
    """
    Visualize an undirected and unweighted graph without signal values.

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

    # 2. Create a graph_tool graph (undirected,unweighted) version 
    gt_graph = gt.Graph(directed=False)  # Undirected graph
    unique_edges = set(tuple(sorted(edge)) for edge in graph.edge_index .t().tolist()) # Remove duplicate edges and add unique edges to the graph
    links = torch.tensor(list(unique_edges)).t().long()
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


def plot_undirected_weighted_without_y(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None, 
        edge_weight_text=True,
        edge_weight_width=True,        
        edge_weight_text_format=".2f", 
        edge_weight_width_scale=1.0
    ):
    """
    Visualize a weighted undirected graph without signal values.

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

    # 2. Create a graph_tool graph -- (undirected,weighted) version
    gt_graph = gt.Graph(directed=False)  # Undirected graph
    unique_edges_dict = {tuple(e.tolist()): w.item() for e, w in zip(graph.edge_index.sort(axis=0)[0].t(), graph.edge_attr)}
    links = torch.stack([torch.tensor(k) for k in unique_edges_dict]).t().long()
    weights = torch.tensor([v for v in unique_edges_dict.values()]).float()
    num_nodes = graph.num_nodes

    # 3. Set node names if provided
    v_text_prop = None
    if node_names:
        v_text_prop = gt_graph.new_vertex_property("string")
        for v, name in enumerate(node_names):
            v_text_prop[v] = name

    # 5. Add edges and set edge weights -- weighted version 
    vertices = {i: gt_graph.add_vertex() for i in range(num_nodes)}
    e_weight = gt_graph.new_edge_property("string")
    e_pen_width = gt_graph.new_edge_property("double")
    for idx, (start, end) in enumerate(links.t().tolist()):
        e = gt_graph.add_edge(vertices[start], vertices[end])
        e_weight[e] = format(weights[idx].item(),edge_weight_text_format)
        e_pen_width[e] = (weights/weights.max())[idx].item() * edge_weight_width_scale

    # 6. Set draw_options
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
    if edge_weight_text: 
        draw_options['edge_text'] = e_weight  # Set edge text property
    if edge_weight_width: 
        draw_options['edge_pen_width'] = e_pen_width  # Use edge weight to adjust edge pen width

    # 7. Perform graph layout using sfdf_layout and draw the graph using graph_draw
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, **draw_options)


def plot_directed_unweighted_without_y(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None
    ):
    """
    Visualize a directed and unweighted graph without signal values.

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

def plot_directed_weighted_without_y(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None,
        edge_weight_text=True,
        edge_weight_width=True,
        edge_weight_text_format=".2f", 
        edge_weight_width_scale=1.0
    ):
    """
    Visualize a directed and unweighted graph without signal values.

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

    # 2. Create a graph_tool graph -- (directed,weighted) version 
    gt_graph = gt.Graph(directed=True)
    links = graph.edge_index
    weights = graph.edge_attr
    num_nodes = graph.num_nodes

    # 3. Set node names if provided
    v_text_prop = None
    if node_names:
        v_text_prop = gt_graph.new_vertex_property("string")
        for v, name in enumerate(node_names):
            v_text_prop[v] = name

    # 5. Add edges and set edge weights -- weighted version 
    vertices = {i: gt_graph.add_vertex() for i in range(num_nodes)}
    e_weight = gt_graph.new_edge_property("string")
    e_pen_width = gt_graph.new_edge_property("double")
    for idx, (start, end) in enumerate(links.t().tolist()):
        e = gt_graph.add_edge(vertices[start], vertices[end])
        e_weight[e] = format(weights[idx].item(),edge_weight_text_format)
        e_pen_width[e] = (weights/weights.max())[idx].item() * edge_weight_width_scale

    # 6. Set draw_options
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
    if edge_weight_text: 
        draw_options['edge_text'] = e_weight  # Set edge text property
    if edge_weight_width: 
        draw_options['edge_pen_width'] = e_pen_width  # Use edge weight to adjust edge pen width

    # 7. Perform graph layout using sfdf_layout and draw the graph using graph_draw
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, vertex_text=v_text_prop, **draw_options)

