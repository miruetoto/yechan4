import graph_tool.all as gt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import torch
import matplotlib.cm as cm

def set_nodes(gt_graph, num_nodes, node_names=None):
    """
    Set node names for a graph_tool graph if provided.

    Parameters:
        gt_graph (graph_tool.Graph): Input graph in graph_tool format.
        num_nodes (int): Number of nodes in the graph.
        node_names (list, optional): List of node names. Default is None.

    Returns:
        v_text_prop (graph_tool.VertexPropertyMap): Vertex property map with node names.
    """
    v_text_prop = None
    if node_names:
        v_text_prop = gt_graph.new_vertex_property("string")
        for v, name in enumerate(node_names):
            v_text_prop[v] = name
    vertices = [gt_graph.add_vertex() for _ in range(num_nodes)]            
    return vertices, v_text_prop, gt_graph

def map_vertex_color(gt_graph, vertices, node_color):
    """
    Map y values to vertex colors for a graph_tool graph.

    Parameters:
        gt_graph (graph_tool.Graph): Input graph in graph_tool format.
        num_nodes (int): Number of nodes in the graph.
        y (torch.Tensor): Tensor with values to map to colors.

    Returns:
        v_color (graph_tool.VertexPropertyMap): Vertex property map with colors.
    """

    if node_color is None:
        return None, gt_graph
    else: 
        y = torch.tensor(node_color)   
        v_color = gt_graph.new_vertex_property("vector<double>")

        # Continuous values or too many unique values
        if y.dtype in (torch.float16, torch.float32, torch.float64) or len(y) > 10:
            y_min, y_max = y.min().item(), y.max().item()
            colormap = mpl.cm.get_cmap('spring')
            for idx, value in enumerate(y):
                normalized_value = (value.item() - y_min) / (y_max - y_min)  # Normalize between [0, 1]
                rgba = list(colormap(normalized_value))  # Convert to RGBA color
                v_color[vertices[idx]] = rgba

        # Categorical or discrete values
        else:
            colors_dict = {
                1: ['#F8766D'],
                2: ['#F8766D', '#00BFC4'],
                3: ['#F8766D', '#00BA38', '#619CFF'],
                4: ['#F8766D', '#7CAE00', '#00BFC4', '#C77CFF'],
                5: ['#F8766D', '#A3A500', '#00BF7D', '#00B0F6', '#E76BF3'],
                6: ['#F8766D', '#B79F00', '#00BA38', '#00BFC4', '#619CFF', '#F564E3'],
                7: ['#F8766D', '#C49A00', '#53B400', '#00C094', '#00B6EB', '#A58AFF', '#FB61D7'],
                8: ['#F8766D', '#CD9600', '#7CAE00', '#00BE67', '#00BFC4', '#00A9FF', '#C77CFF', '#FF61CC'],
                9: ['#F8766D', '#D39200', '#93AA00', '#00BA38', '#00C19F', '#00B9E3', '#619CFF', '#DB72FB', '#FF61C3'],
                10: ['#F8766D', '#A3A500', '#39B600', '#00BF7D', '#00BFC4', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC', '#D89000']
            }
            ggplot_colors = colors_dict[len(set(y.tolist()))]

            def hex_to_rgb_normalized(hex_color):
                rgb = mpl.colors.hex2color(hex_color)  # Gives RGB values between 0 and 1
                return tuple([float(val) for val in rgb])

            ggplot_colors_rgb = [hex_to_rgb_normalized(color) for color in ggplot_colors]
            unique_y = torch.unique(y)
            y_to_color = {int(val.item()): ggplot_colors_rgb[i % len(ggplot_colors_rgb)] for i, val in enumerate(unique_y)}
            for idx, value in enumerate(y):
                rgb = y_to_color[int(value.item())]
                v_color[vertices[idx]] = rgb

        return v_color, gt_graph

def map_vertex_size(gt_graph, vertices, node_size):
    """
    Map y values to vertex sizes for a graph_tool graph.

    Parameters:
        gt_graph (graph_tool.Graph): Input graph in graph_tool format.
        vertices (list): List of vertex indices.
        node_size (torch.Tensor): Tensor with values to map to sizes.

    Returns:
        v_size (graph_tool.VertexPropertyMap): Vertex property map with sizes.
    """

    if node_size is None:
        return None, gt_graph
    else: 
        y = torch.tensor(node_size)
        v_size = gt_graph.new_vertex_property("double")

        y_min, y_max = y.min().item(), y.max().item()
        min_size, max_size = 5, 20  # Define the range of vertex sizes
        for idx, value in enumerate(y):
            normalized_value = (value.item() - y_min) / (y_max - y_min)  # Normalize between [0, 1]
            size = min_size + normalized_value * (max_size - min_size)  # Map to size range
            v_size[vertices[idx]] = size

        return v_size, gt_graph
