import graph_tool.all as gt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import torch
import matplotlib.cm as cm


## 1. 

def _extract_graph_components_uduw(graph):
    """
    Extract important components of a PyTorch Geometric graph. 
    This function is tailored for undirected, unweighted graphs.

    Parameters:
        graph (torch_geometric.data.Data): Input graph in PyTorch Geometric format.

    Returns:
        links (torch.Tensor): Tensor representing graph edges.
        num_nodes (int): Number of nodes in the graph.
        x (torch.Tensor, optional): Node feature matrix.
        y (torch.Tensor, optional): Node labels or target values.
    """
    # Extract unique edges (links)
    unique_edges = set(tuple(sorted(edge)) for edge in graph.edge_index .t().tolist()) 
    links = torch.tensor(list(unique_edges)).t().long()

    # Extract number of nodes
    num_nodes = graph.num_nodes

    # Extract node features (x) and labels/targets (y)
    x = torch.tensor(graph.x) if (hasattr(graph, 'x') and graph.x is not None) else None
    y = torch.tensor(graph.y) if (hasattr(graph, 'y') and graph.y is not None) else None

    return links, num_nodes, x, y


def _extract_graph_components_udw(graph):
    """
    Extract important components of a PyTorch Geometric graph. 
    This function is tailored for undirected, weighted graphs.

    Parameters:
        graph (torch_geometric.data.Data): Input graph in PyTorch Geometric format.

    Returns:
        links (torch.Tensor): Tensor representing graph edges.
        num_nodes (int): Number of nodes in the graph.
        x (torch.Tensor, optional): Node feature matrix.
        y (torch.Tensor, optional): Node labels or target values.
    """
    # Extract unique edges (links)
    unique_edges_dict = {tuple(e.tolist()): w.item() for e, w in zip(graph.edge_index.sort(axis=0)[0].t(), graph.edge_attr)}
    links = torch.stack([torch.tensor(k) for k in unique_edges_dict]).t().long()

    # Extract weights 
    weights = torch.tensor([v for v in unique_edges_dict.values()]).float()

    # Extract number of nodes
    num_nodes = graph.num_nodes

    # Extract node features (x) and labels/targets (y)
    x = torch.tensor(graph.x) if (hasattr(graph, 'x') and graph.x is not None) else None
    y = torch.tensor(graph.y) if (hasattr(graph, 'y') and graph.y is not None) else None

    return links, weights, num_nodes, x, y


## 2. 

def _set_nodes(gt_graph, num_nodes, node_names=None):
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

def _map_y_to_vertex_color(gt_graph, vertices, node_color):
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

## 3. links 

def _set_links_uduw(gt_graph,links):
    gt_graph.add_edge_list(links.t().tolist())
    return gt_graph

def _set_links_udw(gt_graph,vertices,links,weights,edge_weight_text_format,edge_weight_width_scale):
    #vertices = {i: gt_graph.add_vertex() for i in range(num_nodes)}
    e_weight = gt_graph.new_edge_property("string")
    e_pen_width = gt_graph.new_edge_property("double")
    for idx, (start, end) in enumerate(links.t().tolist()):
        e = gt_graph.add_edge(vertices[start], vertices[end])
        e_weight[e] = format(weights[idx].item(),edge_weight_text_format)
        e_pen_width[e] = (weights/weights.max())[idx].item() * edge_weight_width_scale
    return e_weight, e_pen_width, gt_graph

## codes 

def plot_undirected_unweighted(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None,
        node_color=None,
        node_size=None,
    ):
    """
    Visualize an undirected and unweighted graph.

    Parameters:
        graph (torch_geometric.data.Data): Input graph in PyTorch Geometric format.
        node_names (list, optional): List of node names. Default is None.
        output_path (str, optional): File path to save the plot. Default is None (no saving).
        layout_options (dict, optional): Dictionary of layout options for the layout algorithm. Default is None.
        draw_options (dict, optional): Dictionary of drawing options for the graph drawing. Default is None.

    Returns:
        None
    """
    # 1. Prepare layout and drawing options as dictionaries.
    layout_options = layout_options or {}
    draw_options = draw_options or {}

    # 2. extract_graph_components.
    links, num_nodes, x, y = _extract_graph_components_uduw(graph)

    # 3. Create a graph_tool graph.
    gt_graph = gt.Graph(directed=False)  # Undirected graph

    # 4. Set nodes.
    vertices, v_text_prop, gt_graph = _set_nodes(gt_graph, num_nodes, node_names)

    # 5. Map colors and sizes based on y values.
    v_color, gt_graph = _map_y_to_vertex_color(gt_graph, vertices, node_color)

    # 6. Add nodes and edges to the graph.
    gt_graph = _set_links_uduw(gt_graph,links)

    # 7. Set draw_options.
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
    if node_color is not None: 
        draw_options['vertex_fill_color'] = v_color  # Set the vertex color based on y    

    # 8. Perform graph layout using sfdf_layout and draw the graph using graph_draw.
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, vertex_text=v_text_prop, **draw_options)


def plot_undirected_weighted(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None,
        node_color=None,
        node_size=False,        
        edge_weight_text=True,
        edge_weight_width=True,        
        edge_weight_text_format=".2f", 
        edge_weight_width_scale=1.0
    ):
    """
    Visualize an undirected and unweighted graph.

    Parameters:
        graph (torch_geometric.data.Data): Input graph in PyTorch Geometric format.
        node_names (list, optional): List of node names. Default is None.
        output_path (str, optional): File path to save the plot. Default is None (no saving).
        layout_options (dict, optional): Dictionary of layout options for the layout algorithm. Default is None.
        draw_options (dict, optional): Dictionary of drawing options for the graph drawing. Default is None.

    Returns:
        None
    """
    # 1. Prepare layout and drawing options as dictionaries.
    layout_options = layout_options or {}
    draw_options = draw_options or {}

    # 2. extract_graph_components.
    links, weights, num_nodes, x, y = _extract_graph_components_udw(graph)

    # 3. Create a graph_tool graph.
    gt_graph = gt.Graph(directed=False)  # Undirected graph

    # 4. Set nodes.
    vertices, v_text_prop, gt_graph = _set_nodes(gt_graph, num_nodes, node_names)

    # 5. Map colors and sizes based on y values.
    v_color, gt_graph = _map_y_to_vertex_color(gt_graph, vertices, node_color)

    # 6. Set links.
    e_weight, e_pen_width, gt_graph = _set_links_udw(gt_graph,vertices,links,weights,edge_weight_text_format,edge_weight_width_scale)

    # 7. Set draw_options.
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
    if node_color is not None: 
        draw_options['vertex_fill_color'] = v_color  # Set the vertex color based on y    
    if edge_weight_text: 
        draw_options['edge_text'] = e_weight  # Set edge text property
    if edge_weight_width: 
        draw_options['edge_pen_width'] = e_pen_width  # Use edge weight to adjust edge pen width

    # 8. Perform graph layout using sfdf_layout and draw the graph using graph_draw.
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, vertex_text=v_text_prop, **draw_options)
