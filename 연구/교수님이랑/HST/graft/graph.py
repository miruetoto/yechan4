from .graph_extract import * 
from .graph_nodes import * 
from .graph_links import * 

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
    links, num_nodes, x, y = extract_graph_components_uduw(graph)

    # 3. Create a graph_tool graph.
    gt_graph = gt.Graph(directed=False)  # Undirected graph

    # 4. Set nodes.
    vertices, v_text_prop, gt_graph = set_nodes(gt_graph, num_nodes, node_names)

    # 5. Map colors and sizes.
    v_color, gt_graph = map_vertex_color(gt_graph, vertices, node_color)
    v_size, gt_graph = map_vertex_size(gt_graph, vertices, node_size)

    # 6. Add nodes and edges to the graph.
    gt_graph = set_links_uduw(gt_graph,links)

    # 7. Set draw_options.
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
    if node_color is not None:
        draw_options['vertex_fill_color'] = v_color  # Set the vertex color based on y
    if node_size is not None:
        draw_options['vertex_size'] = v_size  # Set the vertex size based on node_size

    # 8. Perform graph layout using sfdf_layout and draw the graph using graph_draw.
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, vertex_text=v_text_prop, **draw_options)


def plot_undirected_weighted(
        graph, 
        node_names=None, 
        layout_options=None, 
        draw_options=None,
        node_color=None,
        node_size=None,        
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
    links, weights, num_nodes, x, y = extract_graph_components_udw(graph)

    # 3. Create a graph_tool graph.
    gt_graph = gt.Graph(directed=False)  # Undirected graph

    # 4. Set nodes.
    vertices, v_text_prop, gt_graph = set_nodes(gt_graph, num_nodes, node_names)

    # 5. Map colors and sizes.
    v_color, gt_graph = map_vertex_color(gt_graph, vertices, node_color)
    v_size, gt_graph = map_vertex_size(gt_graph, vertices, node_size)    

    # 6. Set links.
    e_weight, e_pen_width, gt_graph = set_links_udw(gt_graph,vertices,links,weights,edge_weight_text_format,edge_weight_width_scale)

    # 7. Set draw_options.
    draw_options.setdefault('output_size', (150 + num_nodes, 150 + num_nodes)) # Set default output size if not provided in draw_options
    if node_color is not None:
        draw_options['vertex_fill_color'] = v_color  # Set the vertex color based on y
    if node_size is not None:
        draw_options['vertex_size'] = v_size  # Set the vertex size based on node_size
    if edge_weight_text: 
        draw_options['edge_text'] = e_weight  # Set edge text property
    if edge_weight_width: 
        draw_options['edge_pen_width'] = e_pen_width  # Use edge weight to adjust edge pen width

    # 8. Perform graph layout using sfdf_layout and draw the graph using graph_draw.
    pos = gt.sfdp_layout(gt_graph, **layout_options)
    gt.graph_draw(gt_graph, pos=pos, vertex_text=v_text_prop, **draw_options)