
# def graph_draw(graph, node_names=None, color_by_y=False, vertex_size=50, edge_marker_size=10, edge_font_size=10, edge_pen_width=2.2, vertex_font_size=10, output_size=(400, 400)):
#     # Provided data
#     links = getattr(graph, 'edge_index', None)
#     weights = getattr(graph, 'edge_attr', None)
#     y = getattr(graph, 'y', None)

#     # Create a new graph using graph_tool
#     g = gt.Graph(directed=True)

#     # Add vertices
#     v_list = []
#     for _ in range(links.max().item() + 1):
#         v_list.append(g.add_vertex())

#     # Add edges and set edge weights
#     if weights is not None:
#         e_weight = g.new_edge_property("double")  # Create a double (float) edge weight property map
#         for i in range(links.shape[1]):
#             source, target = links[:, i]
#             e = g.add_edge(v_list[source], v_list[target])
#             e_weight[e] = weights[i]
#     else: 
#         e_weight = None

#     # Convert y values to colors if color_by_y is True
#     if color_by_y:
#         norm = plt.Normalize(y.min().item(), y.max().item())  # Normalize y values between its min and max
#         colormap = mpl.colormaps['bwr']
#         vertex_fill_colors = [colormap(norm(value.item())) for value in y]
#     else:
#         vertex_fill_colors = [(0.5, 0.5, 0.5, 1) for _ in y]  # default gray color

#     v_fill_color_prop = g.new_vertex_property("vector<double>")
#     for i, v in enumerate(g.vertices()):
#         v_fill_color_prop[v] = vertex_fill_colors[i]

#     # Prepare vertex text (node names) if node_names is provided
#     v_text_prop = g.new_vertex_property("string")
#     if node_names:
#         for i, v in enumerate(g.vertices()):
#             v_text_prop[v] = node_names[i]

#     # Draw the graph with adjusted aesthetics
#     pos = gt.sfdp_layout(g, max_iter=0)  # Use force-directed layout

#     gt.graph_draw(g, pos, 
#                 vertex_fill_color=v_fill_color_prop, vertex_size=vertex_size, vertex_pen_width=1.5,
#                 edge_marker_size=edge_marker_size, edge_text=e_weight, edge_font_size=edge_font_size, edge_pen_width=edge_pen_width, 
#                 vertex_font_size=vertex_font_size, vertex_text=v_text_prop if node_names else None, output_size=output_size)

                    
    