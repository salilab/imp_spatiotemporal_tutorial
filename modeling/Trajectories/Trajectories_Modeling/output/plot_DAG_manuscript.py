"""
This script replots the directed acyclic graph of the trajectories model, as shown in Figure 5 of the manuscript.
Specifically, the colors of lettering within the draw_dag_in_graphviz function are modified.
"""
import IMP.spatiotemporal as spatiotemporal
import os
import shutil
from matplotlib import colormaps as cm
from matplotlib import colors as clr
import numpy as np
from graphviz import Digraph


# Rendering DAG
def draw_dag_in_graphviz(nodes, coloring=None, draw_label=True,
                         fontname="Helvetica", fontsize="18", penscale=0.6,
                         arrowsize=1.2, height="0.6", width="0.6"):
    """Draw a DAG representation in graphviz and return the resulting Digraph.
    Takes a list of graphNodes and initializes the nodes and edges.
    Coloring is expected to be a list of RGBA strings specifying how to color
    each node. Expected to be same length as nodes.

    @param nodes: list of graphNode objects
    @param coloring: list of RGBA strings to specify the color of each node.
           Expected to be the same length as nodes
    @param draw_label: bool, whether or not to draw graph labels
    @param fontname: string, name of font for graph labels
    @param fontsize: string, size of font for graph labels
    @param penscale: float, size of pen
    @param arrowsize: float, size of arrows
    @param height: string, height of nodes
    @param width: string, width of nodes
    @return dot: Digraph object to be rendered
    """

    if Digraph is None:
        raise Exception(
            "graphviz not available, will not be able to draw graph")
    else:
        # create a dot object for the graph
        dot = Digraph(format="eps", engine="dot")
        dot.attr(ratio="1.5")
        dot.attr(rotate="0")

        for ni, node in enumerate(nodes):
            if coloring is not None:
                color = coloring[ni]
            else:
                color = "#ffffff"

            if draw_label:
                if node.get_label()=='1':
                    print(node.get_label())
                    dot.node(str(node), label=node.get_label(), style="filled",
                         fillcolor=color, fontname=fontname, fontsize=fontsize,fontcolor="white",
                         height=height, width=width)
                else:
                    dot.node(str(node), label=node.get_label(), style="filled",
                         fillcolor=color, fontname=fontname, fontsize=fontsize,
                         height=height, width=width)
            else:
                dot.node(str(node), label=' ', style="filled",
                         fillcolor=color, fontname=fontname, fontsize=fontsize,
                         height=height, width=width)

        for ni, node in enumerate(nodes):
            edges = node.get_edges()
            for edge in edges:
                dot.edge(str(node),
                         str(edge),
                         arrowsize=str(arrowsize),
                         color="black",
                         penwidth=str(penscale))

    return dot


def draw_dag(dag_fn, nodes, paths, path_prob, keys,
             # 2nd set of parameters are for rendering the heatmap
             heatmap=True, colormap="Purples", penscale=0.6, arrowsize=1.2,
             fontname="Helvetica", fontsize="18", height="0.6", width="0.6",
             draw_label=True):
    """
    Function to render the DAG with heatmap information.
    @param dag_fn: string, filename path
    @param nodes: list of graphNode objects for which the graph will be drawn
    @param paths: list of lists containing all paths visited by the graphNode
           objects
    @param path_prob: list of probabilities for each path, (path_prob from
           score_graph())
    @param keys: states visited in the graph (list of keys to the state_dict)
    @param heatmap: Boolean to determine whether or not to write the dag with
           a heatmap based on the probability of each state (default: True)
    @param colormap: string, colormap used by the dag to represent probability.
           Chooses from those available in matplotlib
           (https://matplotlib.org/stable/users/explain/colors/colormaps.html)
           (default: "Purples").
    @param penscale: float, size of the pen used to draw arrows on the dag
    @param arrowsize: float, size of arrows connecting states on the dag
    @param fontname: string, font used for the labels on the dag
    @param fontsize: string, font size used for the labels on the dag
    @param height: string, height of each node on the dag
    @param width: string, width of each node on the dag
    @param draw_label: Boolean to determine whether or not to draw state
           labels on the dag
    """

    # determines if heatmap will be overlaid on top of DAG
    if heatmap:

        if cm is None or clr is None:
            raise Exception(
                "matplotlib not available, will not be able to draw graph")
        else:

            default_cmap = cm.get_cmap(colormap)

            # make a list of counts for each node to color
            coloring = np.zeros(len(nodes), dtype=float)
            for path, p in zip(paths, path_prob):
                for n in path:
                    coloring[int(n.get_index())] += 1 * p

            # normalize probability
            for t in keys:
                b = np.array([t == n.get_time() for n in nodes])
                coloring[b] /= coloring[b].sum()

            # convert probability to colors
            cmap_colors = [clr.to_hex(default_cmap(color))
                           for color in coloring]

            dot = draw_dag_in_graphviz(
                nodes, coloring=cmap_colors, draw_label=draw_label,
                fontname=fontname, fontsize=fontsize, penscale=penscale,
                arrowsize=arrowsize, height=height, width=width)
            dot.render(dag_fn)

    # no heatmap
    else:
        dot = draw_dag_in_graphviz(
            nodes, coloring=None, draw_label=draw_label, fontname=fontname,
            fontsize=fontsize, penscale=penscale, arrowsize=arrowsize,
            height=height, width=width)
        dot.render(dag_fn)

if __name__ == "__main__":
    # it is important that everything starts from main dir
    main_dir = os.getcwd()
    os.chdir(main_dir)
    state_dict = {'0min': 3, '1min': 3, '2min': 1}

    # then trajectory model is created based on the all copied data
    expected_subcomplexes = ['A', 'B', 'C']
    exp_comp = {'A': 'exp_compA.csv', 'B': 'exp_compB.csv', 'C': 'exp_compC.csv'}
    input = '../data/'
    output = "../output/"

    nodes, graph, graph_prob, graph_scores = spatiotemporal.create_DAG(state_dict,
                                                                       input_dir=input, scorestr='_scores.log',
                                                                       output_dir=output, spatio_temporal_rule=True,
                                                                       out_cdf=False, out_labeled_pdf=False, out_pdf=False,
                                                                       expected_subcomplexes=expected_subcomplexes,
                                                                       score_comp=True, exp_comp_map=exp_comp,
                                                                       draw_dag=False)

    draw_dag('dag_heatmap_manuscript',nodes,graph,graph_prob,state_dict.keys())










