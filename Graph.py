from DotTool import *
import networkx as nx

class Graph:
    """
    The purpose of this class is to store information related to a DAG.This
    information is derived from a dot file located at 'dot_file_path'.

    A dot file is a text file that describes a single (usually) DAG in the
    dot language. The dot language is the language used by the graph
    rendering language GraphViz.

    Attributes
    ----------
    edges: list[(str, str)]
        This is a list of string tuples such as ('a', 'b'), which indicates
        that an arrow points from 'a' to 'b'.
    num_nds: int
        number of nodes
    nx_graph: nx.DiGraph
            networkx graph.
    ord_nds: list[str]
        Ordered nodes. A list of the node names in topological order. Root
        nodes first.
    path:
        path to a dot file. Such files are usually placed in the "dot_atlas'
        directory

    """
    def __init__(self, dot_file_path):
        """
        Constructor

        Parameters
        ----------
        dot_file_path: str
        """
        self.path = dot_file_path
        nx_graph = DotTool.nx_graph_from_dot_file(dot_file_path)
        self.nx_graph = nx_graph
        self.ord_nodes = list(nx.topological_sort(nx_graph))
        self.edges = list(nx_graph.edges)
        self.num_nds = len(self.ord_nodes)

    def draw(self, jupyter):
        """
        This method draws the graph either on the console (jupyter=False) or
        in a jupyter notebook (jupyter=True)

        Parameters
        ----------
        jupyter: bool

        Returns
        -------
        None

        """
        DotTool.draw(self.path, jupyter)


if __name__ == "__main__":

    def main(draw):
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "b->s\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        path = 'tempo13.txt'
        g = Graph(path)
        print("nodes in topological order:")
        print(g.ord_nodes)
        print("edges")
        print(g.edges)
        if draw:
            g.draw(jupyter=False)

    main(True)


