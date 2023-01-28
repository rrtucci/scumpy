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
    amputated_edges: None or list[(str,  str)]
        This is a list of edges to be amputated from edges in input dot
        file. We check that 'amputated_edges' is inside the list
        'full_edges' of edges obtained from reading the input dot file.
        'self.edges' is 'full_edges' minus 'amputated_edges'.
    edges: list[(str, str)]
        This is a list of string tuples such as ('a', 'b'), which indicates
        that an arrow points from 'a' to 'b'.
    num_nds: int
        number of nodes
    nx_graph: nx.DiGraph
            networkx graph.
    ord_nodes: list[str]
        Ordered nodes. A list of the node names in topological order. Root
        nodes first.
    path:
        path to a dot file. Such files are usually placed in the "dot_atlas'
        directory

    """
    def __init__(self, dot_file_path,
                 amputated_edges=None):
        """
        Constructor

        Parameters
        ----------
        dot_file_path: str
        amputated_edges: None or list[(str,str)]
        """
        self.path = dot_file_path
        nodes, all_edges = DotTool.read_dot_file(self.path)
        if amputated_edges is None:
            self.amputated_edges = []
            self.edges = all_edges
        else:
            self.amputated_edges = amputated_edges
            assert set(amputated_edges).issubset(set(all_edges))
            self.edges = [ed for ed in all_edges if ed not in amputated_edges]
        self.nx_graph = nx.DiGraph()
        self.nx_graph.add_edges_from(self.edges)
        self.ord_nodes = list(nx.topological_sort(self.nx_graph))
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
        if len(self.amputated_edges) == 0:
            DotTool.draw(self.path, jupyter)
            return

        new_dot = ""
        with open(self.path) as f:
            in_lines = f.readlines()
            for line in in_lines:
                if "->" not in line:
                    new_dot += line
                else:
                    split_list = line.split(sep="->")
                    # print("ffgg", split_list)
                    pa = split_list[0].strip()
                    ch_list = split_list[1].split(",")
                    ch_list = [x.strip().strip(";").strip() for x in ch_list]
                    # print("ffgg", pa)
                    # print("ffgg", ch_list)
                    for ch in ch_list:
                        ed = (pa, ch)
                        if ed in self.amputated_edges:
                            new_dot += ed[0] + " -> " + ed[1] + \
                                   " [color=red];\n"
                        else:
                            new_dot += ed[0] + " -> " + ed[1] + ";\n"
        with open("tempo1389.dot", "w") as file:
            file.write(new_dot)
        DotTool.draw('tempo1389.dot', jupyter)


if __name__ == "__main__":

    def main(draw):
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "b->s;\n" \
            "b[style = filled, color = yellow];\n"\
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        path = 'tempo13.txt'
        full_g = Graph(path)
        amp_g = Graph(path,
                      amputated_edges=[('b', 's')])
        for g in [full_g, amp_g]:
            print("++++++++++++++++++++++++++++++++++++")
            print("nodes in topological order:")
            print(g.ord_nodes)
            print("edges")
            print(g.edges)
            if draw:
                g.draw(jupyter=False)


    main(True)


