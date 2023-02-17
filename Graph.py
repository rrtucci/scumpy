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
    amputated_arrows: None or list[(str,  str)]
        This is a list of arrows to be amputated from arrows in input dot
        file. We check that 'amputated_arrows' is inside the list
        'full_arrows' of arrows obtained from reading the input dot file.
        'self.arrows' is 'full_arrows' minus 'amputated_arrows'.
    arrows: list[(str, str)]
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
    def __init__(self,
                 dot_file_path,
                 amputated_arrows=None,
                 is_DAG=True):
        """
        Constructor

        Parameters
        ----------
        dot_file_path: str
        amputated_arrows: None or list[(str,str)]
        is_DAG: bool
            Whether the graph is a DAG (directed acyclic graph)
        """
        self.path = dot_file_path
        nodes, all_arrows = DotTool.read_dot_file(self.path)
        self.num_nds = len(nodes)

        if amputated_arrows is None:
            self.amputated_arrows = []
            self.arrows = all_arrows
        else:
            self.amputated_arrows = amputated_arrows
            assert set(amputated_arrows).issubset(set(all_arrows))
            self.arrows = [ed for ed in all_arrows if ed not in amputated_arrows]

        self.nx_graph = None
        self.ord_nodes = None
        if is_DAG:
            self.nx_graph = nx.DiGraph()
            self.nx_graph.add_arrows_from(self.arrows)
            # this bombs if not DAG
            self.ord_nodes = list(nx.topological_sort(self.nx_graph))

    @staticmethod
    def get_pa_and_ch_list(line):
        """
        This is an internal utility function that extracts the name of a
        parent and its children from the string 'line'.

        Parameters
        ----------
        line: str

        Returns
        -------
        str, list[str]

        """
        # get rid of arrow attributes
        line = line.split(sep="[")[0]
        split_list = line.split(sep="->")
        # pa=parent, ch=children
        pa = split_list[0].strip()
        ch_list = split_list[1].split(",")
        ch_list = [x.strip().strip(";").strip() for x in ch_list]
        return pa, ch_list

    def draw(self, jupyter=False):
        """
        This method draws the graph either on the console (jupyter=False) or
        in a jupyter notebook (jupyter=True). Amputated arrows are drawn in
        red, non-amputated ones in black.

        Parameters
        ----------
        jupyter: bool

        Returns
        -------
        None

        """
        new_dot = ""
        with open(self.path) as f:
            in_lines = f.readlines()
            for line in in_lines:
                if "->" not in line:
                    new_dot += line
                else:
                    pa, ch_list = Graph.get_pa_and_ch_list(line)
                    for ch in ch_list:
                        new_dot += pa + " -> " + ch
                        if (pa, ch) in self.amputated_arrows:
                            new_dot += " [color=red];\n"
                        else:
                            new_dot += ";\n"
        with open("tempo1389.dot", "w") as file:
            file.write(new_dot)
        DotTool.draw('tempo1389.dot', jupyter)

    def node_position(self, nd_name):
        """
        This method returns the position of the string 'nd_name' in the
        list 'self.ordered_nds'.

        Parameters
        ----------
        nd_name: str

        Returns
        -------
        int

        """
        for i, nd in enumerate(self.ord_nodes):
            if nd == nd_name:
                return i
        assert False, nd_name + " is not in " + str(self.ord_nodes)


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
                      amputated_arrows=[('b', 's')])
        for g in [full_g, amp_g]:
            print("++++++++++++++++++++++++++++++++++++")
            print("nodes in topological order:")
            print(g.ord_nodes)
            print("arrows")
            print(g.arrows)
            if draw:
                g.draw(jupyter=False)


    main(True)


