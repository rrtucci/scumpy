from DotTool import *
import networkx as nx

class Graph:
    def __init__(self, dot_file_path):
        self.path = dot_file_path
        nx_graph = DotTool.nx_graph_from_dot_file(dot_file_path)
        self.nx_graph = nx_graph
        self.ord_nodes = list(nx.topological_sort(nx_graph))
        self.edges = list(nx_graph.edges)
        self.num_nds = len(self.ord_nodes)

    def draw(self, jupyter):
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


