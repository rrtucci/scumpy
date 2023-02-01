from Graph import *
import numpy as np
from itertools import product
import pandas as pd

class GraphRandomDataMaker:

    def __init__(self, graph, sig_eps=None, alp_mat=None, alp_bound=1.0):
        self.graph = graph
        dim = graph.num_nds
        if sig_eps is None:
            self.sigma_eps = [0]*dim
        else:
            assert len(sig_eps) == dim
            self.sigma_eps = sig_eps
        if alp_mat is None:
            self.alp_mat = self.generate_random_gains(alp_bound)
        else:
            assert alp_mat.shape == (dim, dim)
            self.alp_mat = alp_mat

    def generate_random_gains(self, alp_bound):
        assert alp_bound >0
        dim = self.graph.num_nds
        alp_mat = np.zeros((dim, dim))
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            if row > col and (col_nd, row_nd) in self.graph.edges:
                alp_mat[row, col] = np.random.uniform(-alp_bound, alp_bound)
        return alp_mat

    def generate_one_random_instance(self):

        dim = self.graph.num_nds
        eps_values = [np.random.normal(loc=1.0, scale=self.sigma_eps[i]) \
                      for i in range(dim)]
        nd_values = np.zeros((dim,))
        for i, nd in enumerate(self.graph.ord_nodes):
            for pa_nd in self.graph.nx_graph.predecessors(nd):
                j = self.graph.node_position(pa_nd)
                nd_values[i] += self.alp_mat[i,j]*nd_values[j]
            nd_values[i] += eps_values[i]

        return nd_values

    def generate_pandas_dataframe(self, num_rows, path):
        df = pd.DataFrame(columns=self.graph.ord_nodes)
        for row in range(num_rows):
            df.loc[row] = self.generate_one_random_instance()
        df.to_csv(path, index=False)


if __name__ == "__main__":
    def main(draw):
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        dot_path = 'tempo13.txt'
        # path = 'dot_atlas/good_bad_trols_G1.dot'
        graph = Graph(dot_path)
        if draw:
            graph.draw(jupyter=False)
        dim = graph.num_nds
        sig_eps = [2]*dim
        dmaker = GraphRandomDataMaker(graph, sig_eps=sig_eps)
        data_path = "test_data.csv"
        num_rows = 5
        dmaker.generate_pandas_dataframe(num_rows, data_path)
        print(pd.read_csv(data_path))
        print("------------------------------")
        alp_mat = np.zeros((dim, dim))
        alp_mat[1, 0] = 4
        alp_mat[2, 0], alp_mat[2, 1] = 2, -3
        alp_mat[3, 0], alp_mat[3, 1] = 1, -1
        sig_eps = [0.0]*dim
        dmaker = GraphRandomDataMaker(graph, sig_eps=sig_eps,
                                      alp_mat= alp_mat)
        data_path = "test_data.csv"
        num_rows = 5
        dmaker.generate_pandas_dataframe(num_rows, data_path)
        print(pd.read_csv(data_path))

    main(True)
