from Graph import *
import numpy as np
from itertools import product
import pandas as pd


class RandomDataMaker:
    """
    This purpose of this class is to generate, for a linear SCM WITHOUT
    feedback loops, a synthetic dataset with: (1) column labels= the names
    of the nodes graph.ord_nodes, (2)instances of node values in each row.
    To generate this, we require 'alp_mat', 'graph' and 'sigma_eps'.


    Attributes
    ----------
    alp_mat: np.array of shape=(dim,dim), where dim = number of nodes.
        The matrix of alphas (i.e., gains \alpha_{i|j})
    graph: Graph
        information about DAG structure
    sigma_eps: list[float]
        the standard deviations of \epsilon_j, which is the external root
        node pointing into x_j. The entries in this list are ordered
        according to 'graph.ord_nodes'

    """

    def __init__(self, graph, sig_eps, alp_mat=None, alp_bound=1.0):
        """
        Constructor.

        In this constructor, an alp_mat not equal to 'None' can be
        submitted, or, if alp_mat == None, an alp_mat will be generated
        randomly. If generated randomly, each non-zero gain \alpha_{ i|j} is
        chosen from the uniform distribution over the interval [ -alp_bound,
        alp_bound]

        Parameters
        ----------
        graph: Graph
        sig_eps: list[float]
        alp_mat: np.array of shape=(dim, dim)
        alp_bound: float
            must be a positive number.
        """
        self.graph = graph
        dim = graph.num_nds
        assert len(sig_eps) == dim
        self.sigma_eps = sig_eps
        if alp_mat is None:
            self.alp_mat = self.generate_random_gains(alp_bound)
        else:
            assert alp_mat.shape == (dim, dim)
            self.alp_mat = alp_mat

    def generate_random_gains(self, alp_bound):
        """
        In this internal method, the gains \alpha_{i|j} are generated
        randomly. Each non-zero \alpha_{ i|j} is chosen from the uniform
        distribution over the interval [-alp_bound, alp_bound].

        Parameters
        ----------
        alp_bound: float
            must be a positive number.

        Returns
        -------
        np.array of shape=(dim, dim)

        """
        assert alp_bound > 0
        dim = self.graph.num_nds
        alp_mat = np.zeros((dim, dim))
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            if row > col and (col_nd, row_nd) in self.graph.arrows:
                alp_mat[row, col] = np.random.uniform(-alp_bound, alp_bound)
        return alp_mat

    def generate_one_random_instance(self):
        """
        This internal method returns an array with random values for the
        nodes 'graph.ord_nodes'.

        Returns
        -------
        np.array of shape=(dim,)

        """

        dim = self.graph.num_nds
        eps_values = [np.random.normal(loc=1.0, scale=self.sigma_eps[i])
                      for i in range(dim)]
        nd_values = np.zeros((dim,))
        for i, nd in enumerate(self.graph.ord_nodes):
            for pa_nd in self.graph.nx_graph.predecessors(nd):
                j = self.graph.node_position(pa_nd)
                nd_values[i] += self.alp_mat[i, j]*nd_values[j]
            nd_values[i] += eps_values[i]

        return nd_values

    def generate_dataset_csv(self, num_rows, path):
        """
        This method writes a file which contains a dataset in the
        comma-separated-values (csv) format. The dataset has (1) column
        labels= the names of the nodes graph.ord_nodes, (2)instances of node
        values in each row.

        Parameters
        ----------
        num_rows: int
            number of rows of the dataset
        path: str
            path to the destination of the output file

        Returns
        -------
        None

        """
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
        dmaker = RandomDataMaker(graph, sig_eps=sig_eps)
        data_path = "test_data.csv"
        num_rows = 5
        dmaker.generate_dataset_csv(num_rows, data_path)
        print(pd.read_csv(data_path))
        print("------------------------------")
        alp_mat = np.zeros((dim, dim))
        alp_mat[1, 0] = 4
        alp_mat[2, 0], alp_mat[2, 1] = 2, -3
        alp_mat[3, 0], alp_mat[3, 1] = 1, -1
        sig_eps = [0.0]*dim
        dmaker = RandomDataMaker(graph, sig_eps=sig_eps,
                                 alp_mat=alp_mat)
        data_path = "test_data.csv"
        num_rows = 5
        dmaker.generate_dataset_csv(num_rows, data_path)
        print(pd.read_csv(data_path))

    main(True)
