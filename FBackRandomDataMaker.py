from FBackGraph import *
from RandomDataMaker import *
import numpy as np
from itertools import product
import pandas as pd


class FBackRandomDataMaker(RandomDataMaker):
    """
    This purpose of this class is to generate, for a linear SCM WITH 
    feedback loops, a synthetic dataset with: (1) column labels= the names 
    of the nodes graph.ord_nodes followed by [n] for n=1,2,3, ...n_max, 
    (2) instances of node values in each row. To generate this, we require 
    'alp_mat', 'beta_mat', 'graph' and 'sigma_eps'.

    Attributes
    ----------
    alp_mat: np.array of shape=(dim,dim), where dim = number of nodes.
        The matrix of alphas (i.e., unitime gains \alpha_{i|j})
    beta_mat: np.array of shape=(dim,dim), where dim = number of nodes.
        The matrix of betas (i.e., feedback gains \beta_{i|j})
    n_max: int
        We consider times n=1,2,3,..., n_max
    graph: FBackGraph
        information about DAG structure
    sigma_eps: list[float]
        the standard deviations of \epsilon_j, which is the external root
        node pointing into x_j. The entries in this list are ordered
        according to 'graph.ord_nodes'

    """

    def __init__(self, n_max, graph, sig_eps, alp_mat=None,
                 beta_mat=None, alp_bound=1.0):
        """
        Constructor.

        In this constructor, an alp_mat (resp., beta_mat) not equal to
        'None' can be submitted, or, if alp_mat == None (resp., beta_mat ==
        None), an alp_mat (resp., beta_mat) will be generated randomly. If
        generated randomly, each non-zero gain \alpha_{ i|j} ( resp.,
        \beta_{i|j}) is chosen from the uniform distribution over the
        interval [ -alp_bound, alp_bound]

        Parameters
        ----------
        n_max: int
            We consider times n=1,2,3, ..., n_max
        graph: FBackGraph
        sig_eps: list[float]
        alp_mat: np.array of shape=(dim, dim)
        beta_mat: np.array of shape=(dim, dim)
        alp_bound: float
            must be a positive number.
        """
        self.n_max = n_max
        dim = graph.num_nds
        RandomDataMaker.__init__(self, graph, sig_eps,
                                 alp_mat=np.zeros((dim, dim)),
                                 alp_bound=alp_bound)
        self.alp_mat, self.beta_mat = self.generate_random_gains(alp_bound)
        if beta_mat is not None:
            assert beta_mat.shape == (dim, dim)
            self.beta_mat = beta_mat
        if alp_mat is not None:
            assert alp_mat.shape == (dim, dim)
            self.alp_mat = alp_mat

    @staticmethod
    def get_columns(n_max, graph):
        """
        This method returns the preferred list of column labels. First,
        self.graph.ordered_nds with "[1]" appended to each node name. Next,
        the same thing with "[2]" instead of "[1]". And so on.

        Parameters
        ----------
        n_max: int
        graph: FBackGraph

        Returns
        -------
        list[str]

        """
        dim = graph.ord_nodes
        columns = []
        for n in range(1, n_max + 1):
            x = [graph.ord_nodes[i] + "[" + str(n) + "]" for i in
                 range(dim)]
            columns += x
        return columns

    def generate_random_gains(self, alp_bound):
        """
        In this internal method, the unitime gains \alpha_{i|j} and the
        feedback gains \beta_{i|j} are generated randomly. Each non-zero
        \alpha_{ i|j} and \beta_{ i|j} is chosen from the uniform
        distribution over the interval [-alp_bound, alp_bound].

        Parameters
        ----------
        alp_bound: float
            must be a positive number.

        Returns
        -------
        np.array, np.array
            both arrays of shape=(dim, dim)

        """
        assert alp_bound > 0
        dim = self.graph.num_nds
        alp_mat = np.zeros((dim, dim))
        beta_mat = np.zeros((dim, dim))
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            if row > col and (col_nd, row_nd) in self.graph.unitime_arrows:
                alp_mat[row, col] = np.random.uniform(-alp_bound, alp_bound)
            if (col_nd, row_nd) in self.graph.fback_arrows:
                beta_mat[row, col] = np.random.uniform(-alp_bound, alp_bound)
        return alp_mat, beta_mat

    def generate_one_random_instance(self):
        """
        This internal method returns a dictionary mapping time n to random
        values for the nodes 'graph.ord_nodes'.

        Returns
        -------
        dict(int, np.array)
            np.array of shape=(dim,)

        """

        dim = self.graph.num_nds
        eps_values = [np.random.normal(loc=1.0, scale=self.sigma_eps[i])
                      for i in range(dim)]
        n_to_nd_values = {}
        for n in range(1, self.n_max+1):
            nd_values = np.zeros((dim,))
            for i, nd in enumerate(self.graph.ord_nodes):
                for pa_nd in self.graph.nx_graph.predecessors(nd):
                    j = self.graph.node_position(pa_nd)
                    nd_values[i] += self.alp_mat[i, j]*nd_values[j]
                nd_values[i] += eps_values[i]
                if n >= 2:
                    for j in range(dim):
                        nd_values[i] += self.beta_mat[i, j] *\
                                        n_to_nd_values[n-1][j]

            n_to_nd_values[n] = nd_values

        return n_to_nd_values

    def generate_dataset_csv(self, num_rows, path):
        """
        This method writes a file which contains a dataset in the
        comma-separated-values (csv) format. The dataset has: (1) column
        labels= the names of the nodes graph.ord_nodes followed by [n] for
        n=1,2,3, ...n_max, (2) instances of node values in each row.
        
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
        dim = self.graph.num_nds
        columns = FBackRandomDataMaker.get_columns(self.n_max, self.graph)
        df = pd.DataFrame(columns=columns)
        for row in range(num_rows):
            n_to_nd_values = self.generate_one_random_instance()
            for n, nd_position in product(range(1, self.n_max+1), range(dim)):
                col_name = self.graph.ord_nodes[nd_position] + \
                           "[" + str(n) + "]"
                df.loc[row, col_name] = n_to_nd_values[n][nd_position]
        df.to_csv(path, index=False)


if __name__ == "__main__":
    def main(draw):
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        if draw:
            graph.draw(jupyter=False)
        dim = graph.num_nds
        sig_eps = [2]*dim
        n_max = 4
        dmaker = FBackRandomDataMaker(n_max, graph, sig_eps=sig_eps)
        data_path = "fback_test_data.csv"
        num_rows = 5
        dmaker.generate_dataset_csv(num_rows, data_path)
        print(pd.read_csv(data_path))
        print(dmaker.alp_mat)
        print(dmaker.beta_mat)

    main(False)
