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
    (2) node values in each row. To generate this, we require 
    'alpha_mat', 'beta_mat', 'graph' and 'sigma_eps'.

    Attributes
    ----------
    # alpha_mat is inherited from parent class
    beta_mat: np.array of shape=(dim,dim), where dim = number of nodes.
        The matrix of betas (i.e., feedback gains \beta_{i|j})
    n_max: int
        We consider times n=1,2,3,..., n_max

    """

    def __init__(self, n_max, graph, mean_eps, sig_eps, alpha_mat=None,
                 beta_mat=None, alpha_bound=1, beta_bound=1):
        """
        Constructor.

        In this constructor, an alpha_mat (resp., beta_mat) not equal to
        'None' can be submitted, or, if alpha_mat == None (resp., beta_mat ==
        None), an alpha_mat (resp., beta_mat) will be generated randomly. If
        generated randomly, each non-zero gain \alpha_{ i|j} ( resp.,
        \beta_{i|j}) is chosen from the uniform distribution over the
        interval [ -alpha_bound, alpha_bound]

        Parameters
        ----------
        n_max: int
            We consider times n=1,2,3, ..., n_max
        graph: FBackGraph
        mean_eps: list[float]
            List containing the mean value of each of the gaussian external
            noise variables. This list has length dim, where dim is the
            number of nodes.
        sig_eps: list[float]
            List containing the standard deviation of each of the gaussian
            external noise variables. This list has length dim, where dim is
            the number of nodes.
        alpha_mat: np.array of shape=(dim, dim)
        beta_mat: np.array of shape=(dim, dim)
        alpha_bound: float
            must be a positive number.
        beta_bound: float
            must be a positive number.
        """
        self.n_max = n_max
        dim = graph.num_nds
        RandomDataMaker.__init__(self, graph, mean_eps, sig_eps,
                                 alpha_mat=np.zeros((dim, dim)),
                                 alpha_bound=alpha_bound)
        self.alpha_mat, self.beta_mat = FBackRandomDataMaker.\
            generate_random_alpha_and_beta_mats(graph,
                                               alpha_bound,
                                               beta_bound)
        if beta_mat is not None:
            assert beta_mat.shape == (dim, dim)
            self.beta_mat = beta_mat
        if alpha_mat is not None:
            assert alpha_mat.shape == (dim, dim)
            self.alpha_mat = alpha_mat

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
        dim = graph.num_nds
        columns = []
        for n in range(1, n_max + 1):
            x = [graph.ord_nodes[i] + "[" + str(n) + "]" for i in
                 range(dim)]
            columns += x
        return columns

    @staticmethod
    def generate_random_alpha_and_beta_mats(graph, alpha_bound=1,
                                            beta_bound=1):
        """
        This static method generates randomly and returns the inslice gains
        \alpha_{ i|j} and the feedback gains \beta_{ i|j}. Each non-zero
        \alpha_{ i|j} and \beta_{ i|j} is chosen using my_random()

        Parameters
        ----------
        graph: FBackGraph
        alpha_bound: float
            must be a positive number.
        beta_bound: float
            must be a positive number.

        Returns
        -------
        np.array, np.array
            both arrays of shape=(dim, dim)

        """
        dim = graph.num_nds
        alpha_mat = np.zeros((dim, dim))
        beta_mat = np.zeros((dim, dim))
        for row, col in product(range(dim), range(dim)):
            row_nd = graph.ord_nodes[row]
            col_nd = graph.ord_nodes[col]
            if row > col and (col_nd, row_nd) in graph.inslice_arrows:
                alpha_mat[row, col] = my_random(alpha_bound)
            if (col_nd, row_nd) in graph.fback_arrows:
                beta_mat[row, col] = my_random(beta_bound)
        return alpha_mat, beta_mat

    def generate_one_random_instance(self):
        """
        This internal method returns a dictionary mapping time n to random
        values for the nodes 'graph.ord_nodes', for all n=1,2,3, ..., n_max.

        Returns
        -------
        dict(int, np.array)
            np.array of shape=(dim,)

        """

        dim = self.graph.num_nds
        n_to_nd_values = {}
        for n in range(1, self.n_max+1):
            nd_values = [0]*dim
            for i in range(dim):
                nd_values[i] += np.random.normal(loc=10,
                                                 scale=self.sigma_eps[i])
                if n >= 1:
                    for j in range(dim):
                        if i > j:
                            nd_values[i] += self.alpha_mat[i, j]*nd_values[j]
                if n >= 2:
                    for j in range(dim):
                        nd_values[i] += self.beta_mat[i, j] *\
                                    n_to_nd_values[n-1][j]

            n_to_nd_values[n] = nd_values

        return n_to_nd_values

    def write_dataset_csv(self, num_rows, path):
        """
        This method writes a file which contains a dataset in the
        comma-separated-values (csv) format. The dataset has: (1) column
        labels= the names of the nodes graph.ord_nodes followed by [n] for
        n=1,2,3, ...n_max, (2) node values in each row.
        
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
        mean_eps = [0]*dim
        sig_eps = [10]*dim
        n_max = 4
        alpha_bound = 10
        beta_bound = 1
        dmaker = FBackRandomDataMaker(n_max, graph,
                                      mean_eps=mean_eps,
                                      sig_eps=sig_eps,
                                      alpha_bound=alpha_bound,
                                      beta_bound=beta_bound)
        data_path = "fback_test_data.csv"
        num_rows = 100
        dmaker.write_dataset_csv(num_rows, data_path)
        print(pd.read_csv(data_path))
        print("alpha_mat=\n", dmaker.alpha_mat)
        print("beta_mat=\n", dmaker.beta_mat)

    main(False)
