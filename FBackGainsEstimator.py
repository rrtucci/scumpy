import pandas as pd
import numpy as np
from copy import deepcopy
from itertools import product
import sympy as sp
from FBackGraph import *
from FBackRandomDataMaker import *
from GainsEstimator import *
from FBackGainsCalculator import *


class FBackGainsEstimator(GainsEstimator):
    """
    The goal of this class is to estimate the inslice gain \alpha_{i|j} or
    the feedback gain \beta_{i|j} for each arrow x_j->x_i in a linear SCM
    with feedback loops. The estimation algorithm requires as input a file
    which contains a dataset with, for each time n=1,2, \dots, n_{max},
    the node names (plus string [ n]) as column labels, and with node
    values, at time n, as rows.

    The input dataset column labels must include ALL node names, and nothing
    else, but these column labels need not be in topological order (as they
    are in self.ord_nodes).

    A list of hidden nodes is an argument of the class constructor with None
    as default value. Columns of the input dataset corresponding to hidden
    nodes will be ignored. Hence, these column entries can be any number.
    Correlations <x_i, x_j> where x_i or x_j is a hidden node will be
    expressed symbolically (sb); otherwise, they will be expressed
    numerically (nm).

    Attributes
    ----------
    beta_cum_err: float
        Same as for alpha. See alpha explanation in parent class GainsEstimator
    beta_list: list[sp.Eq]
        Same as for alpha. See alpha explanation in parent class GainsEstimator
    beta_mat: sp.Matrix
        Same as for alpha. See alpha explanation in parent class GainsEstimator
    beta_mat_estimate: np.array
        Same as for alpha. See alpha explanation in parent class GainsEstimator
    cov_mat_list: list[sp.Matrix, sp.Matrix, sp.Matrix]
        [cov_mat0, cov2times, cov_mat1] where cov_mat0=covariance matrix at
        time n, cov2times=the 2-times covariance matrix between times n and
        n+1, and cov_mat1=covariance matrix at time n+1. This is an internal
        variable.
    delta: bool
        see explanation in docstring for class FBackGainsCalculator
    time: None or str or int

    """

    def __init__(self,
                 time,
                 graph,
                 df,
                 solve_symbolically=False,
                 hidden_nds=None,
                 delta=True):
        """
        Constructor

        Parameters
        ----------
        time: None or str or int
        graph: FBackGraph
        df: pd.Dataframe
        solve_symbolically: bool
            solve_symbolically=True if linsolve() is called using a fully
            symbolic covariance matrix, and then the numeric values of the
            covariance matrix are substituted in the solution.
            solve_symbolically=False if linsolve() is called using a hybrid
            covariance matrix, partly symbolic, partly numeric.
        hidden_nds: None or list[str]
        delta: bool
        """
        GainsEstimator.__init__(self, graph, path=None,
                       solve_symbolically=solve_symbolically,
                       hidden_nds=hidden_nds)
        self.time = time
        self.delta = delta
        dim = graph.num_nds
        # alpha version of the following already defined
        # by parent method
        self.beta_mat = None
        self.beta_mat_estimate = np.zeros((dim, dim))
        self.beta_cum_err = 0
        self.beta_list = None

        self.cov_mat_list = None
        self.set_cov_mat(df)
        self.calculate_gains()
        self.fix_alpha_list()
        self.fix_beta_list()

    def set_cov_mat(self, df):
        """
        This method sets the values of the 3 sp.Matrices in
        self.cov_mat_list = [cov_mat0, cov2times, cov_mat1]. Entries of
        these matrix that have hidden nodes in their indices, are symbolic.
        All other entries are numeric.

        Parameters
        ----------
        df: pd.Dataframe

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        columns = df.columns
        assert len(columns) == 2*dim
        cov_mat_nm = df.cov().to_numpy()
        cov_mat_list_nm = [
            cov_mat_nm[np.ix_(range(dim), range(dim))],
            cov_mat_nm[np.ix_(range(dim), range(dim, 2 * dim))],
            cov_mat_nm[np.ix_(range(dim, 2 * dim), range(dim, 2 * dim))]]

        cov_mat0 = cov_sb_mat(dim, time=self.time)
        cov2times = cov2times_sb_mat(dim, time=self.time)
        cov_mat1 = cov_sb_mat(dim, time=self.time + 1)
        self.cov_mat_list = [cov_mat0, cov2times, cov_mat1]
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            symbolic = (row_nd in self.hidden_nds) or \
                       col_nd in self.hidden_nds
            if not symbolic:
                for i in range(3):
                    self.cov_mat_list[i][row, col] = \
                        cov_mat_list_nm[i][row, col]
                    
    def calculate_gains(self):
        """
        This method creates an instance of FBackGainsCalculator and asks it
        to fill self.alpha_list and self.beta_list.

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        calc = FBackGainsCalculator(self.graph, delta=self.delta)
        if self.solve_symbolically:
            cov_mat0 = cov_sb_mat(dim, time=self.time)
            cov2times = cov2times_sb_mat(dim, time=self.time)
            cov_mat1 = cov_sb_mat(dim, time=self.time + 1)
            cov_mat_list_in = [cov_mat0, cov2times, cov_mat1]
        else:
            cov_mat_list_in = self.cov_mat_list
        calc.calculate_gains(cov_mat_list_in=cov_mat_list_in,
                             mat_K=None,
                             time=self.time)

        self.alpha_list = calc.alpha_list
        self.beta_list = calc.beta_list

    def fix_greek_list(self, name,
                       greek_list, greek_mat_estimate, greek_cum_err):
        """
        This method modifies the list "greek_list" (which must be either an
        "alpha_list" or a "beta_list"). For "greek_list": (1) it changes the
        constraint items (2) it inserts numerical values if
        solve_symbolically=True. It also calculates from "greek_list",
        the numpy matrix "greek_mat_estimate" and the float "greek_cum_err".

        Parameters
        ----------
        name: str
            either "alpha" or "beta"
        greek_list: list[sp.Eq]
            either "alpha_list" or "beta_list"
        greek_mat_estimate: np.array
            either "alpha_mat_estimate" or "beta_mat_estimate"
        greek_cum_err: float
            either "alpha_cum_err" or "beta_cum_err"

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        len0 = len(name)
        for i in range(len(greek_list)):
            eq = greek_list[i]
            str0 = str(eq.args[0])
            # last 2 terms in split("_")
            if str0[0:len0] != name:
                # print("hhgd", str0.split("_")[-2:])
                row_str, col_str = str0.split("_")[-2:]
                str1 = "err" + "_" + row_str + "_" + col_str
                eq = sp.Eq(sp.Symbol(str1),
                           eq.args[0] - eq.args[1])
            for row, col in product(range(dim), range(dim)):
                row_nd = self.graph.ord_nodes[row]
                col_nd = self.graph.ord_nodes[col]
                symbolic = (row_nd in self.hidden_nds) or\
                          (col_nd in self.hidden_nds)
                if not symbolic:
                    sb_str = sb_cov_str(row, col, time=self.time)
                    eq = eq.subs(sb_str, self.cov_mat_list[0][row, col])

                    sb_str = sb_cov2times_str(row, col, time=self.time)
                    eq = eq.subs(sb_str, self.cov_mat_list[1][row, col])

                    sb_str = sb_cov_str(row, col, time=self.time+1)
                    eq = eq.subs(sb_str, self.cov_mat_list[2][row, col])

            greek_list[i] = eq

            str1 = str(eq.args[1])
            try:
                xx = float(str1)
            except ValueError:
                xx = np.nan
            if str0[0:len0] == name:
                row_str, col_str = str0[len0+1:len(str0)].split("_L_")
                row, col = int(row_str), int(col_str)
                greek_mat_estimate[row, col] = xx
            else:
                greek_cum_err += abs(xx)

    @staticmethod
    def get_greek_list_comments(name,
                                greek_list,
                                true_greek_mat):
        """
        This method returns a list[str] of the same length as "greek_list".
        The returned list will be used as comments, to be printed to the
        right of each entry, when the entries of "greek_list" are printed.


        Parameters
        ----------
        name: str
            either "alpha" or "beta"
        greek_list: list[sp.Eq]
            either "alpha_list" or "beta_list"
        true_greek_mat: np.array
            the alpha (or beta) matrix used to calculate the synthetic data.

        Returns
        -------
        list[str]

        """
        comments = []
        len0 = len(name)
        for i in range(len(greek_list)):
            str0 = str(greek_list[i].args[0])
            if str0[0:len0] == name:
                row_str, col_str = str0[len0 + 1:].split("_L_")
                row, col = int(row_str), int(col_str)
                comments.append("(true= " +
                                ("%.6f" % true_greek_mat[row, col]) + ")")
            else:
                comments.append("")
        return comments

    def fix_alpha_list(self):
        """
        This method modifies self.alpha_list by calling self.fix_greek_list()

        Returns
        -------
        None

        """
        self.fix_greek_list("alpha",
                       self.alpha_list,
                       self.alpha_mat_estimate,
                       self.alpha_cum_err)

    def fix_beta_list(self):
        """
       This method modifies self.beta_list by calling self.fix_greek_list()


        Returns
        -------
        None

        """
        self.fix_greek_list("beta",
                       self.beta_list, 
                       self.beta_mat_estimate, 
                       self.beta_cum_err)

    def get_alpha_list_comments(self, true_alpha_mat):
        """
        This method returns a list of comments about self.alpha_list. To do
        this, it calls self.get_greek_list_comments().

        Parameters
        ----------
        true_alpha_mat: np.array

        Returns
        -------
        list[str]

        """
        return FBackGainsEstimator.get_greek_list_comments("alpha",
                                                           self.alpha_list,
                                                           true_alpha_mat)

    def get_beta_list_comments(self, true_beta_mat):
        """
        This method returns a list of comments about self.beta_list. To do
        this, it calls self.get_greek_list_comments().

        Parameters
        ----------
        true_beta_mat: np.array

        Returns
        -------
        list[str]

        """
        return FBackGainsEstimator.get_greek_list_comments("beta",
                                                           self.beta_list,
                                                           true_beta_mat)

    def print_alpha_list(self, true_alpha_mat=None, verbose=False):
        """
        This method prints the info in self.alpha_list. It does this by
        calling latexify.print_list_sb()

        Parameters
        ----------
        true_alpha_mat: np.array
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        comments = self.get_alpha_list_comments(true_alpha_mat)
        return print_list_sb(self.alpha_list, self.graph,
                             verbose=verbose, time=self.time,
                             comment_list=comments, rounded=True)
    
    def print_beta_list(self, true_beta_mat=None, verbose=False):
        """
        This method prints the info in self.beta_list. It does this by
        calling latexify.print_list_sb()

        Parameters
        ----------
        true_beta_mat: np.array
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        comments = self.get_beta_list_comments(true_beta_mat)
        return print_list_sb(self.beta_list, self.graph,
                             verbose=verbose, time=self.time,
                             comment_list=comments, rounded=True)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        dim = graph.num_nds
        mean_eps = [0]*dim
        sig_eps = [.001] * dim
        n_max = 2
        alpha_bound = 10
        beta_bound = 1
        dmaker = FBackRandomDataMaker(n_max, graph,
                                      mean_eps=mean_eps,
                                      sig_eps=sig_eps,
                                      alpha_bound=alpha_bound,
                                      beta_bound=beta_bound)
        num_rows = 100
        data_path = "test_data.csv"
        dmaker.write_dataset_csv(num_rows, data_path)
        df = pd.read_csv(data_path)
        for solve_symbolically in [False, True]:
            print("************** solve_symbolically=", solve_symbolically)
            time = 1
            gest = FBackGainsEstimator(time, graph, df,
                                  solve_symbolically=solve_symbolically)
            gest.print_alpha_list(true_alpha_mat=dmaker.alpha_mat,
                                  verbose=True)
            print("alpha_mat_estimate=\n", gest.alpha_mat_estimate)
            print("alpha_cum_err=", gest.alpha_cum_err)

            gest.print_beta_list(true_beta_mat=dmaker.beta_mat,
                                  verbose=True)
            print("beta_mat_estimate=\n", gest.beta_mat_estimate)
            print("beta_cum_err=", gest.beta_cum_err)

    main()
