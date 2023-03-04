import pandas as pd
import numpy as np
from copy import deepcopy
from itertools import product
import sympy as sp
from Graph import *
from RandomDataMaker import *
from GainsCalculator import *


class GainsEstimator:
    """
    The goal of this class is to estimate the gains \alpha_{i|j} from an
    input file that contains a dataset. The dataset has the node names
    graph.ord_nodes as column labels, and instances of node values in each row.

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
    alpha_cum_err: float
        cumulative error, equal to the sum of the absolute values of the
        errors err_i_j in gains_list. This error is calculated only if
        there are no hidden nodes; it's set to zero otherwise.
    alpha_list: list[sp.Eq]
        A list of equations of the form '\alpha_{i|j} = float' if the graph
        has an arrow x_j->x_i, or of the form 'err_{i,j} = float' if that
        arrow is missing from the graph. 'err_i_j' is an error metric equal
        to the difference between both sides of an equation that constrains
        the covariances. Exception: if there are hidden nodes, the right
        hand sides of these equations may contain symbolic expressions
        pertaining to covariances alluding to hidden nodes.
    alpha_mat_estimate: np.array of shape=(dim, dim), where dim=number of nodes
        estimate of the alpha matrix. Contains estimates for the gains
        \alpha_{i|j}. If a particular \alpha_{i|j} estimate can't be
        converted to a float (because, for example, it depends on hidden
        variables), it is set to np.nan. \alpha_{i|j} estimates for
        non-existent arrows are set to 0.
    cov_mat: sp.Matrix
        Let cov_mat_nm be the numpy, numeric (nm) covariance matrix
        calculated from the input dataset. cov_mat is a sp.Matrix of the
        same dimension as cov_mat_nm that coincides with cov_mat_nm on those
        entries that do not have a hidden node as row or column index. Those
        entries of cov_mat that do have hidden nodes in their indices,
        are symbolic (sb).
    graph: Graph
    hidden_nds: list[str] or None
        This is a list of the nodes that are hidden.
    solve_symbolically: bool
    """

    def __init__(self, graph,
                 path,
                 solve_symbolically=False,
                 hidden_nds=None):
        """

        Parameters
        ----------
        graph: Graph
        path: str
            path to input file containing dataset
        solve_symbolically: bool
            solve_symbolically=True if linsolve() is called using a fully
            symbolic covariance matrix, and then the numeric values of the
            covariance matrix are substituted in the solution.
            solve_symbolically=False if linsolve() is called using a hybrid
            covariance matrix, partly symbolic, partly numeric.
        hidden_nds: None or list[str]
        """
        self.graph = graph
        df = None
        if path is not None:
            df = pd.read_csv(path)
            assert set(df.columns) == set(graph.ord_nodes)
            # put columns in same order as graph.ord_nodes
            df = df[graph.ord_nodes]
        self.solve_symbolically = solve_symbolically
        if hidden_nds is None:
            self.hidden_nds = []
        else:
            assert set(hidden_nds).issubset(graph.ord_nodes)
            self.hidden_nds = hidden_nds

        dim = graph.num_nds
        self.alpha_mat_estimate = np.zeros((dim, dim))
        self.alpha_cum_err = 0
        self.alpha_list = None

        self.cov_mat = None
        if df is not None:
            self.set_cov_mat(df)
            self.calculate_gains()
            self.fix_alpha_list()

    def set_cov_mat(self, df):
        """
        This method sets the value of the sp.Matrix called self.cov_mat.
        Entries of that matrix that have hidden nodes in their indices,
        are symbolic. All other entries are numeric.

        Parameters
        ----------
        df: pd.Dataframe

        Returns
        -------
        None

        """
        if df is None:
            assert False
        cov_mat_nm = df.cov().to_numpy()
        dim = self.graph.num_nds
        self.cov_mat = cov_sb_mat(dim, time=None)
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            symbolic = (row_nd in self.hidden_nds) or\
                       col_nd in self.hidden_nds
            if not symbolic:
                self.cov_mat[row, col] = cov_mat_nm[row, col]

    def calculate_gains(self):
        """
        This method sets the value of self.alpha_list.

        Parameters
        ----------

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        calc = GainsCalculator(self.graph)
        if self.solve_symbolically:
            cov_mat_in = cov_sb_mat(dim, time=None)
        else:
            cov_mat_in = self.cov_mat

        calc.calculate_gains(cov_mat_in=cov_mat_in,
                                 mat_K=None, time=None)
        self.alpha_list = calc.alpha_list

    def fix_alpha_list(self):
        """

        Returns
        -------

        """
        dim = self.graph.num_nds
        for i in range(len(self.alpha_list)):
            eq = self.alpha_list[i]
            str0 = str(eq.args[0])
            # last 2 terms in split("_")
            if str0[0:5] != "alpha":
                # print("hhgd", str0.split("_")[-2:])
                row_str, col_str = str0.split("_")[-2:]
                str1 = "err" + "_" + row_str + "_" + col_str
                eq = sp.Eq(sp.Symbol(str1),
                           sp.Symbol(eq.args[0])-sp.Symbol(eq.args[1]))
            for row, col in product(range(dim), range(dim)):
                row_nd = self.graph.ord_nodes[row]
                col_nd = self.graph.ord_nodes[col]
                symbolic = (row_nd in self.hidden_nds) or\
                          (col_nd in self.hidden_nds)
                if not symbolic:
                    sb_str = sb_cov_str(row, col, time=None)
                    eq = eq.subs(sb_str, self.cov_mat[row, col])
            self.alpha_list[i] = eq
            # print("llkkl", eq)
            str1 = str(eq.args[1])
            try:
                xx = float(str1)
            except ValueError:
                xx = np.nan
            if str0[0:5] == "alpha":
                row_str, col_str = str0[6:len(str0)].split("_L_")
                row, col = int(row_str), int(col_str)
                self.alpha_mat_estimate[row, col] = xx
            else:
                self.alpha_cum_err += abs(xx)

    def get_alpha_list_comments(self, true_alpha_mat):
        """

        Parameters
        ----------
        true_alpha_mat

        Returns
        -------

        """
        comments = []
        for i in range(len(self.alpha_list)):
            str0 = str(self.alpha_list[i].args[0])
            if str0[0:5] == "alpha":
                row_str, col_str = str0[6:].split("_L_")
                row, col = int(row_str), int(col_str)
                comments.append("(true= " +
                                ("%.6f" %true_alpha_mat[row, col]) + ")")
            else:
                comments.append("")
        return comments


    def print_alpha_list(self, true_alpha_mat=None, verbose=False):
        """
        This method renders in latex, in a jupyter notebook (but not on the
        console), the estimates of the gains \alpha_{i|j} if arrow $x_j->x_i$
        is present in the DAG, or of err_i_j if that arrow is absent in the
        DAG. It also prints, if it's available, the true \alpha_{i|j} next to
        its estimate. Iff verbose=True, it also prints the same thing in ASCII,
        in both the console and jupyter notebook.

        Parameters
        ----------
        true_alpha_mat: np.array of shape=(dim, dim)
            This input contains the true alpha matrix alpha_mat, if one is
            known. One would be known if the input dataset was generated by
            RandomDataMaker.
        verbose: bool

        Returns
        -------
        None

        """
        comments = self.get_alpha_list_comments(true_alpha_mat)
        return print_list_sb(self.alpha_list, self.graph,
                      verbose=verbose, time=None,
                      comment_list=comments, round=True)


if __name__ == "__main__":
    def main():
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        dot_path = 'tempo13.txt'
        # dot_path = 'dot_atlas/good_bad_trols_G1.dot'
        graph = Graph(dot_path)
        dim = graph.num_nds
        mean_eps = [0]*dim
        sig_eps = [10]*dim
        alpha_bound = 10
        dmaker = RandomDataMaker(graph,
                                 mean_eps=mean_eps,
                                 sig_eps=sig_eps,
                                 alpha_bound=alpha_bound)
        num_rows = 100
        data_path = "test_data.csv"
        dmaker.write_dataset_csv(num_rows, data_path)
        df = pd.read_csv(data_path)
        print(df)
        for solve_symbolically in [False, True]:
            print("************** solve_symbolically=", solve_symbolically)
            gest = GainsEstimator(graph, data_path,
                                  solve_symbolically=solve_symbolically)
            gest.print_alpha_list(true_alpha_mat=dmaker.alpha_mat, verbose=True)
            print("alpha_mat_estimate=\n", gest.alpha_mat_estimate)
            print("alpha_cum_err=", gest.alpha_cum_err)
            gest = GainsEstimator(graph, data_path,
                                  solve_symbolically=solve_symbolically,
                                  hidden_nds=["s"])
            gest.print_alpha_list(true_alpha_mat=dmaker.alpha_mat, verbose=True)
            print("alpha_mat_estimate=\n", gest.alpha_mat_estimate)
            print("alpha_cum_err=", gest.alpha_cum_err)

    main()
