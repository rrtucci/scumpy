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
    alpha_list: list[sp.Equality]
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
    alpha_cum_err: float
        cumulative error, equal to the sum of the absolute values of the
        errors err_i_j in gains_list. This error is calculated only if
        there are no hidden nodes; it's set to zero otherwise.
    graph: FBackGraph
    hidden_nds: list[str] or None
        This is a list of the nodes that are hidden.
    """

    def __init__(self,
                 time,
                 graph,
                 df,
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
        GainsEstimator.__init__(self, graph, path=None,
                       solve_symbolically=solve_symbolically,
                       hidden_nds=hidden_nds)
        self.time = time
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
        dim = self.graph.num_nds
        calc = FBackGainsCalculator(self.graph)
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
        dim = self.graph.num_nds
        len0= len(name)
        for i in range(len(greek_list)):
            eq = greek_list[i]
            str0 = str(eq.args[0])
            # last 2 terms in split("_")
            if str0[0:len0] != name:
                # print("hhgd", str0.split("_")[-2:])
                row_str, col_str = str0.split("_")[-2:]
                str1 = "err" + "_" + row_str + "_" + col_str
                eq = sp.Eq(sp.Symbol(str1), eq.args[0]-eq.args[1])
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
        self.fix_greek_list("alpha",
                       self.alpha_list, 
                       self.alpha_mat_estimate, 
                       self.alpha_cum_err)

    def fix_beta_list(self):
        self.fix_greek_list("beta",
                       self.beta_list, 
                       self.beta_mat_estimate, 
                       self.beta_cum_err)

    def get_alpha_list_comments(self, true_alpha_mat):
        return FBackGainsEstimator.get_greek_list_comments("alpha",
                                                           self.alpha_list,
                                                           true_alpha_mat)
        
    def get_beta_list_comments(self, true_beta_mat):
        return FBackGainsEstimator.get_greek_list_comments("beta",
                                                           self.beta_list,
                                                           true_beta_mat)

    def print_alpha_list(self, true_alpha_mat=None, verbose=False):
        comments = self.get_alpha_list_comments(true_alpha_mat)
        return print_list_sb(self.alpha_list, self.graph,
                             verbose=verbose, time=self.time,
                             comment_list=comments, round=True)
    
    def print_beta_list(self, true_beta_mat=None, verbose=False):
        comments = self.get_beta_list_comments(true_beta_mat)
        return print_list_sb(self.beta_list, self.graph,
                             verbose=verbose, time=self.time,
                             comment_list=comments, round=True)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        dim = graph.num_nds
        sig_eps = [1] * dim
        n_max = 2
        alpha_bound = 1
        dmaker = FBackRandomDataMaker(n_max, graph, sig_eps=sig_eps,
                                      alpha_bound=alpha_bound)
        num_rows = 1000
        data_path = "test_data.csv"
        dmaker.generate_dataset_csv(num_rows, data_path)
        df = pd.read_csv(data_path)
        for solve_symbolically in [False, True]:
            print("************** solve_symbolically=", solve_symbolically)
            time=1
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
