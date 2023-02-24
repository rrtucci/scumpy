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
    numerically (nu).

    Attributes
    ----------
    alp_list: list[sp.Equality]
        A list of equations of the form '\alpha_{i|j} = float' if the graph
        has an arrow x_j->x_i, or of the form 'err_{i,j} = float' if that
        arrow is missing from the graph. 'err_i_j' is an error metric equal
        to the difference between both sides of an equation that constrains
        the covariances. Exception: if there are hidden nodes, the right
        hand sides of these equations may contain symbolic expressions
        pertaining to covariances alluding to hidden nodes.
    alp_mat_estimate: np.array of shape=(dim, dim), where dim=number of nodes
        estimate of the alpha matrix. Contains estimates for the gains
        \alpha_{i|j}. This array is only calculated if there are no hidden
        nodes; it's set to zero otherwise.
    cov_mat: sp.Matrix
        Let cov_mat_nu be the numpy, numeric (nu) covariance matrix
        calculated from the input dataset. cov_mat is a sp.Matrix of the
        same dimension as cov_mat_nu that coincides with cov_mat_nu on those
        entries that do not have a hidden node as row or column index. Those
        entries of cov_mat that do have hidden nodes in their indices,
        are symbolic (sb).
    cum_err: float
        cumulative error, equal to the sum of the absolute values of the
        errors err_i_j in gains_list. This error is calculated only if
        there are no hidden nodes; it's set to zero otherwise.
    graph: Graph
    hidden_nds: list[str] or None
        This is a list of the nodes that are hidden.
    """

    def __init__(self, graph, path, symbolic_solve=False, hidden_nds=None):
        """

        Parameters
        ----------
        graph: Graph
        path: str
            path to input file containing dataset
        symbolic_solve: bool
            symbolic_solve=True if linsolve() is called using a fully
            symbolic covariance matrix, and then the numeric values of the
            covariance matrix are substituted in the solution.
            symbolic_solve=False if linsolve() is called using a hybrid
            symbolic covariance matrix, partly symbolic, partly numeric.
        hidden_nds: None or list[str]
        """
        self.graph = graph
        df = pd.read_csv(path)
        assert set(df.columns) == set(graph.ord_nodes)
        # put columns in same order as graph.ord_nodes
        df = df[graph.ord_nodes]
        if hidden_nds is None:
            self.hidden_nds = []
        else:
            assert set(hidden_nds).issubset(graph.ord_nodes)
            self.hidden_nds = hidden_nds
        self.cov_mat = None
        self.set_cov_mat(df)
        dim = graph.num_nds
        self.alp_mat_estimate = np.zeros((dim, dim))
        self.cum_err = 0
        self.alp_list = None
        self.set_alp_list(symbolic_solve)

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
        cov_mat_nu = df.cov().to_numpy()
        dim = self.graph.num_nds
        self.cov_mat = cov_sb_mat(dim, time=None)
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            symbolic = (row_nd in self.hidden_nds) or\
                       col_nd in self.hidden_nds
            if not symbolic:
                self.cov_mat[row, col] = cov_mat_nu[row, col]

    def set_alp_list(self, symbolic_solve):
        """
        This method sets the value of self.alp_list.

        Parameters
        ----------
        symbolic_solve: bool

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        calc = GainsCalculator(self.graph)
        if not symbolic_solve:
            cov_mat_in = cov_sb_mat(dim, time=None)
        else:
            cov_mat_in = self.cov_mat

        calc.calculate_gains(cov_mat_in=cov_mat_in,
                                 mat_K=None, time=None)
        self.alp_list = calc.alp_list
        assert self.alp_list is not None
        for i in range(len(self.alp_list)):
            eq = self.alp_list[i]
            str0 = str(eq.args[0])

            if str0[0:3] == "cov":
                str1 = "err" + str0[3:len(str0)]
                eq = sp.Eq(sp.Symbol(str1), eq.args[0]-eq.args[1])
            for row, col in product(range(dim), range(dim)):
                row_nd = self.graph.ord_nodes[row]
                col_nd = self.graph.ord_nodes[col]
                symbolic = (row_nd in self.hidden_nds) or\
                          (col_nd in self.hidden_nds)
                if row <= col and not symbolic:
                    cov_sb = sp.Symbol("cov_" +
                         str(row) + "_" + str(col))
                    eq = eq.subs({cov_sb: self.cov_mat[row, col]})
            self.alp_list[i] = eq
            if len(self.hidden_nds) == 0:
                str1 = str(eq.args[1])
                if str0[0:3] == "alp":
                    row_str, col_str = str0[4:len(str0)].split("_L_")
                    row, col = int(row_str), int(col_str)
                    self.alp_mat_estimate[row, col] = float(str1)
                else:
                    self.cum_err += abs(float(str1))

    def print_gains(self, true_alp_mat=None, verbose=False,
                    switch_alp2beta=False):
        """
        This method renders in latex, in a jupyter notebook (but not on the
        console), the estimates of the gains \alp_{i|j} if arrow $x_j->x_i$
        is present in the DAG, or of err_i_j if that arrow is absent in the
        DAG. It also prints, if it's available, the true \alp_{i|j} next to
        its estimate. Iff verbose=True, it also prints the same thing in ASCII,
        in both the console and jupyter notebook.

        Parameters
        ----------
        true_alp_mat: np.array of shape=(dim, dim)
            This input contains the true alpha matrix alp_mat, if one is
            known. One would be known if the input dataset was generated by
            RandomDataMaker.
        verbose: bool
        switch_alp2beta: bool
            This replaces the letter alpha (for unitime gains) by beta (for
            feedback gains).

        Returns
        -------
        None

        """
        full_str = r"\begin{array}{l}" + "\n"
        if not switch_alp2beta:
            mat_name = "alp"
            g_letter = "alpha"
        else:
            mat_name = "beta"
            g_letter = "beta"
        for eq in self.alp_list:
            str0 = str(eq.args[0])
            alp_eq = True
            if str(eq.args[0])[0:3] == mat_name:
                row_str, col_str = str0[4:len(str0)].split("_L_")
            else:
                alp_eq = False
                row_str, col_str = str0[4:len(str0)].split("_")
            row, col = int(row_str), int(col_str)
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            if alp_eq:
                sb_str0 = "\\" + g_letter + r"_{\underline{" + row_nd + "}|"
            else:
                sb_str0 = r"err_{\underline{" + row_nd + "},"
            sb_str0 += r"\underline{" + col_nd + r"}}"
            if len(self.hidden_nds) == 0:
                str1 = "%.6f" % float(eq.args[1])
            else:
                def round_expr(expr, num_digits):
                    return expr.xreplace(
                        {n: round(n, num_digits) for n in expr.atoms(
                            sp.Number)})

                str1 = round_expr(eq.args[1], 6)
                str1 = do_latex_subs(self.graph, str1)
                str1 = sp.latex(str1)
            str01 = sb_str0 + "= " + str1
            if not alp_eq or true_alp_mat is None:
                full_str += str01 + r"\\" + "\n"
            else:
                true = true_alp_mat[row, col]
                true_str = "%.6f" % true
                full_str += str01 + r",\quad true=" + true_str + r"\\" + "\n"
        full_str = full_str[:-3] + "\n" + r"\end{array}"
        if verbose:
            print(full_str)
        return sp.Symbol(full_str)


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
        sig_eps = [0.1]*dim
        dmaker = RandomDataMaker(graph, sig_eps=sig_eps)
        num_rows = 100
        data_path = "test_data.csv"
        dmaker.generate_dataset_csv(num_rows, data_path)
        for symbolic_solve in [False, True]:
            gest = GainsEstimator(graph, data_path,
                                  symbolic_solve=symbolic_solve)
            gest.print_gains(true_alp_mat=dmaker.alp_mat, verbose=True)
            print("alp_mat_estimate=\n", gest.alp_mat_estimate)
            print("cum_err=", gest.cum_err)
            gest = GainsEstimator(graph, data_path,
                                  symbolic_solve=symbolic_solve,
                                  hidden_nds=["s"])
            gest.print_gains(true_alp_mat=dmaker.alp_mat, verbose=True)

    main()
