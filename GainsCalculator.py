from core_matrices import *
from numerical_subs import *
from sympy.solvers.solveset import linsolve
from copy import deepcopy


class GainsCalculator:
    """
    The purpose of this class is to calculate/store the gains \alpha_{i|j}
    expressed symbolically as a function of the covariances <x_i, x_j>. Note
    that the nodes are in topological order so x_i happens after x_j if i>j.
    If the graph has N nodes, and it is fully connected, then gains \alpha_{
    i|j} will be calculated for all i>j, where i, j = 1, 2 , ..., N. If the
    graph is not fully connected, then for each missing arrow x_j->x_i with
    \alpha_{ i|j}=0, instead of the equation \alpha{i|j}=0, there will be a
    constraint among the covariances.

    Attributes
    ----------
    alp_mat: sp.Matrix
        the symbolic solutions for the gains \alpha_{i|j} as an sp.Matrix.
        alpha_{i|j}=0 if arrow x_j->x_i missing.
    alp_list: list[sp.Equality]
        list of symbols, where each symbol in the list is an equation of the
        form:

        \alpha_{i|j} = a function of the covariances

        If the gain \alpha_{i|j} of arrow x_j->x_i is zero because that
        arrow is missing, then \alpha_{i|j} on the left hand side of the
        above equation is replaced by <x_i, x_j>, so the equation becomes a
        constraint on the covariance matrix entries.

    graph: Graph

    """

    def __init__(self, graph):
        """
        Constructor

        Parameters
        ----------
        graph: Graph

        """
        self.graph = graph
        self.alp_list = None
        self.alp_mat = None

    def calculate_gains_sb(self, mat_K=None, time=None):
        """
        This method calculates and stores in 'self.alp_list', a list
        of symbolic equations. Each equation gives either the value of a
        gain \alpha_{i|j}, or a constraint on the covariances.

        Parameters
        ----------
        mat_K: sp.Matrix
            K matrix used only for linear SCM with feedback loops

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        if mat_K is None:
            mat_K = sp.zeros(dim)
        # print('hhgffd', mat_K)
        A = set_to_zero_gains_without_arrows(self.graph,
                                             alp_sb_mat(dim))
        self.alp_list = []
        self.alp_mat = sp.zeros(dim)
        for row in range(1, dim):
            # cov_mat = cov_sb_mat(dim)
            # cov_prod = cov_mat[0:row, 0:row].inv()*cov_mat[0:row, row]
            # cov_prod = sp.simplify(cov_prod)
            # A[row, 0:row] = cov_prod.T

            # sympy can't solve overdetermined system
            # of linear equations so fix it this way
            cov_mat = cov_sb_mat(dim, time=time)
            eqs_mat = cov_mat[0:row, 0:row] * \
                     A[row, 0:row].T - \
                      (cov_mat[0:row, row] - mat_K[0:row, row])
            eqs = [eqs_mat[i, 0] for i in range(row)]
            unknowns = []
            for i in range(row):
                if str(A[row, i]) == '0':
                    # we only use cov_mat[min(i,j), max(i,j)]
                    # because cov_mat[i, j] is symmetric.
                    # Since this system is overdetermined,
                    # make some of the covariances unknowns
                    unknowns.append(cov_mat[min(row, i), max(row, i)])
                else:
                    unknowns.append(A[row, i])
            # the comma does what is called sequence unpacking
            # draws out item from single item list
            sol_list, = linsolve(eqs, unknowns)
            # print(str(sol_list))
            for i in range(row):
                self.alp_list.append(sp.Eq(unknowns[i], sol_list[i]))
                left_str = str(unknowns[i])
                if left_str[0:3] == 'alp':
                    # print("kkkll", left_str)
                    row_str, col_str = left_str[4:].split("_L_")
                    self.alp_mat[int(row_str), int(col_str)] = sol_list[i]


    def print_gains(self, verbose=False):
        """
        This method renders in latex, in a jupyter notebook (but not on the
        console), an equation for the value of each gain \alpha_{i|j} of
        arrow x_j->x_i, or, if that gain is zero, a constraint for <x_i,
        x_j>. Iff verbose=True, it also prints the same thing in ASCII,
        in both the console and jupyter notebook.

        Parameters
        ----------
        verbose: Bool

        Returns
        -------
        sp.Symbol

        """
        return print_list_sb(self.alp_list, self.graph, verbose=verbose)


if __name__ == "__main__":
    def main():
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        path = 'tempo13.txt'
        # path = 'dot_atlas/good_bad_trols_G1.dot'
        graph = Graph(path)
        cal = GainsCalculator(graph)
        cal.calculate_gains_sb()
        cal.print_gains(verbose=True)
        print(cal.alp_mat)


    main()


