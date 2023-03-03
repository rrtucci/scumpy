from core_matrices import *
from numerical_subs import *
from copy import deepcopy


class CovMatCalculator:
    """
    The purpose of this class is to calculate/store a symbolic
    representation of the covariance matrix C and the Jacobian matrix J,
    expressed as a function of the gains \alpha_{i,j}. We define C_{j,
    k} = <x_j, x_k>, where x_j are the internal nodes and \epsilon_j are the
    external ones. We define J_{i, j}= partial of x_i with respect to x_j.

    If self.conditioned_nds = None, we assume <epsilon_i, epsilon_j> =0
    for $i \neq j$.

    If list self.conditioned_nds is non-empty, we allow some non-diagonal
    <\epsilon_i, \epsilon_j> to be nonzero, because if we are conditioning
    on a collider, this can introduce a non-blocked path between two
    of the \epsilon_j.

    Attributes
    ----------
    conditioned_nds: list[str]
        List of the nodes that we want to condition on
    cov_mat_sb: sp.Matrix
        a symbol containing the covariance matrix C.
    graph: Graph
    jacobian_sb: sp.Matrix
        a symbol containing the Jacobian matrix J.
    one_minus_A_inv_sb: sp.Matrix
        (1-A).inv(), where A is the matrix of gains \alpha_{i|j}

    """

    def __init__(self, graph, conditioned_nds=None):
        """
        Constructor

        Parameters
        ----------
        graph: Graph
        conditioned_nds: None or list[str]
            Nodes that are being conditioned on (a.k.a the "controls")
        """
        self.graph = graph
        if conditioned_nds is None:
            self.conditioned_nds = []
        else:
            assert set(conditioned_nds).issubset(graph.ord_nodes)
            self.conditioned_nds = conditioned_nds

        self.cov_mat_sb = None
        self.jacobian_sb = None
        self.one_minus_A_inv_sb = None

    def calculate_cov_mat(self):
        """
        This method calculates and stores in 'self.cov_mat_sb', a symbolic
        expression for each of the entries C_{i,j} =<x_i, x_j>  of the
        covariance matrix C.

        It also calculates and stores in 'self.jacobian_sb', a symbolic
        expression for each of the entries J_{i,j} of the Jacobian matrix J.

        It also calculates and stores in 'self.one_minus_A_inv_sb',
        a symbolic expression for (1-A).inv(), where A is the strictly lower
        diagonal matrix of gains \alpha_{i|j}.


        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        mat_A = set_to_zero_gains_without_arrows(self.graph,
                                             alpha_sb_mat(dim))
        one_minus_A = sp.eye(dim) - mat_A
        self.one_minus_A_inv_sb = one_minus_A.inv()

        eps_cov = ee_sb_mat(dim)
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            if len(self.conditioned_nds) == 0:
                if row_nd != col_nd:
                    eps_cov = eps_cov.subs(eps_cov[row, col], 0)
            else:
                if (row_nd in self.conditioned_nds or
                        col_nd in self.conditioned_nds):
                    eps_cov = eps_cov.subs(eps_cov[row, col], 0)

        cov_mat = sp.simplify(self.one_minus_A_inv_sb * eps_cov *
                              self.one_minus_A_inv_sb.T)
        sigma_nd_sq_inv = sp.zeros(dim)
        for i in range(dim):
            sigma_nd_sq_inv[i, i] = 1 / cov_mat[i, i]
        jacobian = sp.simplify(cov_mat * sigma_nd_sq_inv)
        self.cov_mat_sb = cov_mat
        self.jacobian_sb = jacobian

    def print_cov_mat(self, verbose=False, time=None):
        """
        This method prints the info in self.cov_mat_sb. It does this by
        calling latexify:print_matrix_sb().


        Parameters
        ----------
        verbose: bool
        time: None or str
            must be either None, "one", "n" or "n_plus_one"

        Returns
        -------
        sp.Symbol

        """
        if time is None:
            mat_str = "cov"
        elif time in ["one", "n", "n_plus_one"]:
            mat_str = "cov_" + time
        elif isinstance(time, int):
            mat_str = "cov_n" + str(time)
        else:
            assert False
        return print_matrix_sb(self.cov_mat_sb,
                                mat_str,
                                self.graph,
                                verbose=verbose,
                                time=time)

    def print_jacobian(self, verbose=False, time=None):
        """
        This method renders in latex, in a jupyter notebook (but not in the
        console), the entries, one at a time, of the Jacobian matrix J. Iff
        verbose=True, it also prints the same thing in ASCII, in both the
        console and jupyter notebook.

        Parameters
        ----------
        verbose: bool
        time: None or str or int

        Returns
        -------
        sp.Symbol

        """
        return print_matrix_sb(self.jacobian_sb,
                                "pder",
                                self.graph,
                                verbose=verbose,
                               time=time)


if __name__ == "__main__":
    def main():
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        conditioned = True
        if not conditioned:
            path = 'tempo13.txt'
            conditioned_nds = None
        else:
            path = 'dot_atlas/good_bad_trols_G1.dot'
            conditioned_nds = ["Z"]
        graph = Graph(path)
        cal = CovMatCalculator(graph,
                               conditioned_nds=conditioned_nds)
        cal.calculate_cov_mat()
        cal.print_cov_mat(verbose=True)
        cal.print_jacobian(verbose=True)


    main()


