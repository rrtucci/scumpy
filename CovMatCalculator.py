from core_matrices import *
from numerical_subs import *
from copy import deepcopy


class CovMatCalculator:
    """
    The purpose of this class is to calculate/store the covariance matrix C
    and the Jacobian matrix J, expressed symbolically as a function of the
    gains \alpha_{i,j}. We define C_{j,k} = <x_j, x_k>, where x_j are the
    internal nodes and \epsilon_j are the external ones. We define J_{i,
    j}= partial of x_i with respect to x_j

    Attributes
    ----------
    cov_mat_sym: sp.Matrix
        a symbol containing the covariance matrix C.
    graph: Graph
    jacobian_sym: sp.Matrix
        a symbol containing the Jacobian matrix J.
    
    """
    
    def __init__(self, graph):
        """
        Constructor
        
        Parameters
        ----------
        graph: Graph
        """
        self.graph = graph
        self.cov_mat_sym = None
        self.jacobian_sym = None

    def calculate_cov_mat_sym(self):
        """
        This method calculates and stores in 'self.cov_mat_sym', a symbolic
        expression for each of the entries C_{i,j} =<x_i, x_j>  of the
        covariance matrix C. It also calculates and stores in
        'self.jacobian_sym', a symbolic expression for each of the entries
        J_{i,j} of the Jacobian matrix J.
        
        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        A = set_to_zero_gains_without_arrows(self.graph,
                                             alp_sym_mat(dim))
        one_minus_A = A - sp.eye(dim)
        one_minus_A_inv = sp.simplify(one_minus_A.inv())
        sigma_eps_sq = sigma_eps_sym_mat(dim) ** 2
        cov_mat = one_minus_A_inv * sigma_eps_sq * one_minus_A_inv.T
        sigma_nd_sq_inv = (sigma_nd_sym_mat(dim) ** 2).inv()
        jacobian = cov_mat * sigma_nd_sq_inv
        self.cov_mat_sym = sp.simplify(cov_mat)
        self.jacobian_sym = sp.simplify(jacobian)

    def print_cov_mat_entries(self, verbose=False):
        """
        This method renders in latex, in a jupyter notebook (but not in the
        console), the entries, one at a time, of the covariance matrix
        stored in 'self.cov_mat_sym'. Iff verbose=True, it also prints the
        same thing in ASCII, in both the console and jupyter notebook.
        
        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        x = self.cov_mat_sym
        x_copy = deepcopy(x)
        if verbose:
            print(get_str_for_matrix_entries(x_copy, "cov",
                                       self.graph, latex=False))

        x_copy = do_latex_subs(self.graph, x_copy)
        str0 = get_str_for_matrix_entries(x_copy, "cov",
                                       self.graph, latex=True)
        if verbose:
            print(str0)
        # this return prints nothing on the console, but, if
        # inserted as the last line of a jupyter cell, it renders
        # the latex in str0
        return sp.Symbol(str0)

    def print_jacobian_entries(self, verbose=False):
        """
        This method renders in latex, in a jupyter notebook (but not in the
        console), the entries, one at a time, of the Jacobian matrix J. Iff
        verbose=True, it also prints the same thing in ASCII, in both the
        console and jupyter notebook.
        
        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        x = self.jacobian_sym
        x_copy = deepcopy(x)
        if verbose:
            print(get_str_for_matrix_entries(x_copy, "jacobian",
                                             self.graph, latex=False))

        x_copy = do_latex_subs(self.graph, x_copy)
        str0 = get_str_for_matrix_entries(x_copy, "jacobian",
                                          self.graph, latex=True)
        if verbose:
            print(str0)
        # this return prints nothing on the console, but, if
        # inserted as the last line of a jupyter cell, it renders
        # the latex in str0
        return sp.Symbol(str0)


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
        cal = CovMatCalculator(graph)
        cal.calculate_cov_mat_sym()
        cal.print_cov_mat_entries(verbose=True)
        cal.print_jacobian_entries(verbose=True)

    main()


