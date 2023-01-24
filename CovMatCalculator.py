from core_matrices import *
from numerical_subs import *

class CovMatCalculator:
    
    def __init__(self, graph):
        self.graph = graph
        self.cov_mat_sym = None
        self.jacobian_sym = None

    def calculate_cov_mat_sym(self):
        dim = self.graph.num_nds
        A = set_to_zero_gains_without_arrows(self.graph,
                                             alp_sym_mat(dim))
        one_minus_A = A - sp.eye(dim)
        one_minus_A_inv = sp.simplify(one_minus_A.inv())
        sigma_eps_sq = sigma_eps_sym_mat(dim) ** 2
        cov_mat = one_minus_A_inv * sigma_eps_sq * (one_minus_A_inv.T)
        sigma_nd_sq_inv = (sigma_nd_sym_mat(dim) ** 2).inv()
        jacobian = cov_mat * sigma_nd_sq_inv
        self.cov_mat_sym = sp.simplify(cov_mat)
        self.jacobian_sym = sp.simplify(jacobian)

    def print_cov_mat_entries(self, latex=False,
                              verbose=False):
        x = self.cov_mat_sym
        if latex:
            x = do_latex_subs(self.graph, x)
        str0 = get_str_for_matrix_entries(x, "cov",
                                       self.graph, latex)
        if not latex:
            print(str0)

        if verbose:
            print("\n", sp.Symbol(str0))
        return sp.Symbol(str0)

    def print_jacobian_entries(self, latex=False,
                               verbose=False):
        x = self.jacobian_sym
        if latex:
            x = do_latex_subs(self.graph, x)
        str0 = get_str_for_matrix_entries(x, "jacobian",
                                       self.graph, latex)
        if not latex:
            print(str0)
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
        cal.print_cov_mat_entries(latex=False)
        cal.print_jacobian_entries(latex=False, verbose=True)

    main()


