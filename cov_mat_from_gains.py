from new_matrices import *
from numerical_subs import *


def get_cov_mat_from_gains(graph, verbose=True, latex=False):
    dim = graph.num_nds
    A = set_to_zero_gains_without_arrows(graph, alp_mat(dim))
    one_minus_A = A - sp.eye(dim)
    one_minus_A_inv = sp.simplify(one_minus_A.inv())
    sigma_eps_sq = sigma_eps_mat(dim)**2
    cov_mat = one_minus_A_inv*sigma_eps_sq*(one_minus_A_inv.T)
    sigma_nd_sq_inv = (sigma_nd_mat(dim)**2).inv()
    jacobian = cov_mat*sigma_nd_sq_inv
    cov_mat = sp.simplify(cov_mat)
    jacobian = sp.simplify(jacobian)
    if verbose:
        print_sympy_matrix_entries(cov_mat, "cov_mat", graph, latex)
        print()
        print_sympy_matrix_entries(jacobian, "jacobian", graph, latex)

    return sp.Symbol(
        sp.latex(do_latex_subs(graph,cov_mat)) + "\n" +
            sp.latex(do_latex_subs(graph,jacobian)))

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
        graph = Graph(path)
        get_cov_mat_from_gains(graph, verbose=True, latex=True)

    main()


