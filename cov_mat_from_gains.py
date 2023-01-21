from new_matrices import *
from numerical_subs import *

def get_cov_mat_from_gains(graph):
    dim = graph.num_nds
    A = set_to_zero_gains_without_arrows(graph, alp_mat(dim))
    one_plus_A = A + sp.eye(dim)
    one_plus_A_inv = sp.simplify(one_plus_A.inv())
    sigma_eps_sq = sigma_eps_mat(dim)**2
    cov_mat = one_plus_A_inv*sigma_eps_sq*(one_plus_A_inv.T)
    sigma_nd_sq_inv = (sigma_nd_mat(dim)**2).inv()
    jacobian = cov_mat*sigma_nd_sq_inv
    return sp.simplify(cov_mat), sp.simplify(jacobian)

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
        cov_mat, jacobian = get_cov_mat_from_gains(graph)

        print(cov_mat)
        print()
        print(jacobian)

    main()


