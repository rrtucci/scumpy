from core_matrices import *
from numerical_subs import *


class GainsCalculator:

    def __init__(self, graph):
        self.graph = graph
        self.gains_sym = None

    def calculate_gains_sym(self):
        dim = self.graph.num_nds
        A = set_to_zero_gains_without_arrows(self.graph,
                                             alp_sym_mat(dim))
        for row in range(1, dim):
            cov_mat = cov_sym_mat(dim)
            #print("dfgh", cov_mat)
            cov_prod = cov_mat[0:row, 0:row].inv()*cov_mat[0:row, row]
            # print("zxcv", cov_mat[0:row, row])
            cov_prod = sp.sympify(cov_prod)
            # print("llkj", cov_prod)
            A[row, 0:row] = cov_prod.T
            # print("eedr", A)
        self.gains_sym= A

    def print_gains(self, latex=False):
        x = self.gains_sym
        print("lllj", type(x))
        if latex:
            x = do_latex_subs(self.graph, x)
        print("xcfg", x)
        str0 = get_str_for_matrix_entries(x, "gains",
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
        cal = GainsCalculator(graph)
        cal.calculate_gains_sym()
        cal.print_gains(latex=False)


    main()


