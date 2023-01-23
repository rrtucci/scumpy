from core_matrices import *
from numerical_subs import *
from sympy.solvers.solveset import linsolve



class GainsCalculator:

    def __init__(self, graph):
        self.graph = graph
        self.gains_sym_list = None

    def calculate_gains_sym(self):
        dim = self.graph.num_nds
        A = set_to_zero_gains_without_arrows(self.graph,
                                             alp_sym_mat(dim))
        self.gains_sym_list = []
        for row in range(1, dim):
            # cov_mat = cov_sym_mat(dim)
            # #print("dfgh", cov_mat)
            # cov_prod = cov_mat[0:row, 0:row].inv()*cov_mat[0:row, row]
            # # print("zxcv", cov_mat[0:row, row])
            # cov_prod = sp.simplify(cov_prod)
            # # print("llkj", cov_prod)
            # A[row, 0:row] = cov_prod.T
            # # print("eedr", A)
            cov_mat = cov_sym_mat(dim)
            eqs_mat = cov_mat[0:row, 0:row]*\
                      (A[row, 0:row].T)-cov_mat[0:row, row]
            eqs = [eqs_mat[i, 0] for i in range(row)]
            unknowns = []
            # sympy can't solve overdetermined system
            # of linear equations so fix it this way
            for i in range(row):
                if str(A[row, i]) == '0':
                    # we only use cov_mat[min(i,j), max(i,j)]
                    # because cov_mat[i, j] is symmetric
                    # Since this system is overdetermined
                    # make some of the correlations unknowns
                    unknowns.append(cov_mat[min(row, i), max(row, i)])
                else:
                    unknowns.append(A[row, i])
            sol_list, = linsolve(eqs, unknowns)
            # print(str(sol_list))
            for i in range(row):
               self.gains_sym_list.append(
                   sp.Eq(unknowns[i], sol_list[i]))

    def print_gains(self, latex=False):
        str0 = ""
        x = self.gains_sym_list
        # print("lllj", type(x))
        if not latex:
            for i in range(len(x)):
                print(str(x[i]), "\n")
        else:
            str0 += r"\begin{array}{l}" + "\n"
            for i in range(len(x)):
                x[i] = sp.latex(do_latex_subs(self.graph, x[i]))
                str0 += x[i] + "\n" + r"\\" + "\n"
            str0 = str0[:-3]
            str0 += r"\end{array}"
            # print("lluj", str0)
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


