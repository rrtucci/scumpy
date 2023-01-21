from Graph import *
import sympy as sp
from itertools import product


class NewSymbols:

    def __init__(self, graph):
        self.graph = graph
        self.define_symbols()

    def define_symbols(self):
        '''
        "sigma_eps_" + str(i)
        "sigma_" + str(i)
        "alp_" + str(i2) + "_L_" + str(i1)
        "cov_" + str(i1) + "_" + str(i2)
        "rho_" + str(i1) + "_" + str(i2)
        "pder_" + str(i2) + "_wrt_" + str(i1)
        '''

        num_nds = len(self.graph.ord_nodes)
        for i in range(num_nds):
            nd = self.graph.ord_nodes[i]

            sigma_eps_latex = r"\sigma_{\underline{\epsilon}" + \
                r"_{\underline{" + nd + r"}}}"
            sigma_eps_sym = "sigma_eps_" + str(i)
            globals()[sigma_eps_sym] = sp.Symbol(sigma_eps_latex)

            sigma_nd_latex = r"\sigma_{\underline{" + nd + r"}}"
            sigma_nd_sym = "sigma_" + str(i)
            globals()[sigma_nd_sym] = sp.Symbol(sigma_nd_latex)

        for i1, i2 in product(range(num_nds), range(num_nds)):
            nd1 = self.graph.ord_nodes[i1]
            nd2 = self.graph.ord_nodes[i2]

            alp_latex = r"\alpha_{\underline{" + nd2 + \
                        r"}\;|\;\underline{" + nd1 + r"}}"
            alp_sym = "alp_" + str(i2) + "_L_" + str(i1)
            globals()[alp_sym] = sp.Symbol(alp_latex)
            # if i2 > i1:
            #     if (nd1, nd2) not in self.graph.edges:  # (nd1, nd2) is
            #         # nd1->nd2
            #         # gain=0 for missing arrows
            #         eval(alp_sym + "= 0")
            # else:
            #     # alpha matrix is strictly lower diagonal
            #     exec(alp_sym + "= 0")

            cov_latex = r"< \underline{" + nd1 + \
                        r"},\underline{" + nd2  + "}>"
            cov_sym = "cov_" + str(i1) + "_" + str(i2)
            globals()[cov_sym] = sp.Symbol(cov_latex)

            rho_latex = r"\rho_{\underline{" + nd1 + \
                        r"},\underline{" + nd2  + r"}}"
            rho_sym = "rho_" + str(i1) + "_" + str(i2)
            globals()[rho_sym] = sp.Symbol(rho_latex)

            pder_latex = r"\partial_{\underline{" + nd2 + \
                        r"}}\underline{" + nd1  + r"}"
            pder_sym = "pder_" + str(i1) + "_wrt_" + str(i2)
            globals()[pder_sym] = sp.Symbol(pder_latex)

    def print_symbols(self):
        def print_latex(str0):
            exec("print(" + str0 + ")")
            # print(sp.sympify(str0)) # doesn't work

        num_nds = len(self.graph.ord_nodes)

        print("\nsigma_for_eps_for_i:")
        for i in range(num_nds):
            str0 = "sigma_eps_" + str(i)
            print(str0 + "=")
            print_latex(str0)

        print("\nsigma_for_i:")
        for i in range(num_nds):
            str0 = "sigma_" + str(i)
            print(str0 + "=")
            print_latex(str0)

        print("\nalpha_i2_L_i1:")
        for i1, i2 in product(range(num_nds), range(num_nds)):
            str0 = "alp_" + str(i2) + "_L_" + str(i1)
            print(str0 + "=")
            print_latex(str0)

        print("\ncov_i2_i1:")
        for i1, i2 in product(range(num_nds), range(num_nds)):
            str0 = "cov_" + str(i2) + "_" + str(i1)
            print(str0 + "=")
            print_latex(str0)

        print("\nrho_i2_i1:")
        for i1, i2 in product(range(num_nds), range(num_nds)):
            str0 = "rho_" + str(i2) + "_" + str(i1)
            print(str0 + "=")
            print_latex(str0)

        print("\npder_i2_wrt_i1:")
        for i1, i2 in product(range(num_nds), range(num_nds)):
            str0 = "pder_" + str(i2) + "_wrt_" + str(i1)
            print(str0 + "=")
            print_latex(str0)


if __name__ == "__main__":

    def main():
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "b->s\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        path = 'tempo13.txt'
        graph = Graph(path)
        nsyms = NewSymbols(graph)
        nsyms.print_symbols()

    main()


