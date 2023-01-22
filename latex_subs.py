from Graph import *
import sympy as sp
from itertools import product
from new_matrices import *

def do_latex_subs(graph, x):
    num_nds = graph.num_nds
    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_latex = r"\sigma_{\underline{\epsilon}" + \
                          r"_{\underline{" + nd + r"}}}"
        sigma_eps_sym = "sigma_eps_" + str(i)
        x = x.subs(sp.Symbol(sigma_eps_sym), sp.Symbol(sigma_eps_latex))

    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_latex = r"\sigma_{\underline{\epsilon}" + \
                          r"_{\underline{" + nd + r"}}}"
        sigma_eps_sym = "sigma_eps_" + str(i)
        x = x.subs(sp.Symbol(sigma_eps_sym),
                   sp.Symbol(sigma_eps_latex))

        sigma_nd_latex = r"\sigma_{\underline{" + nd + r"}}"
        sigma_nd_sym = "sigma_" + str(i)
        x = x.subs(sp.Symbol(sigma_nd_sym),
                   sp.Symbol(sigma_nd_latex))

    for i1, i2 in product(range(num_nds), range(num_nds)):
        nd1 = graph.ord_nodes[i1]
        nd2 = graph.ord_nodes[i2]

        if i2 > i1:
            alp_latex = r"\alpha_{\underline{" + nd2 + \
                        r"}|\underline{" + nd1 + r"}}"
            alp_sym = "alp_" + str(i2) + "_L_" + str(i1)
            x = x.subs(sp.Symbol(alp_sym),
                       sp.Symbol(alp_latex))

        cov_latex = r"< \underline{" + nd1 + \
                    r"},\underline{" + nd2  + "}>"
        cov_sym = "cov_" + str(i1) + "_" + str(i2)
        x = x.subs(sp.Symbol(cov_sym),
                   sp.Symbol(cov_latex))

        rho_latex = r"\rho_{\underline{" + nd1 + \
                    r"},\underline{" + nd2  + r"}}"
        rho_sym = "rho_" + str(i1) + "_" + str(i2)
        x = x.subs(sp.Symbol(rho_sym),
                   sp.Symbol(rho_latex))

        pder_latex = r"\partial_{\underline{" + nd2 + \
                     r"}}\underline{" + nd1  + r"}"
        pder_sym = "pder_" + str(i1) + "_wrt_" + str(i2)
        x = x.subs(sp.Symbol(pder_sym),
                   sp.Symbol(pder_latex))

    return x

def print_all_mats_after_latex_subs(graph):
    dim = graph.num_nds
    x = sigma_eps_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = sigma_nd_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = alp_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = cov_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = rho_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = pder_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))


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
        print_all_mats_after_latex_subs(graph)

    main()