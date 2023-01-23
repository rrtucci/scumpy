from Graph import *
import sympy as sp
from itertools import product
from core_matrices import *
from numerical_subs import *

def do_latex_subs(graph, x):
    num_nds = graph.num_nds
    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_latex_str = r"\sigma_{\underline{\epsilon}" + \
                          r"_{\underline{" + nd + r"}}}"
        sigma_eps_sym = "sigma_eps_" + str(i)
        x = x.subs(sp.Symbol(sigma_eps_sym), sp.Symbol(sigma_eps_latex_str))

    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_latex_str = r"\sigma_{\underline{\epsilon}" + \
                          r"_{\underline{" + nd + r"}}}"
        sigma_eps_sym = "sigma_eps_" + str(i)
        x = x.subs(sp.Symbol(sigma_eps_sym),
                   sp.Symbol(sigma_eps_latex_str))

        sigma_nd_latex_str = r"\sigma_{\underline{" + nd + r"}}"
        sigma_nd_sym = "sigma_" + str(i)
        x = x.subs(sp.Symbol(sigma_nd_sym),
                   sp.Symbol(sigma_nd_latex_str))

    for row, col in product(range(num_nds), range(num_nds)):
        row_nd = graph.ord_nodes[row]
        col_nd = graph.ord_nodes[col]

        if row > col:
            alp_latex_str = r"\alpha_{\underline{" + row_nd + \
                        r"}|\underline{" + col_nd + r"}}"
            alp_sym = "alp_" + str(row) + "_L_" + str(col)
            x = x.subs(sp.Symbol(alp_sym),
                       sp.Symbol(alp_latex_str))

        cov_latex_str = r"\left\langle\underline{" + row_nd + \
                    r"},\underline{" + col_nd  + r"}\right\rangle"
        cov_sym = "cov_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(cov_sym),
                   sp.Symbol(cov_latex_str))

        rho_latex_str = r"\rho_{\underline{" + row_nd + \
                    r"},\underline{" + col_nd  + r"}}"
        rho_sym = "rho_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(rho_sym),
                   sp.Symbol(rho_latex_str))

        pder_latex_str = r"\partial_{\underline{" + row_nd + \
                     r"}}\underline{" + col_nd  + r"}"
        pder_sym = "pder_" + str(row) + "_wrt_" + str(col)
        x = x.subs(sp.Symbol(pder_sym),
                   sp.Symbol(pder_latex_str))

    return x

def print_all_mats_after_latex_subs(graph):
    dim = graph.num_nds
    x = sigma_eps_sym_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = sigma_nd_sym_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = alp_sym_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = cov_sym_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = rho_sym_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

    x = pder_sym_mat(dim)
    x = do_latex_subs(graph, x)
    print("\n", x)
    print(sp.latex(x))

def get_str_for_matrix_entries(mat,
                               mat_name,
                               graph,
                               latex=False):
    dim = mat.shape[0]
    str0 = ""
    if latex:
        str0 += r"\begin{array}{l}"
    for row, col in product(range(dim), range(dim)):
        row_nd = graph.ord_nodes[row]
        col_nd = graph.ord_nodes[col]
        if mat_name == "gains" and col >= row:
            continue

        if mat_name == "cov" and latex:
            str0 += "\n" + r"\left\langle\underline{" + row_nd + "}" + \
                   r", \underline{" + col_nd + r"}\right\rangle="
        elif mat_name == "jacobian" and latex:
            str0 += "\n" + r"\frac{\partial\underline{" + row_nd + \
                   r"}}{\partial\underline{" + col_nd + r"}}="
        elif mat_name == "gains" and latex:
            str0 += "\n" + r"\alpha_{\underline{"+ row_nd +\
                r"}| \underline{" + col_nd + r"}}="
        else:
            str0 += "\n" + mat_name + "[" + str(row) + ":" + row_nd + ", " + \
                   str(col) + ":" + col_nd + "]="

        if latex:
            x = mat[row, col]
            x = do_latex_subs(graph, x)
            str0 += sp.latex(x) + "\n" + r"\\"
        else:
            str0 += str(mat[row, col]) + "\n"
    if latex:
        str0 = str0[:-2]  # remove last 2 characters (\\)
        str0 += r"\end{array}"

    return str0

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