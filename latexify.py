from Graph import *
import sympy as sp
from itertools import product
from core_matrices import *
from numerical_subs import *
from copy import deepcopy

"""

The functions in this file return their input x after performing some 
symbolic substitutions of some of the symbols in x. x is usually a sp.Matrix.

"""


def do_latex_subs(graph, x):
    """
    This method substitutes

    sp.Symbol("sigma_eps_" + str(i))
    sp.Symbol("sigma_" + str(i))
    sp.Symbol("alp_" + str(row) + "_L_" + str(col))
    sp.Symbol("beta_" + str(row) + "_L_" + str(col))
    sp.Symbol("cov_two_time_" + str(row) + "_" + str(col))
    sp.Symbol("cov_" + str(row) + "_" + str(col))
    sp.Symbol("cov_one_" + str(row) + "_" + str(col))
    sp.Symbol("cov_n_" + str(row) + "_" + str(col))
    sp.Symbol("cov_n_plus_one_" + str(row) + "_" + str(col))
    sp.Symbol("eps_" + str(row) + "_" + str(col))
    sp.Symbol("rho_" + str(row) + "_" + str(col))
    sp.Symbol("pder_" + str(row) + "_wrt_" + str(col))

    by their latex counterparts, str(i), str(row) and str(col) are all
    replaced by a node name from the list graph.ord_nodes

    Parameters
    ----------
    graph: Graph
    x: sp.Symbol or sp.Matrix

    Returns
    -------
    type(x)

    """
    num_nds = graph.num_nds
    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_latex_str = r"\sigma_{\underline{\epsilon}" + \
                          r"_{\underline{" + nd + r"}}}"
        sigma_eps_sb = "sigma_eps_" + str(i)
        x = x.subs(sp.Symbol(sigma_eps_sb),
                   sp.Symbol(sigma_eps_latex_str))

        sigma_nd_latex_str = r"\sigma_{\underline{" + nd + r"}}"
        sigma_nd_sb = "sigma_nd_" + str(i)
        x = x.subs(sp.Symbol(sigma_nd_sb),
                   sp.Symbol(sigma_nd_latex_str))

    for row, col in product(range(num_nds), range(num_nds)):
        row_nd = graph.ord_nodes[row]
        col_nd = graph.ord_nodes[col]

        if row > col:
            alp_latex_str = r"\alpha_{\underline{" + row_nd + \
                        r"}|\underline{" + col_nd + r"}}"
            alp_sb = "alp_" + str(row) + "_L_" + str(col)
            x = x.subs(sp.Symbol(alp_sb),
                       sp.Symbol(alp_latex_str))

        beta_latex_str = r"\beta_{\underline{" + row_nd + \
                    r"}|\underline{" + col_nd + r"}}"
        beta_sb = "beta_" + str(row) + "_L_" + str(col)
        x = x.subs(sp.Symbol(beta_sb),
                   sp.Symbol(beta_latex_str))

        cov2_latex_str = r"\left\langle\underline{" + row_nd + \
                    r"}^{[n]},\underline{" + col_nd + r"}^{[n+1]}\right\rangle"
        cov2_sb = "cov2times_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(cov2_sb),
                   sp.Symbol(cov2_latex_str))

        for time in [None, "one", "n", "n_plus_one"]:
            superscript = time_superscript(time)
            if row_nd == col_nd:
                cov_latex_str = r"\sigma^2_{\underline{" + row_nd + \
                                r"}" + superscript + r"}"
            else:
                cov_latex_str = r"\left\langle\underline{" + row_nd + \
                        r"}" + superscript + r",\underline{" + col_nd + \
                        r"}" + superscript + r"\right\rangle"
            if time is None:
                cov_sb = "cov_"
            else:
                cov_sb = "cov_" + time + "_"
            cov_sb += str(row) + "_" + str(col)
            x = x.subs(sp.Symbol(cov_sb),
                       sp.Symbol(cov_latex_str))

        if row_nd == col_nd:
            ee_latex_str = r"\sigma^2_{\underline{\epsilon}" + \
              r"_{\underline{" + row_nd + r"}}}"
        else:
            ee_latex_str = r"\left\langle\underline{\epsilon}_\underline{" + \
                            row_nd + r"},\underline{\epsilon}_\underline{" + \
                            col_nd + r"}\right\rangle"
        ee_sb = "ee_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(ee_sb),
                   sp.Symbol(ee_latex_str))

        rho_latex_str = r"\rho_{\underline{" + row_nd + \
                    r"},\underline{" + col_nd + r"}}"
        rho_sb = "rho_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(rho_sb),
                   sp.Symbol(rho_latex_str))

        pder_latex_str = r"\partial_{\underline{" + row_nd + \
                     r"}}\underline{" + col_nd + r"}"
        pder_sb = "pder_" + str(row) + "_wrt_" + str(col)
        x = x.subs(sp.Symbol(pder_sb),
                   sp.Symbol(pder_latex_str))

    return x


def print_all_mats_after_latex_subs(graph):
    """
    This method is for debugging 'do_latex_subs()'. It creates the core
    matrices built by the functions in file core_matrices.py. Then it passes
    those core matrices through 'do_latex_subs()'. Finally, it prints the
    changed core matrices.

    Parameters
    ----------
    graph: Graph

    Returns
    -------
    None

    """
    dim = graph.num_nds

    def sb_mat_print(x):
        x = do_latex_subs(graph, x)
        print("\n", x)
        print(sp.latex(x))

    sb_mat_print(sigma_eps_sb_mat(dim))
    sb_mat_print(sigma_nd_sb_mat(dim))
    sb_mat_print(alp_sb_mat(dim))
    sb_mat_print(beta_sb_mat(dim))
    sb_mat_print(cov2times_sb_mat(dim))
    sb_mat_print(cov_sb_mat(dim))
    sb_mat_print(cov_sb_mat(dim, time="one"))
    sb_mat_print(ee_sb_mat(dim))
    sb_mat_print(rho_sb_mat(dim))
    sb_mat_print(jacobian_sb_mat(dim))


def time_superscript(time):
    if time is None:
        superscript = ""
    elif time == "one":
        superscript = r"^{[1]}"
    elif time == "n":
        superscript = r"^{[n]}"
    elif time == "n_plus_one":
        superscript = r"^{[n+1]}"
    else:
        assert False
    return superscript


def get_str_for_matrix_entries(mat,
                               mat_name,
                               graph,
                               latex=False):
    """
    If latex=True, this method returns a string which is a well-formed latex
    expression for a latex array with one column. For each i, j, there is
    one column entry of the form

    mat[i, j]= EXP

    where mat is the input sp.Matrix and EXP is a string coming from a
    symbolic expression.  Both mat[i,j] and EXP have been fully latexified
    and the integers i, j have been replaced by node names.

    If latex=False, this method outputs a line mat[i, j]= EXP for each i,
    j. The lines are not fully latexified and don't have the integers i,
    j replaced by node names.

    'mat_name' can be either "cov", "jacobian", "gains", or anything else.
    The left hand side mat[i,j] of the equation in each line depends on what
    we submit for 'mat_name'

    Parameters
    ----------
    mat: sp.Matrix
    mat_name: str
    graph: Graph
    latex: bool

    Returns
    -------
    str

    """
    dim = mat.shape[0]
    str0 = ""
    if latex:
        str0 += r"\begin{array}{l}"
    for row, col in product(range(dim), range(dim)):
        row_nd = graph.ord_nodes[row]
        col_nd = graph.ord_nodes[col]

        if mat_name == "gains" and col >= row:
            continue

        if mat_name[0:3] == "cov" and latex:
            time = mat_name[3:]
            if time != "":
                # this removes underscore at beginning
                time = time[1:]
            else:
                time = None
            superscript = time_superscript(time)
            str0 += "\n" + r"\left\langle\underline{" + row_nd + r"}" + \
                    superscript + r", \underline{" + col_nd + r"}" + \
                    superscript + r"\right\rangle="
        elif mat_name == "jacobian" and latex:
            str0 += "\n" + r"\frac{\partial\underline{" + row_nd + \
                   r"}}{\partial\underline{" + col_nd + r"}}="
        elif mat_name == "gains" and latex:
            str0 += "\n" + r"\alpha_{\underline{" + row_nd +\
                r"}| \underline{" + col_nd + r"}}="
        else:
            str0 += "\n" + mat_name + r"_{" + str(row) + r"=\underline{" + \
                    row_nd + r"}," + str(col) + r"=\underline{" + \
                    col_nd + r"}}="

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


def print_matrix_sb_entries(mat, mat_name, graph, verbose=False):
    """
    This method renders in latex, in a jupyter notebook (but not in the
    console), the entries, one at a time, of the symbolic matrix 'mat' named
    'mat_name'. Iff verbose=True, it also prints the same thing in ASCII,
    in both the console and jupyter notebook.

    Parameters
    ----------
    mat: sp.Matrix
    mat_name: str
    graph: Graph or FBackGraph
    verbose: bool

    Returns
    -------
    sp.Symbol

    """
    x = mat
    x_copy = deepcopy(x)
    if verbose:
        print(get_str_for_matrix_entries(x_copy, mat_name,
                                         graph, latex=False))

    x_copy = do_latex_subs(graph, x_copy)
    str0 = get_str_for_matrix_entries(x_copy, mat_name,
                                      graph, latex=True)
    if verbose:
        print(str0)
    # this return prints nothing on the console, but, if
    # inserted as the last line of a jupyter cell, it renders
    # the latex in str0
    return sp.Symbol(str0)


def print_list_sb(eq_list, graph, verbose=False):
    """
    This method renders in latex, in a jupyter notebook (but not on the
    console), a list 'eq_list' of sp.Equality. Iff verbose=True, it also
    prints the same thing in ASCII, in both the console and jupyter notebook.

    Parameters
    ----------
    graph: Graph or FBackGraph
    eq_list: list[sp.Equality]
    verbose: Bool

    Returns
    -------
    sp.Symbol

    """
    str0 = ""
    x = eq_list
    x_copy = deepcopy(x)
    # print("lllj", type(x))
    if verbose:
        for i in range(len(x)):
            print(str(x[i]), "\n")
    str0 += r"\begin{array}{l}" + "\n"
    for i in range(len(x)):
        x_copy[i] = sp.latex(do_latex_subs(graph, x_copy[i]))
        str0 += x_copy[i] + "\n" + r"\\" + "\n"
    str0 = str0[:-3]
    str0 += r"\end{array}"
    # print("lluj", str0)
    if verbose:
        print("\n", str0)
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
              "b->s\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        path = 'tempo13.txt'
        graph = Graph(path)
        print_all_mats_after_latex_subs(graph)

    main()
