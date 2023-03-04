from Graph import *
import sympy as sp
from itertools import product
from core_matrices import *
from numerical_subs import *
from copy import deepcopy

"""

The file contains a very complete substitution function do_latex_subs() that 
returns any input x after performing symbolic substitutions of some of the 
symbols in x. x is usually an sp.Matrix.

This file also contains various functions for printing matrices (sp.Matrix) 
and lists (list[sp.Matrix]).

"""


def round_expr(expr, num_digits):
    """
    This function rounds the numerical parts of any SymPy expression.

    Parameters
    ----------
    expr: most Sympy expressions
    num_digits: int

    Returns
    -------
    type(expr)

    """
    return expr.xreplace(
        {n: round(n, num_digits) for n in expr.atoms(sp.Number)})


def latex_time_superscript(time):
    """
    This method returns a superscript for various time inputs

    Parameters
    ----------
    time: None or str or int
        Must belong to list [None, "one", "n", "n_plus_one"] or be an int

    Returns
    -------
    str

    """
    if time is None:
        superscript = ""
    elif time == "one":
        superscript = r"^{[1]}"
    elif time == "n":
        superscript = r"^{[n]}"
    elif time == "n_plus_one":
        superscript = r"^{[n+1]}"
    elif isinstance(time, int):
        superscript = r"^{[" + str(time) + r"]}"
    else:
        assert False
    return superscript


def sb_cov_str(row, col, time):
    """
    This method returns a symbolic string for the entry of the covariance
    matrix at time "time" and position (row, col).

    Parameters
    ----------
    row: int
    col: int
    time: None or str or int

    Returns
    -------
    str

    """
    if time is None:
        sb_str = "cov"
    elif time in ["one", "n", "n_plus_one"]:
        sb_str = "cov_" + time
    elif isinstance(time, int):
        sb_str = "cov_n" + str(time)
    else:
        assert False
    sb_str += "_" + str(row) + "_" + str(col)
    return sb_str


def latex_cov_str(row_nd, col_nd, time):
    """
    This method returns a latex string for the entry of the covariance
    matrix at time "time" and position  (row_nd, col_nd). Here row_nd and
    col_nd are node names.


    Parameters
    ----------
    row_nd: str
    col_nd: str
    time: None or str or int

    Returns
    -------
    str

    """
    superscript = latex_time_superscript(time)
    if row_nd == col_nd:
        latex_str = r"\sigma^2_{\underline{" + row_nd + \
                        r"}" + superscript + r"}"
    else:
        latex_str = r"\left\langle\underline{" + row_nd + \
                r"}" + superscript + r",\underline{" + col_nd + \
                r"}" + superscript + r"\right\rangle"
    return latex_str


def sb_cov2times_str(row, col, time, delta=False):
    """
    This method returns a symbolic string for the entry of the following
    matrices between times time  and time+1 and at position (row, col).

    if delta=False, 2-times covariance matrix cov2times_n

    if delta=True, 2-times covariance matrix d_cov2times_n

    Parameters
    ----------
    row: int
    col: int
    time: None or str or int
    delta: bool

    Returns
    -------
    str

    """
    xtra_str = ""
    if delta:
        xtra_str = "d_"
    if time == "n":
        sb_str = xtra_str + "cov2times_n"
    elif isinstance(time, int):
        sb_str = xtra_str + "cov2times_n" + str(time)
    else:
        assert False
    sb_str += "_" + str(row) + "_" + str(col)
    return sb_str


def latex_cov2times_str(row_nd, col_nd, time, delta=False):
    """
    This method returns a latex string for the entry of the following
    matrices between times time and time+1 and at position (row_nd,
    col_nd), where row_nd and col_nd are node names.

    if delta=False, 2-times covariance matrix cov2times_n

    if delta=True, 2-times covariance matrix d_cov2times_n


    Parameters
    ----------
    row_nd: str
    col_nd: str
    time: None or str or int
    delta: bool

    Returns
    -------
    str

    """
    if time == "n":
        superscript = latex_time_superscript("n")
        superscript_plus = latex_time_superscript("n_plus_one")
    elif isinstance(time, int):
        superscript = latex_time_superscript(time)
        superscript_plus = latex_time_superscript(time + 1)
    else:
        assert False
    xtra_str = ""
    if delta:
        xtra_str = r"\Delta "
    latex_str = r"\left\langle\underline{" + row_nd + \
            r"}" + superscript + r"," + xtra_str + \
             r"\underline{" + col_nd + \
            r"}" + superscript_plus + r"\right\rangle"
    return latex_str


def do_latex_subs(graph, x, time=None):
    """
    This method substitutes

    sp.Symbol("sigma_eps_" + str(i))
    sp.Symbol("sigma_" + str(i))
    sp.Symbol("alpha_" + str(row) + "_L_" + str(col))
    sp.Symbol("beta_" + str(row) + "_L_" + str(col))
    sp.Symbol("cov_" + str(row) + "_" + str(col))
    sp.Symbol("cov_one_" + str(row) + "_" + str(col))
    sp.Symbol("cov_n_" + str(row) + "_" + str(col))
    sp.Symbol("cov_n_plus_one_" + str(row) + "_" + str(col))
    sp.Symbol("cov_n" + str(time) + "_" + str(row) + "_" + str(col))
    sp.Symbol("cov2times_n" + str(row) + "_" + str(col))
    sp.Symbol("cov2times_n" + str(time) + str(row) + "_" + str(col))
    sp.Symbol("d_cov2times_n" + str(row) + "_" + str(col))
    sp.Symbol("d_cov2times_n" + str(time) + str(row) + "_" + str(col))
    sp.Symbol("ee_" + str(row) + "_" + str(col))
    sp.Symbol("rho_" + str(row) + "_" + str(col))
    sp.Symbol("pder_" + str(row) + "_wrt_" + str(col))
    sp.Symbol("G_" + str(row) + "_" + str(col))
    sp.Symbol("err_" + str(row) + "_" + str(col))

    by their latex counterparts. str(i), str(row) and str(col) are all
    replaced by a node name from the list graph.ord_nodes

    Parameters
    ----------
    graph: Graph
    x: sp.Symbol or sp.Matrix
    time: None or str or int

    Returns
    -------
    type(x)

    """
    num_nds = graph.num_nds
    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        latex_str = r"\sigma_{\underline{\epsilon}" + \
                          r"_{\underline{" + nd + r"}}}"
        sb_str = "sigma_eps_" + str(i)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        latex_str = r"\sigma_{\underline{" + nd + r"}}"
        sb_str = "sigma_nd_" + str(i)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

    for row, col in product(range(num_nds), range(num_nds)):
        row_nd = graph.ord_nodes[row]
        col_nd = graph.ord_nodes[col]

        latex_str = r"\alpha_{\underline{" + row_nd + \
                    r"}|\underline{" + col_nd + r"}}"
        sb_str = "alpha_" + str(row) + "_L_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        latex_str = r"\beta_{\underline{" + row_nd + \
                    r"}|\underline{" + col_nd + r"}}"
        sb_str = "beta_" + str(row) + "_L_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        for time0 in [None, "one", "n", "n_plus_one"]:
            latex_str = latex_cov_str(row_nd, col_nd, time0)
            sb_str = sb_cov_str(row, col, time0)
            x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        for delta in [True, False]:
            latex_str = latex_cov2times_str(row_nd, col_nd, time="n",
                                            delta=delta)
            sb_str = sb_cov2times_str(row, col, time="n",
                                      delta=delta)
            x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        # print("vvbgh", time)
        if isinstance(time, int):
            latex_str = latex_cov_str(row_nd, col_nd, time)
            sb_str = sb_cov_str(row, col, time)
            # print("vvbg", sb_str, latex_str)
            x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

            latex_str = latex_cov2times_str(row_nd, col_nd, time)
            sb_str = sb_cov2times_str(row, col, time)
            x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        if row_nd == col_nd:
            latex_str = r"\sigma^2_{\underline{\epsilon}" + \
              r"_{\underline{" + row_nd + r"}}}"
        else:
            latex_str = r"\left\langle\underline{\epsilon}_\underline{" + \
                    row_nd + r"},\underline{\epsilon}_\underline{" + \
                    col_nd + r"}\right\rangle"
        sb_str = "ee_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        latex_str = r"\rho_{\underline{" + row_nd + \
                    r"},\underline{" + col_nd + r"}}"
        sb_str = "rho_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        latex_str = r"\frac{\partial\underline{" + row_nd + \
                     r"}}{\partial\underline{" + col_nd + r"}}"
        sb_str = "pder_" + str(row) + "_wrt_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        latex_str = r"G_{\underline{" + row_nd + \
                    r"},\underline{" + col_nd + r"}}"
        sb_str = "G_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

        latex_str = r"err_{\underline{" + row_nd + \
                    r"},\underline{" + col_nd + r"}}"
        sb_str = "err_" + str(row) + "_" + str(col)
        x = x.subs(sp.Symbol(sb_str), sp.Symbol(latex_str))

    return x


def print_all_core_mats_after_latex_subs(graph):
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

    def sb_mat_print(x, time=None):
        x = do_latex_subs(graph, x, time)
        print("\n", x)
        print(sp.latex(x))

    sb_mat_print(sigma_eps_sb_mat(dim))
    sb_mat_print(sigma_nd_sb_mat(dim))
    sb_mat_print(alpha_sb_mat(dim))
    sb_mat_print(beta_sb_mat(dim))
    sb_mat_print(cov_sb_mat(dim))
    sb_mat_print(cov_sb_mat(dim, time="one"))
    sb_mat_print(cov_sb_mat(dim, time=5))
    sb_mat_print(cov2times_sb_mat(dim))
    sb_mat_print(cov2times_sb_mat(dim, time=5))
    sb_mat_print(cov2times_sb_mat(dim, delta=True))
    sb_mat_print(cov2times_sb_mat(dim, time=5, delta=True))
    sb_mat_print(ee_sb_mat(dim))
    sb_mat_print(rho_sb_mat(dim))
    sb_mat_print(jacobian_sb_mat(dim))


def create_eq_list_from_matrix(mat, mat_name, graph, time):
    """
    This method takes as input an sp.Matrix mat and creates a list[sp.Eq]
    from it. The latter list can then be printed using latexify.print_list_sb()

    Parameters
    ----------
    mat: sp.Matrix
    mat_name: str
    graph: Graph
    time: None or str or int

    Returns
    -------
    list[sp.Eq]

    """
    eq_list = []
    dim = graph.num_nds
    for row, col in product(range(dim), range(dim)):
        if row <= col and mat_name == "alpha":
            continue
        mat_str = mat_name
        sep = "_"
        if mat_name in ["alpha", "beta"]:
            sep = "_L_"
        if mat_name == "pder":
            sep = "_wrt_"
        if isinstance(time, int):
            mat_str += str(time)
        mat_str += "_" + str(row) + sep + str(col)
        eq_list.append(sp.Eq(sp.Symbol(mat_str), mat[row, col]))
    return eq_list


def print_matrix_sb(mat, mat_name, graph, verbose=False, time=None):
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
    time: None or str or int

    Returns
    -------
    sp.Symbol

    """
    eq_list = create_eq_list_from_matrix(mat, mat_name, graph, time)
    return print_list_sb(eq_list, graph, verbose=verbose, time=time)


def print_list_sb(eq_list, graph, verbose=False,
                  time=None, comment_list=None, rounded=True):
    """
    This method renders in latex, in a jupyter notebook (but not on the
    console), a list 'eq_list' of sp.Eq. Iff verbose=True, it also
    prints the same thing in ASCII, in both the console and jupyter notebook.

    Parameters
    ----------
    eq_list: list[sp.Eq]
    graph: Graph or FBackGraph
    verbose: Bool
    time: None or str or int
    comment_list: list[str]
        This List[str] should be of the same length as eq_list
    rounded: bool
        This is True iff the numerical parts of the answer are to be rounded.

    Returns
    -------
    sp.Symbol

    """
    if comment_list is None:
        comment_list = [""]*len(eq_list)
    assert len(eq_list) == len(comment_list)
    str0 = ""
    x = eq_list
    x_copy = deepcopy(x)
    # print("lllj", type(x))
    if verbose:
        for i in range(len(x)):
            print(str(x[i]) + "\t" + comment_list[i] + "\n")
    str0 += r"\begin{array}{l}" + "\n"
    for i in range(len(x)):
        if rounded:
            x_copy[i] = round_expr(x_copy[i], 6)
        x_copy[i] = do_latex_subs(graph, x_copy[i], time)
        x_copy[i] = sp.latex(x_copy[i])
        str0 += x_copy[i] + r"\quad" + comment_list[i] + "\n" + r"\\" + "\n"
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
        print_all_core_mats_after_latex_subs(graph)

    main()
