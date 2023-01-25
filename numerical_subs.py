from Graph import *
import sympy as sp
from itertools import product

"""

The functions in this file substitute (i.e. replace) a sp.Symbol by a 
numerical value in an input x, where x is usually a sp.Matrix.

"""

def set_to_one_sigma_eps_mat(graph, x):
    """
    This method substitutes \sigma_{\epsilon_j}$ by 0 in x.

    Parameters
    ----------
    graph: Graph
    x: sp.Matrix or sp.Symbol

    Returns
    -------
    type(x)

    """
    num_nds = graph.num_nds
    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_sym = "sigma_eps_" + str(i)
        x = x.subs({sp.Symbol(sigma_eps_sym): 1})
    return x

def set_to_zero_gains_without_arrows(graph, x):
    """
    This method substitutes \sigma_{x_j}$ by 0 in x.


    Parameters
    ----------
    graph: Graph
    x: sp.Matrix or sp.Symbol

    Returns
    -------
    type(x)
    """
    num_nds = graph.num_nds
    for row, col in product(range(num_nds), range(num_nds)):
        row_nd = graph.ord_nodes[row]
        col_nd = graph.ord_nodes[col]

        if row > col and (col_nd, row_nd) not in graph.edges:
            alp_sym = "alp_" + str(row) + "_L_" + str(col)
            x = x.subs({sp.Symbol(alp_sym):0})
    return x

if __name__ == "__main__":
    from core_matrices import *

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
        dim = graph.num_nds

        x = sigma_eps_sym_mat(dim)
        print("\n", x)
        x = set_to_one_sigma_eps_mat(graph, x)
        print(x)

        print()
        # dot has 5 arrows out of maximum of 6 possible
        # for 4 nodes. Hence one of the gains is zero.
        # arrow 2->3 (i.e., "b"->"s") is missing
        # so alp_3_L_2=0
        print(graph.ord_nodes)
        print(graph.edges)

        x = alp_sym_mat(dim)
        print("\n", x)

        x = set_to_zero_gains_without_arrows(graph, x)
        print(x)

    main()
