from Graph import *
import sympy as sp
from itertools import product

def set_to_one_sigma_eps_mat(graph, x):
    num_nds = graph.num_nds
    for i in range(num_nds):
        nd = graph.ord_nodes[i]

        sigma_eps_sym = "sigma_eps_" + str(i)
        x = x.subs({sp.Symbol(sigma_eps_sym): 1})
    return x

def set_to_zero_gains_without_arrows(graph, x):
    num_nds = graph.num_nds
    for i1, i2 in product(range(num_nds), range(num_nds)):
        nd1 = graph.ord_nodes[i1]
        nd2 = graph.ord_nodes[i2]

        if i2 > i1 and (nd1, nd2) not in graph.edges:
            alp_sym = "alp_" + str(i2) + "_L_" + str(i1)
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
