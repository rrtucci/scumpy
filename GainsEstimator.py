import pandas as pd
import numpy as np
from copy import deepcopy
from itertools import product
import sympy as sp
from Graph import *
from GraphRandomDataMaker import *
from GainsCalculator import *

class GainsEstimator:

    def __init__(self, graph, path, gains_calculator):
        df = pd.read_csv(path)
        self.cov_mat = df.cov().to_numpy()
        # print("ddfgh", df.cov())
        self.gains_sb_list = deepcopy(gains_calculator.gains_sb_list)
        assert self.gains_sb_list is not None
        dim = graph.num_nds
        self.alp_mat_estimate = np.zeros((dim, dim))
        self.cum_err = 0
        for i in range(len(self.gains_sb_list)):
            eq = self.gains_sb_list[i]
            str0 = str(eq.args[0])
            if str0[0:3] == "cov":
                str1 = "err" + str0[3:len(str0)]
                eq = sp.Eq(sp.Symbol(str1), eq.args[0]-eq.args[1])
            for row, col in product(range(dim), range(dim)):
                if row <= col:
                    cov_sb = sp.Symbol("cov_" + \
                         str(row) + "_" + str(col))
                    alp_sb = sp.Symbol("alp_" + \
                         str(row) + "_L_" + str(col))
                    err_sb = sp.Symbol("err_" + \
                                       str(row) + "_" + str(col))
                    eq = eq.subs({cov_sb :self.cov_mat[row, col]})
            self.gains_sb_list[i] = eq
            str1 = str(eq.args[1])
            if str0[0:3] == "alp":
                row_str, col_str = str0[4:len(str0)].split("_L_")
                row, col = int(row_str), int(col_str)
                self.alp_mat_estimate[row, col] = float(str1)
            else:
                self.cum_err += abs(float(str1))


    def print_gains(self, true_alp_mat=None):
        for eq in self.gains_sb_list:
            str0 = str(eq.args[0])
            str1 = "%.6f"%float(eq.args[1])
            str01 = str0 + "= " + str1
            if str(eq.args[0])[0:3] != "alp" or true_alp_mat is None:
                print("\n", str01 )
            else:

                # print("kkjg", str0[4:len(str0)].split("_L_"))
                row_str, col_str = str0[4:len(str0)].split("_L_")
                true = true_alp_mat[int(row_str), int(col_str)]
                true_str = "%.6f"%true
                print("\n", str01, ", true=", true_str)


if __name__ == "__main__":
    def main():
        dot = "digraph G {\n" \
              "a->b;\n" \
              "a->s;\n" \
              "n->s,a,b;\n" \
              "}"
        with open("tempo13.txt", "w") as file:
            file.write(dot)
        dot_path = 'tempo13.txt'
        # dot_path = 'dot_atlas/good_bad_trols_G1.dot'
        graph = Graph(dot_path)
        dim = graph.num_nds
        sig_eps = [0.1]*dim
        dmaker = GraphRandomDataMaker(graph, sig_eps=sig_eps)
        num_rows = 1000
        data_path = "test_data.csv"
        dmaker.generate_pandas_dataframe(num_rows, data_path)
        calc = GainsCalculator(graph)
        calc.calculate_gains_sb()
        gest = GainsEstimator(graph, data_path, calc)
        gest.print_gains(true_alp_mat=dmaker.alp_mat)
        print("alp_mat_estimate=\n", gest.alp_mat_estimate)
        print("cum_err=", gest.cum_err)

    main()