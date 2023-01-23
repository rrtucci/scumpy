import sympy as sp
from itertools import product
from latexify import *

def make_sym_mat(dim, mat_str, mat_type="full"):
    rows = []
    for i in range(dim):
        col = []
        for j in range(dim):
            if mat_type == "general":
                col.append(sp.Symbol(mat_str % (i, j)))
            elif mat_type == "symmetric":
                # we only use cov_mat[min(i,j), max(i,j)]
                # because cov_mat[i, j] is symmetric
                col.append(sp.Symbol(mat_str % (min(i, j), max(i, j))))
            elif mat_type == "strict_lower_triangular":
                if i > j:
                    col.append(sp.Symbol(mat_str% (i, j)))
                else:
                    col.append(0)
            elif mat_type == "diagonal":
                if i == j:
                    col.append(sp.Symbol(mat_str% i ))
                else:
                    col.append(0)
            else:
                assert(False)

        rows.append(col)
    return sp.Matrix(rows)

    '''
    "sigma_eps_" + str(i)
    "sigma_" + str(i)
    "alp_" + str(row) + "_L_" + str(col)
    "cov_" + str(row) + "_" + str(col)
    "rho_" + str(row) + "_" + str(col)
    "pder_" + str(row) + "_wrt_" + str(col)
    '''

def sigma_eps_sym_mat(dim):
    return make_sym_mat(dim, "sigma_eps_%d",
                        mat_type="diagonal")

def sigma_nd_sym_mat(dim):
    return make_sym_mat(dim, "sigma_%d",
                        mat_type="diagonal")

def alp_sym_mat(dim):
    return make_sym_mat(dim, "alp_%d_L_%d",
                        mat_type="strict_lower_triangular")

def cov_sym_mat(dim):
    return make_sym_mat(dim, "cov_%d_%d",
                        mat_type="symmetric")

def rho_sym_mat(dim):
    return make_sym_mat(dim, "rho_%d_%d",
                        mat_type="symmetric")

def pder_sym_mat(dim):
    return make_sym_mat(dim, "pder_%d_wrt_%d",
                        mat_type="general")


if __name__ == "__main__":

    def main():
        dim = 3
        print(sigma_eps_sym_mat(dim))
        print(sigma_nd_sym_mat(dim))
        print(alp_sym_mat(dim))
        print(cov_sym_mat(dim))
        print(rho_sym_mat(dim))
        print(pder_sym_mat(dim))
        print()
        print((sp.eye(dim) - alp_sym_mat(dim)).inv())
    main()