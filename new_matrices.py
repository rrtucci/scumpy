import sympy as sp

def make_mat(dim, mat_str, mat_type="full"):
    rows = []
    for i in range(dim):
        col = []
        for j in range(dim):
            if mat_type == "full":
                col.append(sp.Symbol(mat_str % (i, j)))
            elif mat_type == "lower_triangular":
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
    "alp_" + str(i2) + "_L_" + str(i1)
    "cov_" + str(i1) + "_" + str(i2)
    "rho_" + str(i1) + "_" + str(i2)
    "pder_" + str(i2) + "_wrt_" + str(i1)
    '''

def sigma_eps_mat(dim):
    return make_mat(dim, "sigma_eps_%d",
        mat_type="diagonal")

def sigma_nd_mat(dim):
    return make_mat(dim, "sigma_%d",
        mat_type="diagonal")

def alp_mat(dim):
    return make_mat(dim, "alp_%d_L_%d",
        mat_type="lower_triangular")

def cov_mat(dim):
    return make_mat(dim, "cov_%d_%d",
        mat_type="full")

def rho_mat(dim):
    return make_mat(dim, "rho_%d_%d",
        mat_type="full")

def pder_mat(dim):
    return make_mat(dim, "pder_%d_wrt_%d",
        mat_type="full")


if __name__ == "__main__":

    def main():
        dim = 3
        print(sigma_eps_mat(dim))
        print(sigma_nd_mat(dim))
        print(alp_mat(dim))
        print(cov_mat(dim))
        print(rho_mat(dim))
        print(pder_mat(dim))
    main()