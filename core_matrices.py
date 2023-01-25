import sympy as sp
from itertools import product
from latexify import *

"""

The functions in this file return various core (i.e., fundamental for 
scumpy) matrices of the type sp.Matrix. The entries of those matrices are 
sp.Symbol instances.

The names of the entries of the matrices created by this file, are as follows:

    "sigma_eps_" + str(i)
    "sigma_" + str(i)
    "alp_" + str(row) + "_L_" + str(col)
    "cov_" + str(row) + "_" + str(col)
    "rho_" + str(row) + "_" + str(col)
    "pder_" + str(row) + "_wrt_" + str(col)

"""

def make_sym_mat(dim, mat_str, mat_type="full"):
    """
    This method returns a matrix of type sp.Matrix.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.
    mat_str: str
        the name of the matrix being returned. For example, if mat_str='cov',
        the returned matrix has entries cov[i, j].
    mat_type: str
        This flag must be one of these: "general", "symmetric",
        "strict_lower_triangular", "diagonal"


    Returns
    -------
    sp.Matrix

    """
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

def sigma_eps_sym_mat(dim):
    """
    This method returns a diagonal matrix (of type sp.Matrix) with the
    standard deviations \sigma_{\epsilon_j} of \epsilon_j along its diagonal.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sym_mat(dim, "sigma_eps_%d",
                        mat_type="diagonal")

def sigma_nd_sym_mat(dim):
    """
    This method returns a diagonal matrix (of type sp.Matrix) with the
    standard deviations \sigma_{ x_j} of node x_j along its diagonal.


    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sym_mat(dim, "sigma_%d",
                        mat_type="diagonal")

def alp_sym_mat(dim):
    """
    This method returns a matrix (of type sp.Matrix) of gains M with M_{ i,
    j} = \alpha_{i|j}

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.


    Returns
    -------
    sp.Matrix

    """
    return make_sym_mat(dim, "alp_%d_L_%d",
                        mat_type="strict_lower_triangular")

def cov_sym_mat(dim):
    """
    This method returns the covariance matrix C (of type sp.Matrix)with C_{
    i,j} = <x_i, x_j>

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.


    Returns
    -------
    sp.Matrix

    """
    return make_sym_mat(dim, "cov_%d_%d",
                        mat_type="symmetric")

def rho_sym_mat(dim):
    """
    This method returns the correlation matrix (of type sp.Matrix) \rho
    with \rho_{i, j}.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.


    Returns
    -------
    sp.Matrix

    """
    return make_sym_mat(dim, "rho_%d_%d",
                        mat_type="symmetric")

def jacobian_sym_mat(dim):
    """
    This method returns the Jacobian matrix J matrix (of type sp.Matrix)
    with the J_{i,j} = partial derivative of x_i with respect to x_j.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.


    Returns
    -------
    sp.Matrix

    """
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
        print(jacobian_sym_mat(dim))
        print()
        print((sp.eye(dim) - alp_sym_mat(dim)).inv())
    main()