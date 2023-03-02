import sympy as sp
from itertools import product
from latexify import *

"""

The functions in this file return various core (i.e., fundamental for 
scumpy) matrices of the type sp.Matrix. The entries of those matrices are 
sp.Symbol instances.

The names of the entries of the matrices created by this file, are as follows:

    "sigma_eps_" + str(i)
    "sigma_nd_" + str(i)
    "alpha_" + str(row) + "_L_" + str(col)
    "beta_" + str(row) + "_L_" + str(col)
    "K_" + str(row) + "_" + str(col)
    "cov_" + str(row) + "_" + str(col)
    "cov_one_" + str(row) + "_" + str(col)
    "cov_n_" + str(row) + "_" + str(col)
    "cov_n_plus_one" + str(row) + "_" + str(col)
    "cov_n" + str(n) + "_" + str(row) + "_" + str(col)
    "ee_" + str(row) + "_" + str(col)
    "rho_" + str(row) + "_" + str(col)
    "pder_" + str(row) + "_wrt_" + str(col)

"""


def make_sb_mat(dim, mat_str, mat_type="general"):
    """
    This method returns a symbolic (sb) matrix of type sp.Matrix.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.
    mat_str: str
        the name of the matrix being returned. For example, if mat_str='cov',
        the returned matrix has entries cov[i, j].
    mat_type: str
        This flag must be one of the following: "general", "symmetric",
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
                    col.append(sp.Symbol(mat_str % (i, j)))
                else:
                    col.append(0)
            elif mat_type == "diagonal":
                if i == j:
                    col.append(sp.Symbol(mat_str % i))
                else:
                    col.append(0)
            else:
                assert False

        rows.append(col)
    return sp.Matrix(rows)


def sigma_eps_sb_mat(dim):
    """
    This method returns a diagonal matrix (of type sp.Matrix) with diagonal
    entries equal to the standard deviations sigma_eps_j=\sigma_{\epsilon_j}
    of \epsilon_j for each j.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "sigma_eps_%d",
                        mat_type="diagonal")


def sigma_nd_sb_mat(dim):
    """
    This method returns a diagonal matrix (of type sp.Matrix) with diagonal
    entries equal to the standard deviations sigma_nd_j = \sigma_{x_j} of
    node x_j for each j.


    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "sigma_nd_%d",
                        mat_type="diagonal")


def alpha_sb_mat(dim):
    """
    This method returns a matrix (of type sp.Matrix) of gains A with entries
    A_{ i, j} = alpha_i_L_j=\alpha_{i|j}

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "alpha_%d_L_%d",
                        mat_type="strict_lower_triangular")


def beta_sb_mat(dim):
    """
    This method returns a matrix (of type sp.Matrix) of feedback gains B
    with  entries B_{ i, j} = beta_i_L_j = \beta_{i|j}

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "beta_%d_L_%d",
                        mat_type="general")


def cov_sb_mat(dim, time=None):
    """
    This method returns the covariance matrix at time t, C^t (of type
    sp.Matrix) with entries C^t_{i,j}=<x^t_i, x^t_j> = cov_t_i_j. $t$ can be
    None, "one", "n", "n_plus_one" or an int

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.
    time: None or str or int


    Returns
    -------
    sp.Matrix

    """
    assert time in [None, "one", "n", "n_plus_one"] or \
           isinstance(time, int)
    if time is None:
        mat_str = "cov_%d_%d"
    elif isinstance(time, int):
        mat_str = "cov_n" + str(time) + "_%d_%d"
    else:
        mat_str = "cov_" + time + "_%d_%d"
    return make_sb_mat(dim, mat_str, mat_type="symmetric")


def cov2times_sb_mat(dim, time=None, delta=False):
    """
    This method returns 2-times covariance matrix C^{n,n+1} (of type
    sp.Matrix) with entries C^{n,n+1}_{i,j}=<x^{n}_i, x^{n+1}_j> =
    cov2times_i_j.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.
    time: "n" or int
        evaluate n at time

    Returns
    -------
    sp.Matrix

    """
    xtra_str = ""
    if delta:
        xtra_str = "d_"
    if isinstance(time, int):
        str0 = xtra_str + "cov2times_n" + str(time)
    elif time == "n":
        str0 = xtra_str + "cov2times_n"
    else:
        assert False, "time=" + str(time)
    str0 += "_%d_%d"

    return make_sb_mat(dim, str0,
                       mat_type="general")


def ee_sb_mat(dim):
    """
    This method returns the epsilon covariance matrix ee (of type sp.Matrix)
    with entries ee_i_j = <eps_i, eps_j>

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "ee_%d_%d",
                        mat_type="symmetric")


def rho_sb_mat(dim):
    """
    This method returns the correlation matrix \rho (of type sp.Matrix) with
    entries rho_i_j=\rho_{i, j}.

    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "rho_%d_%d",
                        mat_type="symmetric")


def jacobian_sb_mat(dim):
    """
    This method returns the Jacobian matrix J (of type sp.Matrix) with
    enties J_{i, j} = jacobian_i_j = partial derivative of x_i with respect
    to x_j


    Parameters
    ----------
    dim: int
        dimension of square matrix = number of nodes in graph.

    Returns
    -------
    sp.Matrix

    """
    return make_sb_mat(dim, "pder_%d_wrt_%d",
                        mat_type="general")


if __name__ == "__main__":

    def main():
        dim = 3
        print(sigma_eps_sb_mat(dim))
        print(sigma_nd_sb_mat(dim))
        print(alpha_sb_mat(dim))
        print(beta_sb_mat(dim))
        print(ee_sb_mat(dim))
        print(cov_sb_mat(dim, time=None))
        print(cov_sb_mat(dim, time="one"))
        print(cov_sb_mat(dim, time=2))
        print(cov2times_sb_mat(dim, time=None))
        print(cov2times_sb_mat(dim, time=5))
        print(ee_sb_mat(dim))
        print(rho_sb_mat(dim))
        print(jacobian_sb_mat(dim))
        print((sp.eye(dim) - alpha_sb_mat(dim)).inv())
    main()

