import sympy

from CovMatCalculator import *
from copy import deepcopy


class FBackCovMatCalculator(CovMatCalculator):
    """
    This class is a subclass of 'CovMatCalculator'. Whereas the parent class
    is for analyzing DAGs without feedback loops, this class can handle
    feedback loops.

    This class implements formulae derived in my book "Bayesuvius", in the
    chapter entitled "LDEN diagrams with feedback loops". In that chapter,
    I show that

    C^n = G^{n-1} C_1 (G^T)^{n-1}

    where n=1,2,3,... corresponds to time, G is called the growth matrix,
    and C_1 is the covariance matrix when n=1

    The class calculates/stores a symbolic representation of the covariance
    matrix C^t= <x^t_i, x^t_j> at time t=1, and of the growth matrix G,
    both expressed as a function of the gains \alpha_{i,j} and the feedback
    gains \beta{i|j}.

    Attributes
    ----------
    growth_mat_sb: sp.Matrix

    """

    def __init__(self, graph, conditioned_nds=None):
        """
        Constructor

        Parameters
        ----------
        graph: FBackGraph
        conditioned_nds: list[str]
        """
        CovMatCalculator.__init__(self, graph, conditioned_nds=conditioned_nds)
        self.growth_mat_sb = None

    def calculate_cov_mat_sb(self):
        """
        This method overrides CovMatCalculator.calculate_cov_mat_sb(self).
        It calls that parent method, plus, in addition, it calculates/stores
        self.one_minus_A_inv for future use.

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        CovMatCalculator.calculate_cov_mat_sb(self)
        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        self.growth_mat_sb = sp.simplify(self.one_minus_A_inv_sb*mat_B)

    def print_cov_mat_entries(self, verbose=False, time=None):
        """
        This method prevents the user from using the parent method that it
        overrides.

        Parameters
        ----------
        verbose: bool
        time: str or None

        Returns
        -------
        None

        """
        print("Do you mean 'print_initial_cov_mat_entries()'?")

    def print_jacobian_entries(self, verbose=False):
        """
        This method prevents the user from using the parent method that it
        overrides.

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        None

        """
        print("Don't use this when the graph has feedback loops.")

    def print_initial_cov_mat_entries(self, verbose=False):
        """
        This method renders in latex, in a jupyter notebook (but not in the
        console), the entries, one at a time, of the initial (time=1)
        covariance matrix stored in 'self.cov_mat_sb'. Iff verbose=True,
        it also prints the same thing in ASCII, in both the console and
        jupyter notebook.

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return CovMatCalculator.print_cov_mat_entries(self,
                                                      verbose=verbose,
                                                      time="one")

    def print_growth_mat_entries(self, verbose=False):
        """
        This method renders in latex, in a jupyter notebook (but not in the
        console), the entries, one at a time, of the growth matrix
        stored in 'self.growth_mat_sb'. Iff verbose=True, it also prints the
        same thing in ASCII, in both the console and jupyter notebook.

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return print_matrix_sb_entries(
                                self.growth_mat_sb,
                                "G",
                                self.graph,
                                verbose=verbose)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        cal = FBackCovMatCalculator(graph)
        cal.calculate_cov_mat_sb()
        cal.print_initial_cov_mat_entries(verbose=True)
        cal.print_growth_mat_entries(verbose=True)


    main()
