import sympy

from CovMatCalculator import *
from copy import deepcopy

class FBackCovMatCalculator(CovMatCalculator):

    def __init__(self, graph, conditioned_nds=None):
        """

        Parameters
        ----------
        graph: FBackGraph
        conditioned_nds: list[str]
        """
        CovMatCalculator.__init__(self, graph, conditioned_nds=None)
        self.mat_L = None
        self.mat_U = None
        self.growth_mat_sb = None

    def calculate_cov_mat_sb(self):
        dim = self.graph.num_nds
        CovMatCalculator.calculate_cov_mat_sb(self)
        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        self.growth_mat_sb = sp.simplify(self.one_minus_A_inv_sb*mat_B)

    def print_cov_mat_entries(self, verbose=False, time=None):
        print("Do you mean 'print_initial_cov_mat_entries()'?")

    def print_initial_cov_mat_entries(self, verbose=False):
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