from GainsCalculator import *

class FBackGainsCalculator(GainsCalculator):

    def __init__(self, graph):
        GainsCalculator.__init__(self, graph)
        self.mat_K = None

    def calculate_gains_sb(self, mat_K=None, time="n_plus_one"):
        dim = self.graph.num_nds
        GainsCalculator.calculate_gains_sb(self, k_sb_mat(dim), time=time)
        B = beta_sb_mat(dim)
        C_n = cov_sb_mat(dim, time="n")
        mat_A = set_to_zero_gains_without_arrows(self.graph,
                                                 alp_sb_mat(dim))
        one_minus_A = sp.eye(dim) - mat_A
        one_minus_A_inv = one_minus_A.inv()
        self.mat_K = sp.simplify(B*C_n*B.T*one_minus_A_inv)

    def print_K_mat_entries(self, verbose=False):
        return print_matrix_sb_entries(
                                self.mat_K,
                                "K",
                                self.graph,
                                verbose=verbose)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        cal = FBackGainsCalculator(graph)
        cal.calculate_gains_sb()
        cal.print_gains(verbose=True)
        cal.print_K_mat_entries(verbose=True)

    main()