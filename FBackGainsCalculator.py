from GainsCalculator import *

class FBackGainsCalculator(GainsCalculator):
    """
    This class is a subclass of 'GainsCalculator'. Whereas the parent class
    is for analyzing DAGs without feedback loops, this class can handle
    feedback loops.

    This class implements formulae derived in my book "Bayesuvius", in the
    chapter entitled "LDEN diagrams with feedback loops". In that chapter,
    I show that the arrow gains \alp_{i|j} for a graph with feedback loops
    can be calculated in terms of

     V1= the covariance matrix entries C^t = <x^t_i, x^t_j> at times t=n and
     t=n+1.

     V2= the feedback gains \beta_{i|j}

     A separate formula can be used to calculate the feedback gains in terms
     of V1. The final result is that we can calculate the arrow gains \alp_{
     i|j} in terms of V1.

    So, the goal of this class is to calculate a symbolic representation of
    the arrow gains \alp_{i|j} for a graph with feedback loops, expressed in
    terms of V1.
    """

    def __init__(self, graph):
        GainsCalculator.__init__(self, graph)

        self.alp_list_with_betas = None
        self.alp_mat_with_betas = None

        self.beta_list = None
        self.beta_mat = None


    def calculate_gains_sb(self, mat_K=None, time=None):
        dim = self.graph.num_nds
        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        mat_K = mat_B*cov2time_sb_mat(dim)

        GainsCalculator.calculate_gains_sb(self,
                                           mat_K,
                                           time="n_plus_one")
        self.alp_mat_with_betas = deepcopy(self.alp_mat)
        self.alp_mat = None

        self.alp_list_with_betas = deepcopy(self.alp_list)
        self.alp_list = None

        self.calculate_betas()
        self.calculate_alphas()

    def calculate_betas(self):
        dim = self.graph.num_nds
        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        eq_mat = (sp.eye(dim) - self.alp_mat_with_betas) * \
                 cov2time_sb_mat(dim).T - mat_B * cov_sb_mat(dim, time="n")

        unknowns = []
        eq_list = []
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            eq_list.append(eq_mat[row, col])
            if (col_nd, row_nd) in self.graph.fback_arrows:
                beta_str = "beta_" + str(row) + "_L_" + str(col)
                unknowns.append(sp.Symbol(beta_str))
            else:
                cov2time_str = "cov2time_" + str(row) + "_" + str(col)
                unknowns.append(sp.Symbol(cov2time_str))
        # the comma does what is called sequence unpacking.
        # draws out item from single item list
        sol_list, = linsolve(eq_list, unknowns)
        # print(str(sol_list))

        self.beta_list = []
        self.beta_mat = sp.zeros(dim)
        for i in range(len(sol_list)):
            self.beta_list.append(sp.Eq(unknowns[i], sol_list[i]))
            left_str = str(unknowns[i])
            if left_str[0:3] == 'beta':
                # print("kkkll", left_str)
                row_str, col_str = left_str[4:].split("_L_")
                self.beta_mat[int(row_str), int(col_str)] = sol_list[i]

    def calculate_alphas(self):
        dim = self.graph.num_nds
        self.alp_mat = deepcopy(self.alp_mat_with_betas)
        self.alp_list = deepcopy(self.alp_list_with_betas)
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            beta_str = "beta_" + str(row) + "_L_" + str(col)
            self.alp_mat =sp.simplify(self.alp_mat.subs(sp.Symbol(beta_str),
                                           self.beta_mat[row, col]))
            for i in range(len(self.alp_list)):
                self.alp_list[i] = sp.simplify(
                    self.alp_list[i].subs(sp.Symbol(beta_str),
                            self.beta_mat[row, col]))
        # print("ccvvf", self.alp_mat)

    def print_gains(self, verbose=False):
        assert False

    def print_alp_list_with_betas(self, verbose=False):
        return print_list_sb(self.alp_list_with_betas,
                             self.graph,
                            verbose=verbose)

    def print_beta_list(self, verbose=False):
        return print_list_sb(self.beta_list,
                             self.graph,
                             verbose=verbose)

    def print_alp_list(self, verbose=False):
        return print_list_sb(self.alp_list,
                             self.graph,
                             verbose=verbose)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        cal = FBackGainsCalculator(graph)
        cal.calculate_gains_sb()
        cal.print_alp_list_with_betas(verbose=True)
        cal.print_beta_list(verbose=True)
        cal.print_alp_list(verbose=True)

    main()