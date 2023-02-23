from GainsCalculator import *


class FBackGainsCalculator(GainsCalculator):
    """
    This class is a subclass of 'GainsCalculator'. Whereas the parent class
    is for analyzing DAGs without feedback loops, this class can handle
    feedback loops.

    This class implements formulae derived in my book "Bayesuvius", in the
    chapter entitled "LDEN diagrams with feedback loops". In that chapter,
    I consider, for a graph with feedback loops,

    the matrix A with entries A_{i,j}=\alp_{ i|j}= unitime arrow gains

    the matrix B with entries B_{i,j}=\beta_{i|j}= feedback arrow gains

    I show that A and B satisfy a system of 2 linear equations with two
    unknowns A,B.

    Let CM_info = the single-time covariance matrices C^n, C^{n+1} with
    entries C^t_{i,i}=<x^t_i, x^t_j> at times t=n and t=n+1, and the 2-times
    covariance matrix C2times with entries C2times_{i,j}= <x^{ n}_i, x^{n+1}_j>

    To solve that system of 2 equations in (A, B), this class first solves
    one equation for A in terms of B and CMinfo, thus obtaining A(B,
    CMinfo). Then it substitutes A(B, CMinfo) into the remaining equation to
    obtain B(CMinfo). Finally, it substitutes B(CMinfo) into A(B, CMinfo) to
    get A(B(CMinfo), CMinfo).
    
    self.alp_list_with_betas and self.alp_mat_with_betas are used to store 
    A(B, CM_info).
    
    self.beta_list and self.beta_mat are used to store B(CM_info).
    
    self.alp_list and self.alp_mat are used to store A(B(CM_info), CM_info).
    
    Attributes
    ----------
    # alp_list and alp_mat are inherited from parent class
    alp_list_with_betas: list[sp.Equality]
    alp_mat_with_betas: sp.Matrix
    beta_list: list[sp.Equality]
    beta_mat: sp.Matrix

    """

    def __init__(self, graph):
        """
        Constructor

        Parameters
        ----------
        graph: FBackGraph
        """
        GainsCalculator.__init__(self, graph)

        # self.alp_list and self.alp_mat are inherited from parent class.

        self.alp_list_with_betas = None
        self.alp_mat_with_betas = None

        self.beta_list = None
        self.beta_mat = None

    def calculate_gains(self, mat_K=None, time=None):
        """
        This method overrides the parent method. It calls the parent method
        within itself. It fills in

        self.alp_list and self.alp_mat

        self.alp_list_with_betas and self.alp_mat_with_betas

        self.beta_list and self.beta_mat

        Parameters
        ----------
        mat_K: sp.Matrix
        time: None or str

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        mat_K = mat_B*cov2times_sb_mat(dim)

        GainsCalculator.calculate_gains(self,
                                           mat_K=mat_K,
                                           time="n_plus_one")
        self.alp_mat_with_betas = deepcopy(self.alp_mat)
        self.alp_list_with_betas = deepcopy(self.alp_list)

        self.alp_mat = None
        self.alp_list = None

        self.calculate_betas()
        self.calculate_alphas()

    def calculate_betas(self):
        """
        This method fills in self.beta_list and self.beta_mat

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        eq_mat = (sp.eye(dim) - self.alp_mat_with_betas) * \
                 cov2times_sb_mat(dim).T - mat_B * cov_sb_mat(dim, time="n")

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
                cov2times_str = "cov2times_" + str(col) + "_" + str(row)
                unknowns.append(sp.Symbol(cov2times_str))
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
        """
        This method fills in self.alp_list and self.alp_mat

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        self.alp_mat = deepcopy(self.alp_mat_with_betas)
        self.alp_list = deepcopy(self.alp_list_with_betas)
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            beta_str = "beta_" + str(row) + "_L_" + str(col)
            self.alp_mat = sp.simplify(self.alp_mat.subs(sp.Symbol(beta_str),
                                           self.beta_mat[row, col]))
            for i in range(len(self.alp_list)):
                self.alp_list[i] = sp.simplify(
                    self.alp_list[i].subs(sp.Symbol(beta_str),
                            self.beta_mat[row, col]))
        # print("ccvvf", self.alp_mat)

    def print_gains(self, verbose=False):
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
        assert False

    def print_alp_list_with_betas(self, verbose=False):
        """
        This method prints the info in self.alp_list_with_betas. It does
        this by calling latexify:print_list_sb().

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return print_list_sb(self.alp_list_with_betas,
                             self.graph,
                            verbose=verbose)

    def print_beta_list(self, verbose=False):
        """
        This method prints the info in self.beta_list. It does this by
        calling latexify:print_list_sb().

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return print_list_sb(self.beta_list,
                             self.graph,
                             verbose=verbose)

    def print_alp_list(self, verbose=False):
        """
        This method prints the info in self.alp_list. It does this by
        calling latexify:print_list_sb().

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return print_list_sb(self.alp_list,
                             self.graph,
                             verbose=verbose)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        cal = FBackGainsCalculator(graph)
        cal.calculate_gains()
        cal.print_alp_list_with_betas(verbose=True)
        cal.print_beta_list(verbose=True)
        cal.print_alp_list(verbose=True)

    main()
