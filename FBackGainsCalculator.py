from GainsCalculator import *


class FBackGainsCalculator(GainsCalculator):
    """
    This class is a subclass of 'GainsCalculator'. Whereas the parent class
    is for analyzing DAGs without feedback loops, this class can handle
    feedback loops.

    This class implements formulae derived in my book "Bayesuvius", in the
    chapter entitled "LDEN diagrams with feedback loops". In that chapter,
    I consider, for a graph with feedback loops,

    the matrix A with entries A_{i,j}=\alpha_{ i|j}= inslice arrow gains

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
    
    self.alpha_list_with_betas and self.alpha_mat_with_betas are used to store 
    A(B, CM_info).
    
    self.beta_list and self.beta_mat are used to store B(CM_info).
    
    self.alpha_list and self.alpha_mat are used to store A(B(CM_info), CM_info).
    
    Attributes
    ----------
    # alpha_list and alpha_mat are inherited from parent class
    alpha_list_with_betas: list[sp.Equality]
    alpha_mat_with_betas: sp.Matrix
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

        # self.alpha_list and self.alpha_mat are inherited from parent class.
        self.alpha_list_with_betas = None
        self.alpha_mat_with_betas = None

        self.beta_list = None
        self.beta_mat = None

    def calculate_gains(self, cov_mat_list_in=None, mat_K=None, time="n"):
        """
        This method overrides the parent method. It calls the parent method
        within itself. It fills in

        self.alpha_list and self.alpha_mat

        self.alpha_list_with_betas and self.alpha_mat_with_betas

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

        if time == "n":
            time0 = "n"
            time1 = "n_plus_one"
        elif isinstance(time, int):
            time0 = time
            time1 = time + 1
        else:
            assert False

        if cov_mat_list_in is None:
            cov_mat0 = cov_sb_mat(dim, time=time0)
            cov2times = cov2times_sb_mat(dim, time=time0)
            cov_mat1 = cov_sb_mat(dim, time=time1)
        else:
            cov_mat0, cov2times, cov_mat1 = cov_mat_list_in

        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        mat_K = mat_B * cov_mat0

        calc = GainsCalculator(self.graph)
        calc.calculate_gains(cov_mat_in=cov_mat1, mat_K=mat_K, time=time1)
        self.alpha_mat_with_betas = deepcopy(calc.alpha_mat)
        self.alpha_list_with_betas = deepcopy(calc.alpha_list)

        self.alpha_mat = None
        self.alpha_list = None

        self.calculate_betas(cov_mat0=cov_mat0, cov2times=cov2times)
        self.calculate_alphas()

    def calculate_betas(self, cov_mat0=None, cov2times=None, time="n"):
        """
        This method fills in self.beta_list and self.beta_mat

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        if time == "n":
            time0 = "n"
        elif isinstance(time, int):
            time0 = time
        else:
            assert False
        if cov_mat0 is None:
            cov_mat0 = cov2times_sb_mat(dim, time=time0)
        if cov2times is None:
            cov2times = cov2times_sb_mat(dim, time=time0)

        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        eq_mat = (sp.eye(dim) - self.alpha_mat_with_betas) * \
                 cov2times.T - \
                 mat_B * cov_mat0

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
                if time == "n":
                    xtra_str = "n"
                elif isinstance(time, int):
                    xtra_str = "n" + str(time)
                else:
                    assert False
                cov2times_str = "cov2times_" + xtra_str +\
                                 "_" + str(col) + "_" + str(row)
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
        This method fills in self.alpha_list and self.alpha_mat

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        self.alpha_mat = deepcopy(self.alpha_mat_with_betas)
        self.alpha_list = deepcopy(self.alpha_list_with_betas)
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            beta_str = "beta_" + str(row) + "_L_" + str(col)
            self.alpha_mat = sp.simplify(self.alpha_mat.subs(sp.Symbol(beta_str),
                                           self.beta_mat[row, col]))
            for i in range(len(self.alpha_list)):
                self.alpha_list[i] = sp.simplify(
                    self.alpha_list[i].subs(sp.Symbol(beta_str),
                            self.beta_mat[row, col]))
        # print("ccvvf", self.alpha_mat)

    def print_alpha_list(self, verbose=False, time=None):
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

    def print_alpha_list_with_betas(self, verbose=False, time="n"):
        """
        This method prints the info in self.alpha_list_with_betas. It does
        this by calling latexify:print_list_sb().

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return print_list_sb(self.alpha_list_with_betas,
                             self.graph,
                            verbose=verbose,
                             time=time)

    def print_beta_list(self, verbose=False, time="n"):
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
                             verbose=verbose,
                             time=time)

    def print_alpha_list(self, verbose=False, time="n"):
        """
        This method prints the info in self.alpha_list. It does this by
        calling latexify:print_list_sb().

        Parameters
        ----------
        verbose: bool

        Returns
        -------
        sp.Symbol

        """
        return print_list_sb(self.alpha_list,
                             self.graph,
                             verbose=verbose,
                             time=time)


if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        cal = FBackGainsCalculator(graph)
        cal.calculate_gains()
        cal.print_alpha_list_with_betas(verbose=True)
        cal.print_beta_list(verbose=True)
        cal.print_alpha_list(verbose=True)

    main()
