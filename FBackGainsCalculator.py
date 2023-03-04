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
    covariance matrix cov2times_n with entries cov2times_n_{i,j}= <x^{ n}_i,
    x^{n+1}_j>

    To solve that system of 2 equations in (A, B), this class first solves
    one equation for A in terms of B and CMinfo, thus obtaining A(B,
    CMinfo). Then it substitutes A(B, CMinfo) into the remaining equation to
    obtain B(CMinfo). Finally, it substitutes B(CMinfo) into A(B, CMinfo) to
    get A(B(CMinfo), CMinfo).
    
    self.alpha_list_with_betas and self.alpha_mat_with_betas are used to store 
    A(B, CM_info).
    
    self.beta_list and self.beta_mat are used to store B(CM_info).
    
    self.alpha_list and self.alpha_mat are used to store A(B(CM_info),
    CM_info).
    
    The value of random variable x at time n will be denoted by x^{[ n]}. We
    will also use the notation

    \Delta x^{[n]} = x^{[n+1]}- x^{[n]}

    Set "delta=False" if you want 2-times correlations < x_i^{[n]},
    x_j^{[n+1]}> in the final result to be expressed as themselves. Set
    "delta=True" (recommended) if you want 2-times correlations < x_i^{[
    n]}, x_j^{[n+1]}> in the final result to be replaced by 2 terms,
    using the identity

    < x_i^{[n]}, x_j^{[n+1]}>= < x_i^{[n]}, x_j^{[n]}> + < x_i^{[n]},
    \Delta x_j^{[n]}>


    Attributes
    ----------
    # alpha_list and alpha_mat are inherited from parent class
    alpha_list_with_betas: list[sp.Equality]
    alpha_mat_with_betas: sp.Matrix
    beta_list: list[sp.Equality]
    beta_mat: sp.Matrix
    delta: bool

    """

    def __init__(self, graph, delta=True):
        """
        Constructor

        Parameters
        ----------
        graph: FBackGraph
        delta: bool
        """
        GainsCalculator.__init__(self, graph)
        self.delta = delta

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

        The inputs to cov_mat_list_in and mat_K do not matter as these
        variables are reassigned internally.

        Parameters
        ----------
        cov_mat_list_in: list[sp.Matrix]
        mat_K: sp.Matrix
        time: None or str or int

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
            d_cov2times = cov2times_sb_mat(dim, time=time0, delta=True)
        else:
            cov_mat0, cov2times, cov_mat1 = cov_mat_list_in
            d_cov2times = cov2times - cov_mat0

        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))
        mat_K = mat_B * cov_mat0

        calc = GainsCalculator(self.graph)
        calc.calculate_gains(cov_mat_in=cov_mat1, mat_K=mat_K, time=time1)
        self.alpha_mat_with_betas = deepcopy(calc.alpha_mat)
        self.alpha_list_with_betas = deepcopy(calc.alpha_list)

        self.alpha_mat = None
        self.alpha_list = None

        self.calculate_betas(cov_mat0, cov2times, d_cov2times, time=time0)
        self.calculate_alphas()

    def calculate_betas(self, cov_mat0, cov2times, d_cov2times, time):
        """
        This method fills in self.beta_list and self.beta_mat. It's an
        internal method called by calculate_gains().

        Parameters
        ----------
        cov_mat0: sp.Matrix
        cov2times: sp.Matrix
        d_cov2times: sp.Matrix
        time: None or str or int


        Returns
        -------
        None

        """
        dim = self.graph.num_nds

        mat_B = set_to_zero_fback_gains_without_arrows(self.graph,
                                             beta_sb_mat(dim))

        eq_list = eq_list0 = eq_list1 = None
        eq_mat = eq_mat0 = eq_mat1 = None
        unknowns = unknowns0 = unknowns1 = None

        if not self.delta:
            eq_mat = (sp.eye(dim) - self.alpha_mat_with_betas) * \
                     cov2times.T - \
                     mat_B * cov_mat0
        else:
            # delta method consists in solving
            # b= A x
            # by solving two systems of linear  equations:
            # b_0 = A x_0 solved for x_0
            # b_1 = A x_1 solved for x_1
            # so
            # b= b_0 + b_1 = A(x_0 + x_1)

            eq_mat0 = (sp.eye(dim) - self.alpha_mat_with_betas) * \
                     cov_mat0.T - mat_B * cov_mat0
            eq_mat1 = (sp.eye(dim) - self.alpha_mat_with_betas) * \
                      d_cov2times.T - mat_B * cov_mat0
        if not self.delta:
            unknowns = []
            eq_list = []
        else:
            unknowns0 = []
            unknowns1 = []
            eq_list0 = []
            eq_list1 = []
        for row, col in product(range(dim), range(dim)):
            row_nd = self.graph.ord_nodes[row]
            col_nd = self.graph.ord_nodes[col]
            if not self.delta:
                eq_list.append(eq_mat[row, col])
            else:
                eq_list0.append(eq_mat0[row, col])
                eq_list1.append(eq_mat1[row, col])
            if (col_nd, row_nd) in self.graph.fback_arrows:
                beta_str = "beta_" + str(row) + "_L_" + str(col)
                if not self.delta:
                    unknowns.append(sp.Symbol(beta_str))
                else:
                    unknowns0.append(sp.Symbol(beta_str))
                    unknowns1.append(sp.Symbol(beta_str))
            else:
                if time == "n":
                    xtra_str = "n"
                elif isinstance(time, int):
                    xtra_str = "n" + str(time)
                else:
                    assert False
                cov2times_str = "cov2times_" + xtra_str +\
                                 "_" + str(col) + "_" + str(row)
                if not self.delta:
                    unknowns.append(sp.Symbol(cov2times_str))
                else:
                    # do nothing for unknowns0
                    unknowns1.append(sp.Symbol("d_" + cov2times_str))

        # the comma does what is called sequence unpacking.
        # draws out item from single item list
        if not self.delta:
            sol_list, = linsolve(eq_list, unknowns)
            # print(str(sol_list))
            sol_list = sp.factor(sol_list)
        else:
            sol_list0, = linsolve(eq_list0, unknowns0)
            sol_list1, = linsolve(eq_list1, unknowns1)
            sol_list = []
            unknowns = unknowns1
            for i in range(len(sol_list1)):
                if i < len(unknowns0):
                    sol_list.append(sol_list0[i] + sol_list1[i])
                else:
                    sol_list.append(sol_list1[i])

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
        This method fills in self.alpha_list and self.alpha_mat. It's an
        internal method called by calculate_gains().

        Returns
        -------
        None

        """
        dim = self.graph.num_nds
        self.alpha_mat = deepcopy(self.alpha_mat_with_betas)
        self.alpha_list = deepcopy(self.alpha_list_with_betas)
        for row, col in product(range(dim), range(dim)):
            beta_str = "beta_" + str(row) + "_L_" + str(col)
            self.alpha_mat = sp.simplify(
                self.alpha_mat.subs(sp.Symbol(beta_str),
                self.beta_mat[row, col]))
            for i in range(len(self.alpha_list)):
                self.alpha_list[i] = sp.simplify(
                    self.alpha_list[i].subs(sp.Symbol(beta_str),
                            self.beta_mat[row, col]))
        # print("ccvvf", self.alpha_mat)

    def print_alpha_list_with_betas(self, verbose=False, time="n"):
        """
        This method prints the info in self.alpha_list_with_betas. It does
        this by calling latexify:print_list_sb().

        Parameters
        ----------
        verbose: bool
        time: None or str or int

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
        time: None or str or int

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
        time: None or str or int

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
