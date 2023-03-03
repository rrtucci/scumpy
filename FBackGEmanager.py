import numpy as np

from FBackGainsEstimator import *

class FBackGEmanager:
    """

    Attributes
    ----------
    graph: FBackGraph
    hidden_nds:
    n_max: int
    n_to_estimator:
    solve_symbolically: bool

    """

    def __init__(self, n_max,
                 graph,
                 path,
                 solve_symbolically=False,
                 hidden_nds=None,
                 delta=True):
        """

        Parameters
        ----------
        n_max
        graph
        path
        solve_symbolically
        hidden_nds
        delta
        """
        self.n_max = n_max
        self.graph = graph
        self.solve_symbolically = solve_symbolically
        if hidden_nds is None:
            self.hidden_nds = []
        else:
            assert set(hidden_nds).issubset(graph.ord_nodes)
            self.hidden_nds = hidden_nds

        df = pd.read_csv(path)
        columns = FBackRandomDataMaker.get_columns(n_max, graph)
        assert set(df.columns) == set(columns)
        # put columns in same order as graph.ord_nodes
        df = df[columns]
        dim = self.graph.num_nds
        self.n_to_estimator = {}
        for time in range(1, self.n_max):
            slice_n = columns[(time-1)*dim : time*dim]
            slice_n_plus_one = columns[time * dim: (time+1) * dim]
            two_slices = slice_n + slice_n_plus_one
            df_two_slices = df[two_slices]
            self.n_to_estimator[time] = FBackGainsEstimator(
                time,
                graph,
                df_two_slices,
                solve_symbolically,
                hidden_nds,
                delta
            )
            self.n_to_estimator[time].calculate_gains()

    def print_greek_lists(self, name, true_greek_mat=None, verbose=False):
        """

        Parameters
        ----------
        name
        true_greek_mat
        verbose

        Returns
        -------

        """
        assert name in ["alpha", "beta"]
        str0 = r"\begin{array}{l}" + "\n"
        for time in range(1, self.n_max):
            str0 += r"\text{******** time=" + str(time) + r"}\\" + "\n"
            gest = self.n_to_estimator[time]
            if name == "alpha":
                mat = gest.alpha_mat_estimate
            else:
                mat = gest.beta_mat_estimate
            eq_list = create_eq_list_from_matrix(mat, name, self.graph,
                                                 time=None)
            comments = gest.get_greek_list_comments(
                name , eq_list, true_greek_mat=true_greek_mat)
            x = print_list_sb(eq_list,
                              self.graph,
                              comment_list=comments)
            str0 += str(x) + "\n" + r"\\" + "\n"
        str0 = str0[:-3]
        str0 += r"\end{array}"
        if verbose:
            print(str0)
        return sp.Symbol(str0)

    def print_mean_greek_list(self, name, true_greek_mat=None,
                              verbose=False):
        """

        Parameters
        ----------
        name
        true_greek_mat
        verbose

        Returns
        -------

        """
        assert name in ["alpha", "beta"]
        if name == "alpha":
            av_np = np.mean([self.n_to_estimator[
                               time].alpha_mat_estimate for
                           time in range(1, self.n_max)], axis=0)
        else:
            av_np = np.mean([self.n_to_estimator[
                                 time].beta_mat_estimate for
                             time in range(1, self.n_max)], axis=0)
        mat =  sp.Matrix(av_np)
        eq_list = create_eq_list_from_matrix(mat, name, self.graph,
                                             time=None)
        comments = self.n_to_estimator[1].get_greek_list_comments(
            name, eq_list, true_greek_mat=true_greek_mat)

        return print_list_sb(eq_list, self.graph,
                             comment_list=comments, verbose=verbose)

    
    def print_alpha_lists(self, true_alpha_mat=None, verbose=False):
        """

        Parameters
        ----------
        true_alpha_mat
        verbose

        Returns
        -------

        """
        return self.print_greek_lists("alpha", 
                                      true_greek_mat=true_alpha_mat,
                                      verbose=verbose)
    
    def print_mean_alpha_list(self, true_alpha_mat=None, verbose=False):
        """

        Parameters
        ----------
        true_alpha_mat
        verbose

        Returns
        -------

        """
        return self.print_mean_greek_list("alpha", 
                                      true_greek_mat=true_alpha_mat,
                                      verbose=verbose)

    def print_beta_lists(self, true_beta_mat=None, verbose=False):
        """

        Parameters
        ----------
        true_beta_mat
        verbose

        Returns
        -------

        """
        return self.print_greek_lists("beta",
                                      true_greek_mat=true_beta_mat,
                                      verbose=verbose)

    def print_mean_beta_list(self, true_beta_mat=None, verbose=False):
        """

        Parameters
        ----------
        true_beta_mat
        verbose

        Returns
        -------

        """
        return self.print_mean_greek_list("beta",
                                          true_greek_mat=true_beta_mat,
                                          verbose=verbose)

if __name__ == "__main__":
    def main():
        path = 'dot_atlas/fback-2node.dot'
        graph = FBackGraph(path)
        dim = graph.num_nds
        mean_eps = [0]*dim
        sig_eps = [10] * dim
        n_max = 4
        alpha_bound = 10
        beta_bound = 1
        alpha_mat, beta_mat = \
            FBackRandomDataMaker.generate_random_alpha_and_beta_mats(
            graph, alpha_bound=alpha_bound, beta_bound=beta_bound)
        dmaker = FBackRandomDataMaker(n_max, graph,
                                      alpha_mat=alpha_mat,
                                      beta_mat=beta_mat,
                                      mean_eps=mean_eps,
                                      sig_eps=sig_eps)
        num_rows = 100
        data_path = "test_data.csv"
        dmaker.write_dataset_csv(num_rows, data_path)
        df = pd.read_csv(data_path)
        print(df)
        mger = FBackGEmanager(n_max, graph, data_path)
        mger.print_alpha_lists(true_alpha_mat=dmaker.alpha_mat, verbose=True)
        mger.print_mean_alpha_list(true_alpha_mat=dmaker.alpha_mat,
                                  verbose=True)
        mger.print_beta_lists(true_beta_mat=dmaker.beta_mat, verbose=True)
        mger.print_mean_beta_list(true_beta_mat=dmaker.beta_mat,
                                  verbose=True)


    main()