from FBackGainsEstimator import *

class FBackGEmanager:

    def __init__(self, n_max,
                 graph,
                 path,
                 solve_symbolically=False,
                 hidden_nds=None):
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
                hidden_nds
            )
            self.n_to_estimator[time].calculate_gains()

    def print_alpha_lists(self, true_alpha_mat=None, verbose=False):
        str0 = r"\begin{array}{l}" + "\n"
        for time in range(1, self.n_max):
            str0 += r"\text{******** time=" + str(time) + r"}\\" + "\n"
            gest = self.n_to_estimator[time]
            x = gest.print_alpha_list(true_alpha_mat=true_alpha_mat,
                                  verbose=verbose)
            str0 += str(x) + "\n" + r"\\" + "\n"
        str0 = str0[:-3]
        str0 += r"\end{array}"
        return sp.Symbol(str0)
    
    def print_beta_lists(self, true_beta_mat=None, verbose=False):
        str0 = r"\begin{array}{l}" + "\n"
        for time in range(1, self.n_max):
            str0 += r"\text{******** time=" + str(time) + r"}\\" + "\n"
            gest = self.n_to_estimator[time]
            x = gest.print_beta_list(true_beta_mat=true_beta_mat,
                                  verbose=verbose)
            str0 += str(x) + "\n" + r"\\" + "\n"
        str0 = str0[:-3]
        str0 += r"\end{array}"
        return sp.Symbol(str0)