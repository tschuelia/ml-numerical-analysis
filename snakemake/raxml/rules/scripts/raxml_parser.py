import dataclasses

from snakelib.custom_types import *
from snakelib.utils import get_parameter_value, read_file_contents

# from iqtree_statstest_parser import get_iqtree_results
from raxml_utils import (
    get_all_raxml_llhs,
    get_all_raxml_seeds,
    get_best_raxml_llh,
    get_cleaned_rf_dist,
    get_raxml_abs_rf_distance,
    get_raxml_num_unique_topos,
    get_raxml_rel_rf_distance,
    get_raxml_run_param_values_from_file,
    get_raxml_elapsed_time,
    get_raxml_treesearch_elapsed_time_entire_run,
    read_rfdistances,
)


@dataclasses.dataclass
class Raxml:
    # Raxmlng stuff
    blmin: float
    blmax: float
    num_pars_trees: int
    num_rand_trees: int
    best_treesearch_llh: float
    best_evaluation_llh: float
    treesearch_total_time: float
    avg_abs_rfdist_treesearch: float
    avg_rel_rfdist_treesearch: float
    num_unique_topos_treesearch: int

    # RaxmlTreesearchTree stuff
    best_tree_newick: Newick
    treeseach_seeds: TreeIndexed[int]
    treesearch_trees: TreeIndexed[Newick]
    treesearch_llhs: TreeIndexed[float]
    treesearch_compute_times: TreeIndexed[float]
    # iqtree_statstests_results: TreeIndexed[IqTreeMetrics]

    # RaxmlEvalTree stuff
    best_eval_tree_newick: Newick
    eval_blmins: TreeIndexed[float]
    eval_blmaxs: TreeIndexed[float]
    eval_trees: TreeIndexed[Newick]
    eval_llhs: TreeIndexed[float]
    eval_compute_times: TreeIndexed[float]

    # RFDistTreesearchTree stuff
    all_treesearch_tree_rfdists: TreeTreeIndexed

    def get_num_of_trees(self) -> int:
        return self.num_pars_trees + self.num_rand_trees

    def tree_for_index_is_best(self, i: TreeIndex) -> bool:
        return self.get_newick_tree_for_tree_index(i) == self.best_tree_newick

    def get_treesearch_llh_for_tree_index(self, i: TreeIndex) -> float:
        return self.treesearch_llhs[i]

    def get_treesearch_compute_time_for_tree_index(self, i: TreeIndex) -> float:
        return self.treesearch_compute_times[i]

    def get_newick_tree_for_tree_index(self, i: TreeIndex) -> Newick:
        return self.treesearch_trees[i]

    def get_treesearch_seed_for_tree_index(self, i: TreeIndex) -> int:
        return self.treeseach_seeds[i]

    """ def get_iqtree_llh_for_tree_index(self, i: TreeIndex) -> float:
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["logL"]

    def get_iqtree_deltaL_for_tree_index(self, i: TreeIndex) -> float:
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["deltaL"]

    def get_iqtree_test_results_for_tree_index(self, i: TreeIndex) -> Dict:
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["tests"] """

    def get_plain_rfdist_for_trees(self, tree_indices: Tuple[TreeIndex, TreeIndex]):
        plain_rfdist, _ = self.all_treesearch_tree_rfdists[tree_indices]
        return plain_rfdist

    def get_normalized_rfdist_for_trees(
        self, tree_indices: Tuple[TreeIndex, TreeIndex]
    ):
        _, norm_rfdist = self.all_treesearch_tree_rfdists[tree_indices]
        return norm_rfdist

    def get_num_of_eval_trees(self) -> int:
        return len(self.eval_trees)

    def get_eval_llh_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_llhs[i]

    def get_eval_compute_time_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_compute_times[i]

    def get_newick_eval_tree_for_tree_index(self, i: TreeIndex) -> Newick:
        return self.eval_trees[i]

    def eval_tree_for_index_is_best(self, i: TreeIndex) -> bool:
        return self.get_newick_eval_tree_for_tree_index(i) == self.best_eval_tree_newick

    def get_eval_blmin_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_blmins[i]

    def get_eval_blmax_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_blmaxs[i]


# fmt: off
def create_raxml(
    parameter_file_path         : FilePath,
    treesearch_log_file_path    : FilePath,
    eval_log_file_path          : FilePath,
    rfdist_log_file_path        : FilePath,
    best_tree_file_path         : FilePath,
    all_treesearch_trees_file_path: FilePath,
    #iqtree_statstest_results_file_path: FilePath,
    best_eval_tree_file_path    : FilePath,
    command                     : str,
    all_eval_trees_file_path    : FilePath,
    rfdistances_file_path       : FilePath,
):  
    return Raxml(
        # Raxmlng stuff
        blmin                       = get_parameter_value(parameter_file_path, "blmin"),
        blmax                       = get_parameter_value(parameter_file_path, "blmax"),
        num_pars_trees              = get_parameter_value(parameter_file_path, "num_pars_trees"),
        num_rand_trees              = get_parameter_value(parameter_file_path, "num_rand_trees"),
        best_treesearch_llh         = get_best_raxml_llh(treesearch_log_file_path),
        best_evaluation_llh         = get_best_raxml_llh(eval_log_file_path),
        treesearch_total_time       = get_raxml_treesearch_elapsed_time_entire_run(treesearch_log_file_path),
        avg_abs_rfdist_treesearch   = get_raxml_abs_rf_distance(rfdist_log_file_path),
        avg_rel_rfdist_treesearch   = get_raxml_rel_rf_distance(rfdist_log_file_path),
        num_unique_topos_treesearch = get_raxml_num_unique_topos(rfdist_log_file_path),

        # RaxmlTreesearchTree stuff
        best_tree_newick            = read_file_contents(best_tree_file_path)[0],
        treeseach_seeds             = get_all_raxml_seeds(treesearch_log_file_path),
        treesearch_trees            = read_file_contents(all_treesearch_trees_file_path),
        treesearch_llhs             = get_all_raxml_llhs(treesearch_log_file_path),
        treesearch_compute_times    = get_raxml_elapsed_time(treesearch_log_file_path),
        # iqtree_statstests_results   = get_iqtree_results(iqtree_statstest_results_file_path),

        # RaxmlEvalTree stuff
        best_eval_tree_newick   = read_file_contents(best_eval_tree_file_path)[0],
        eval_blmins             = get_raxml_run_param_values_from_file(eval_log_file_path, command, "blmin"),
        eval_blmaxs             = get_raxml_run_param_values_from_file(eval_log_file_path, command, "blmax"),
        eval_trees              = read_file_contents(all_eval_trees_file_path),
        eval_llhs               = get_all_raxml_llhs(eval_log_file_path),
        eval_compute_times      = get_raxml_elapsed_time(eval_log_file_path),

        # RFDistTreesearchTree stuff
        all_treesearch_tree_rfdists = read_rfdistances(rfdistances_file_path),
    )
# fmt: on


@dataclasses.dataclass
class Experiment:
    # fmt: on
    raxml_best_trees            : RunIndexed[Newick]
    raxml_best_eval_trees       : RunIndexed[Newick]
    rfdist_raxml_best_trees     : RunRunIndexed
    rfdist_raxml_best_eval_trees: RunRunIndexed
    # fmt: off


    def get_plain_rfdist_for_raxml_trees(
        self, run_indices: Tuple[RunIndex, RunIndex]
    ) -> float:
        rf_dist_plain, _ = self.rfdist_raxml_best_trees[run_indices]
        return rf_dist_plain

    def get_normalized_rfdist_for_raxml_trees(
        self, run_indices: Tuple[RunIndex, RunIndex]
    ) -> float:
        _, rf_dist_norm = self.rfdist_raxml_best_trees[run_indices]
        return rf_dist_norm

    def get_plain_rfdist_for_raxml_eval_trees(
        self, run_indices: Tuple[RunIndex, RunIndex]
    ) -> float:
        rf_dist_plain, _ = self.rfdist_raxml_best_eval_trees[run_indices]
        return rf_dist_plain

    def get_normalized_rfdist_for_raxml_eval_trees(
        self, run_indices: Tuple[RunIndex, RunIndex]
    ) -> float:
        _, rf_dist_norm = self.rfdist_raxml_best_eval_trees[run_indices]
        return rf_dist_norm



# fmt: off
def create_Experiment(
    raxml_best_trees_path       : FilePath,
    raxml_best_eval_trees_path  : FilePath,
    rfdist_raxml_best_trees_path      : FilePath,
    rfdist_raxml_best_eval_trees_path : FilePath,
):
    return Experiment(
        raxml_best_trees              = read_file_contents(raxml_best_trees_path),
        raxml_best_eval_trees         = read_file_contents(raxml_best_eval_trees_path),
        rfdist_raxml_best_trees       = read_rfdistances(rfdist_raxml_best_trees_path),
        rfdist_raxml_best_eval_trees  = read_rfdistances(rfdist_raxml_best_eval_trees_path),
    )
# fmt: on
