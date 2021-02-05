import dataclasses

from custom_types import *
from utils import (
    get_all_iqtree_llhs,
    get_parameter_value,
    get_best_iqtree_llh,
    read_file_contents,
    get_iqtree_treesearch_wallclock_time_entire_run,
    get_iqtree_wallclock_time,
    get_iqtree_run_param_values_from_file,
)


@dataclasses.dataclass
class Iqtree:
    num_pars_trees: int
    num_rand_trees: int
    best_treesearch_llh: float
    best_evaluation_llh: float
    treesearch_total_time: float

    # IqtreeTreesearchTree stuff
    best_tree_newick: Newick
    treesearch_seeds: TreeIndexed[int]
    treesearch_trees: TreeIndexed[Newick]
    treesearch_llhs: TreeIndexed[float]
    treesearch_compute_times: TreeIndexed[float]

    # IqtreeEvalTree
    best_eval_tree_newick: Newick
    eval_blmins: TreeIndexed[float]
    eval_blmaxs: TreeIndexed[float]
    eval_trees: TreeIndexed[float]
    eval_llhs: TreeIndexed[float]
    eval_compute_times: TreeIndexed[float]

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
        return self.treesearch_seeds[i]

    def get_num_of_eval_trees(self) -> int:
        return len(self.eval_trees)

    def get_eval_llh_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_llhs[i]

    def get_newick_eval_tree_for_tree_index(self, i: TreeIndex) -> Newick:
        return self.eval_trees[i]

    def get_eval_compute_time_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_compute_times[i]

    def eval_tree_for_index_is_best(self, i: TreeIndex) -> bool:
        return self.get_newick_eval_tree_for_tree_index(i) == self.best_eval_tree_newick

    def get_eval_blmin_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_blmins[i]

    def get_eval_blmax_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_blmaxs[i]


def create_iqtree(
    parameter_file_path: FilePath,
    treesearch_log_file_path: FilePath,
    eval_log_file_path: FilePath,
    best_tree_file_path: FilePath,
    all_treesearch_trees_file_path: FilePath,
    best_eval_tree_file_path: FilePath,
    all_eval_trees_file_path: FilePath,
):
    # fmt: off
    return Iqtree(
        num_pars_trees          = get_parameter_value(parameter_file_path, "num_pars_trees"),
        num_rand_trees          = get_parameter_value(parameter_file_path, "num_rand_trees"),
        best_treesearch_llh     = get_best_iqtree_llh(treesearch_log_file_path),
        best_evaluation_llh     = get_best_iqtree_llh(eval_log_file_path),
        treesearch_total_time   = get_iqtree_treesearch_wallclock_time_entire_run(treesearch_log_file_path),
        
        # IqtreeTreesearchTree stuff
        best_tree_newick    = read_file_contents(best_tree_file_path)[0],
        treesearch_seeds    = get_iqtree_run_param_values_from_file(treesearch_log_file_path, "seed"),
        treesearch_trees    = read_file_contents(all_treesearch_trees_file_path),
        treesearch_llhs     = get_all_iqtree_llhs(treesearch_log_file_path),
        treesearch_compute_times = get_iqtree_wallclock_time(treesearch_log_file_path),
        
        # IqtreeEvalTree
        best_eval_tree_newick = read_file_contents(best_eval_tree_file_path)[0],
        eval_blmins         = get_iqtree_run_param_values_from_file(eval_log_file_path, "blmin"),
        eval_blmaxs         = get_iqtree_run_param_values_from_file(eval_log_file_path, "blmax"),
        eval_trees          = read_file_contents(all_eval_trees_file_path),
        eval_llhs           = get_all_iqtree_llhs(eval_log_file_path),
        eval_compute_times  = get_iqtree_wallclock_time(eval_log_file_path),
    )
    # fmt: on
