import dataclasses

from snakelib.custom_types import *
from snakelib.utils import (
    get_parameter_value,
    read_file_contents,
    get_number_of_taxa_for_tree,
    get_total_branch_length_for_tree,
    get_average_branch_length_for_tree
)

from snakelib.iqtree_statstest_parser import get_iqtree_results

from iqtree_utils import (
    get_all_iqtree_llhs,
    get_best_iqtree_llh,
    get_iqtree_treesearch_cpu_time_entire_run,
    get_iqtree_cpu_time,
    get_iqtree_run_param_values_from_file,
)


@dataclasses.dataclass
class Iqtree:
    blmin: float
    blmax: float
    model_param_epsilon: float
    lh_epsilon: float

    num_pars_trees: int
    # num_rand_trees: int
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
    eval_model_param_epsilon: TreeIndexed[float]
    eval_lh_epsilon: TreeIndexed[float]

    eval_trees: TreeIndexed[float]
    eval_llhs: TreeIndexed[float]
    eval_compute_times: TreeIndexed[float]

    def get_num_of_trees(self) -> int:
        return self.num_pars_trees  # + self.num_rand_trees

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

    def get_number_of_taxa_for_tree_index(self, i: TreeIndex) -> int:
        newick_str = self.get_newick_tree_for_tree_index(i)
        return get_number_of_taxa_for_tree(newick_str)

    def get_total_branch_length_for_tree_index(self, i: TreeIndex) -> float:
        newick_str = self.get_newick_tree_for_tree_index(i)
        return get_total_branch_length_for_tree(newick_str)

    def get_average_branch_length_for_tree_index(self, i: TreeIndex) -> float:
        newick_str = self.get_newick_tree_for_tree_index(i)
        return get_average_branch_length_for_tree(newick_str)

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

    def get_number_of_taxa_for_eval_tree_index(self, i: TreeIndex) -> int:
        newick_str = self.get_newick_eval_tree_for_tree_index(i)
        return get_number_of_taxa_for_tree(newick_str)

    def get_total_branch_length_for_eval_tree_index(self, i: TreeIndex) -> float:
        newick_str = self.get_newick_eval_tree_for_tree_index(i)
        return get_total_branch_length_for_tree(newick_str)

    def get_average_branch_length_for_eval_tree_index(self, i: TreeIndex) -> float:
        newick_str = self.get_newick_eval_tree_for_tree_index(i)
        return get_average_branch_length_for_tree(newick_str)

    def get_eval_blmin_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_blmins[i]

    def get_eval_blmax_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_blmaxs[i]

    def get_eval_model_param_epsilon_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_model_param_epsilon[i]

    def get_eval_lh_epsilon_for_tree_index(self, i: TreeIndex) -> float:
        return self.eval_lh_epsilon[i]


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
        blmin                   = get_parameter_value(parameter_file_path, "blmin"),
        blmax                   = get_parameter_value(parameter_file_path, "blmax"),
        model_param_epsilon     = get_parameter_value(parameter_file_path, "model_param_epsilon"),
        lh_epsilon              = get_parameter_value(parameter_file_path, "lh_eps"),

        num_pars_trees          = get_parameter_value(parameter_file_path, "num_pars_trees"),
        #num_rand_trees          = get_parameter_value(parameter_file_path, "num_rand_trees"),
        best_treesearch_llh     = get_best_iqtree_llh(treesearch_log_file_path),
        best_evaluation_llh     = get_best_iqtree_llh(eval_log_file_path),
        treesearch_total_time   = get_iqtree_treesearch_cpu_time_entire_run(treesearch_log_file_path),
        
        # IqtreeTreesearchTree stuff
        best_tree_newick    = read_file_contents(best_tree_file_path)[0],
        treesearch_seeds    = get_iqtree_run_param_values_from_file(treesearch_log_file_path, "seed"),
        treesearch_trees    = read_file_contents(all_treesearch_trees_file_path),
        treesearch_llhs     = get_all_iqtree_llhs(treesearch_log_file_path),
        treesearch_compute_times = get_iqtree_cpu_time(treesearch_log_file_path),
        
        # IqtreeEvalTree
        best_eval_tree_newick = read_file_contents(best_eval_tree_file_path)[0],
        eval_blmins         = get_iqtree_run_param_values_from_file(eval_log_file_path, "blmin"),
        eval_blmaxs         = get_iqtree_run_param_values_from_file(eval_log_file_path, "blmax"),
        eval_model_param_epsilon = get_iqtree_run_param_values_from_file(eval_log_file_path, "me"),
        eval_lh_epsilon          =get_iqtree_run_param_values_from_file(eval_log_file_path, "eps"),
        eval_trees          = read_file_contents(all_eval_trees_file_path),
        eval_llhs           = get_all_iqtree_llhs(eval_log_file_path),
        eval_compute_times  = get_iqtree_cpu_time(eval_log_file_path),
    )
    # fmt: on

@dataclasses.dataclass
class Experiment:
    # fmt: off
    best_eval_trees                  : RunIndexed[Newick]
    best_overall_eval_tree      : Newick
    iqtree_statstests_results   : TreeIndexed[IqTreeMetrics]
    # fmt: on

    def _get_idx_for_newick(self, newick_str: Newick) -> int:
        return self.best_eval_trees.index(newick_str)

    def eval_tree_is_overall_best(self, newick_str: Newick) -> bool:
        return newick_str == self.best_overall_eval_tree

    def get_iqtree_llh_for_eval_tree(self, newick_str: Newick) -> float:
        i = self._get_idx_for_newick(newick_str)
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["logL"]

    def get_iqtree_deltaL_for_eval_tree(self, newick_str: Newick) -> float:
        i = self._get_idx_for_newick(newick_str)
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["deltaL"]

    def get_iqtree_test_results_for_eval_tree(self, newick_str: Newick) -> Dict:
        i = self._get_idx_for_newick(newick_str)
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["tests"]

# fmt: off
def create_Experiment(
        best_trees_file_path                : FilePath,
        best_overall_eval_tree_file_path    : FilePath,
        iqtree_statstest_results_file_path  : FilePath,
):
    return Experiment(
        best_eval_trees=read_file_contents(best_trees_file_path),
        best_overall_eval_tree=read_file_contents(best_overall_eval_tree_file_path)[0],
        iqtree_statstests_results=get_iqtree_results(iqtree_statstest_results_file_path),
    )
# fmt: on
