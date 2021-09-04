import dataclasses

from snakelib.custom_types import *
from snakelib.utils import (
    NewickTree,
    get_parameter_value,
    read_file_contents,
    parse_newick_string
)

from snakelib.program_parser import Program
from snakelib.iqtree_statstest_parser import get_iqtree_results

from iqtree_utils import (
    get_all_iqtree_llhs,
    get_best_iqtree_llh,
    get_iqtree_treesearch_cpu_time_entire_run,
    get_iqtree_cpu_time,
    get_iqtree_run_param_values_from_file,
)

def create_iqtree(
    parameter_file_path: FilePath,
    treesearch_log_file_path: FilePath,
    eval_log_file_path: FilePath,
    best_tree_file_path: FilePath,
    all_treesearch_trees_file_path: FilePath,
    best_eval_tree_file_path: FilePath,
    all_eval_trees_file_path: FilePath,
):
    treesearch_file_content = read_file_contents(all_treesearch_trees_file_path)
    treesearch_trees = [parse_newick_string(newick_str) for newick_str in treesearch_file_content]

    eval_file_content = read_file_contents(all_eval_trees_file_path)
    eval_trees = [parse_newick_string(newick_str) for newick_str in eval_file_content]

    # fmt: off
    iqtree = Program(
        blmin                   = get_parameter_value(parameter_file_path, "blmin"),
        blmax                   = get_parameter_value(parameter_file_path, "blmax"),
        lh_epsilon              = get_parameter_value(parameter_file_path, "lh_eps"),
        model_param_epsilon     = get_parameter_value(parameter_file_path, "model_param_epsilon"),
        branch_length_smoothing = None,
        spr_lh_epsilon          = None,
        bfgs_factor             = None,

        num_pars_trees          = int(get_parameter_value(parameter_file_path, "num_pars_trees")),
        num_rand_trees          = 0,
        best_treesearch_llh     = get_best_iqtree_llh(treesearch_log_file_path),
        best_evaluation_llh     = get_best_iqtree_llh(eval_log_file_path),
        treesearch_total_time   = get_iqtree_treesearch_cpu_time_entire_run(treesearch_log_file_path),

        # Tree search
        best_treesearch_tree        = parse_newick_string(read_file_contents(best_tree_file_path)[0]),
        treeseach_seeds             = get_iqtree_run_param_values_from_file(treesearch_log_file_path, "seed"),
        treesearch_trees            = treesearch_trees,
        treesearch_llhs             = get_all_iqtree_llhs(treesearch_log_file_path),
        treesearch_compute_times    = get_iqtree_cpu_time(treesearch_log_file_path),

        # Eval
        best_eval_tree              = parse_newick_string(read_file_contents(best_eval_tree_file_path)[0]),
        eval_blmins                 = get_iqtree_run_param_values_from_file(eval_log_file_path, "blmin"),
        eval_blmaxs                 = get_iqtree_run_param_values_from_file(eval_log_file_path, "blmax"),
        eval_lh_epsilons            = get_iqtree_run_param_values_from_file(eval_log_file_path, "eps"),
        eval_model_param_epsilons   = get_iqtree_run_param_values_from_file(eval_log_file_path, "me"),
        eval_raxml_brlen_smoothings = None,
        eval_spr_lh_epsilons        = None,
        eval_bfgs_factors           = None,

        eval_trees          = eval_trees,
        eval_llhs           = get_all_iqtree_llhs(eval_log_file_path),
        eval_compute_times  = get_iqtree_cpu_time(eval_log_file_path),
    )
    # fmt: on

    return iqtree

@dataclasses.dataclass
class Experiment:
    # fmt: off
    best_eval_trees                  : RunIndexed[NewickString]
    best_overall_eval_tree      : NewickString
    iqtree_statstests_results   : TreeIndexed[IqTreeMetrics]
    # fmt: on

    def _get_idx_for_newick(self, newick_str: NewickString) -> int:
        return self.best_eval_trees.index(newick_str)

    def eval_tree_is_overall_best(self, newick_str: NewickString) -> bool:
        return newick_str == self.best_overall_eval_tree

    def get_iqtree_llh_for_eval_tree(self, newick_str: NewickString) -> float:
        i = self._get_idx_for_newick(newick_str)
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["logL"]

    def get_iqtree_deltaL_for_eval_tree(self, newick_str: NewickString) -> float:
        i = self._get_idx_for_newick(newick_str)
        results_for_tree_index = self.iqtree_statstests_results[i]
        return results_for_tree_index["deltaL"]

    def get_iqtree_test_results_for_eval_tree(self, newick_str: NewickString) -> Dict:
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
