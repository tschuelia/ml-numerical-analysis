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

from raxml_utils import (
    get_all_raxml_llhs,
    get_all_raxml_seeds,
    get_best_raxml_llh,
    get_raxml_run_param_values_from_file,
    get_raxml_elapsed_time,
    get_raxml_treesearch_elapsed_time_entire_run,
    read_rfdistances,
)


def create_raxml(
        parameter_file_path: FilePath,
        treesearch_log_file_path: FilePath,
        eval_log_file_path: FilePath,
        best_tree_file_path: FilePath,
        all_treesearch_trees_file_path: FilePath,
        best_eval_tree_file_path: FilePath,
        raxml_command: str,
        all_eval_trees_file_path: FilePath,
):
    treesearch_trees_file_content = read_file_contents(all_treesearch_trees_file_path)
    treesearch_trees = [parse_newick_string(newick_str) for newick_str in treesearch_trees_file_content]

    eval_trees_file_content = read_file_contents(all_eval_trees_file_path)
    eval_trees = [parse_newick_string(newick_str) for newick_str in eval_trees_file_content]

    raxml = Program(
        blmin                   = get_parameter_value(parameter_file_path, "blmin"),
        blmax                   = get_parameter_value(parameter_file_path, "blmax"),
        lh_epsilon              = get_parameter_value(parameter_file_path, "lh_eps"),
        model_epsilon     = get_parameter_value(parameter_file_path, "model_epsilon"),
        branch_length_smoothing = int(get_parameter_value(parameter_file_path, "raxml_brlen_smoothings")),
        spr_lh_epsilon          = get_parameter_value(parameter_file_path, "spr_lh_epsilon"),
        bfgs_factor             = get_parameter_value(parameter_file_path, "bfgs_factor"),

        num_pars_trees          = int(get_parameter_value(parameter_file_path, "num_pars_trees")),
        num_rand_trees          = int(get_parameter_value(parameter_file_path, "num_rand_trees")),
        best_treesearch_llh     = get_best_raxml_llh(treesearch_log_file_path),
        best_evaluation_llh     = get_best_raxml_llh(eval_log_file_path),
        treesearch_total_time   = get_raxml_treesearch_elapsed_time_entire_run(treesearch_log_file_path),

        # Tree search
        best_treesearch_tree        = parse_newick_string(read_file_contents(best_tree_file_path)[0]),
        treeseach_seeds             = get_all_raxml_seeds(treesearch_log_file_path),
        treesearch_trees            = treesearch_trees,
        treesearch_llhs             = get_all_raxml_llhs(treesearch_log_file_path),
        treesearch_compute_times    = get_raxml_elapsed_time(treesearch_log_file_path),

        # Eval
        best_eval_tree              = parse_newick_string(read_file_contents(best_eval_tree_file_path)[0]),
        eval_blmins                 = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "blmin"),
        eval_blmaxs                 = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "blmax"),
        eval_lh_epsilons            = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "lh-epsilon"),
        eval_model_param_epsilons   = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "param-eps"),
        eval_raxml_brlen_smoothings = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "brlen-smoothings"),
        eval_spr_lh_epsilons        = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "spr-lheps"),
        eval_bfgs_factors           = get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "bfgs-factor"),

        eval_trees          = eval_trees,
        eval_llhs           = get_all_raxml_llhs(eval_log_file_path),
        eval_compute_times  = get_raxml_elapsed_time(eval_log_file_path),
    )
    return raxml
# fmt: on

@dataclasses.dataclass
class Experiment:
    # fmt: ff
    raxml_best_trees        : RunIndexed[NewickString]
    raxml_best_eval_trees   : RunIndexed[NewickString]
    rfdist_raxml_best_trees : RunRunIndexed
    rfdist_raxml_best_eval_trees    : RunRunIndexed
    #best_overall_eval_tree          : NewickString
    #iqtree_statstests_results       : TreeIndexed[IqTreeMetrics]

    # fmt: on

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

    def _get_idx_for_newick(self, newick_str: NewickString) -> int:
        return self.raxml_best_eval_trees.index(newick_str)

    # def eval_tree_is_overall_best(self, newick_str: NewickString) -> bool:
    #     return newick_str == self.best_overall_eval_tree
    #
    # def get_iqtree_llh_for_eval_tree(self, newick_str: NewickString) -> float:
    #     i = self._get_idx_for_newick(newick_str)
    #     results_for_tree_index = self.iqtree_statstests_results[i]
    #     return results_for_tree_index["logL"]
    #
    # def get_iqtree_deltaL_for_eval_tree(self, newick_str: NewickString) -> float:
    #     i = self._get_idx_for_newick(newick_str)
    #     results_for_tree_index = self.iqtree_statstests_results[i]
    #     return results_for_tree_index["deltaL"]
    #
    # def get_iqtree_test_results_for_eval_tree(self, newick_str: NewickString) -> Dict:
    #     i = self._get_idx_for_newick(newick_str)
    #     results_for_tree_index = self.iqtree_statstests_results[i]
    #     return results_for_tree_index["tests"]


# fmt: off
def create_Experiment(
        raxml_best_trees_path               : FilePath,
        raxml_best_eval_trees_path          : FilePath,
        rfdist_raxml_best_trees_path        : FilePath,
        rfdist_raxml_best_eval_trees_path   : FilePath,
        #best_overall_eval_tree_file_path    : FilePath,
        #iqtree_statstest_results_file_path  : FilePath,
):
    return Experiment(
        raxml_best_trees=read_file_contents(raxml_best_trees_path),
        raxml_best_eval_trees=read_file_contents(raxml_best_eval_trees_path),
        rfdist_raxml_best_trees=read_rfdistances(rfdist_raxml_best_trees_path),
        rfdist_raxml_best_eval_trees=read_rfdistances(rfdist_raxml_best_eval_trees_path),

        # Iqtree significance tests stuff
        #best_overall_eval_tree=read_file_contents(best_overall_eval_tree_file_path)[0],
        #iqtree_statstests_results=get_iqtree_results(iqtree_statstest_results_file_path),

    )
# fmt: on
