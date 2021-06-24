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

from fasttree_utils import (
    get_all_fasttree_llhs,
    get_best_fasttree_llh,
    get_fasttree_treesearch_entire_run,
    get_fasttree_runtimes,
    get_fasttree_run_param_values_from_file,
)


@dataclasses.dataclass
class Fasttree:
    blmin: float
    lh_eps: float

    num_trees: int
    best_treesearch_llh: float
    treesearch_total_time: float

    # fasttreeTreesearchTree stuff
    best_tree_newick: Newick
    treesearch_seeds: TreeIndexed[int]
    treesearch_trees: TreeIndexed[Newick]
    treesearch_llhs: TreeIndexed[float]
    treesearch_compute_times: TreeIndexed[float]

    def get_num_of_trees(self) -> int:
        return self.num_trees

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


def create_fasttree(
        parameter_file_path: FilePath,
        treesearch_log_file_path: FilePath,
        best_tree_file_path: FilePath,
        all_treesearch_trees_file_path: FilePath,
):
    # fmt: off
    return Fasttree(
        blmin=get_parameter_value(parameter_file_path, "blmin"),
        lh_eps=get_parameter_value(parameter_file_path, "lh_eps"),

        num_trees=get_parameter_value(parameter_file_path, "num_trees"),

        best_treesearch_llh=get_best_fasttree_llh(treesearch_log_file_path),
        treesearch_total_time=get_fasttree_treesearch_entire_run(treesearch_log_file_path),

        # FasttreeTreesearchTree stuff
        best_tree_newick=read_file_contents(best_tree_file_path)[0],
        treesearch_seeds=get_fasttree_run_param_values_from_file(treesearch_log_file_path, "seed"),
        treesearch_trees=read_file_contents(all_treesearch_trees_file_path),
        treesearch_llhs=get_all_fasttree_llhs(treesearch_log_file_path),
        treesearch_compute_times=get_fasttree_runtimes(treesearch_log_file_path),

    )
    # fmt: on

@dataclasses.dataclass
class Experiment:
    # fmt: off
    best_trees                  : RunIndexed[Newick]
    best_overall_eval_tree      : Newick
    iqtree_statstests_results   : TreeIndexed[IqTreeMetrics]
    # fmt: on

    def _get_idx_for_newick(self, newick_str: Newick) -> int:
        return self.best_trees.index(newick_str)

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
        best_trees=read_file_contents(best_trees_file_path),
        best_overall_eval_tree=read_file_contents(best_overall_eval_tree_file_path)[0],
        iqtree_statstests_results=get_iqtree_results(iqtree_statstest_results_file_path),
    )
# fmt: on