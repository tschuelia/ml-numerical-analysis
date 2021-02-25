import dataclasses

from snakelib.custom_types import *
from snakelib.utils import get_parameter_value, read_file_contents

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


def create_fasttree(
        parameter_file_path: FilePath,
        treesearch_log_file_path: FilePath,
        best_tree_file_path: FilePath,
        all_treesearch_trees_file_path: FilePath,
):
    # fmt: off
    return Fasttree(
        blmin=get_parameter_value(parameter_file_path, "blmin"),

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
