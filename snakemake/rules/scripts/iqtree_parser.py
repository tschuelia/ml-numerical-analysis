import dataclasses

from custom_types import *
from utils import (
    get_all_iqtree_llhs,
    get_parameter_value,
    get_best_iqtree_llh,
    read_file_contents,
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


def create_iqtree(
    parameter_file_path: FilePath,
    treesearch_log_file_path: FilePath,
    eval_log_file_path: FilePath,
    best_tree_file_path: FilePath,
    all_treesearch_trees_file_path: FilePath,
    best_eval_tree_file_path: FilePath,
    all_eval_trees_file_path: FilePath,
):
    return Iqtree(
        num_pars_trees=get_parameter_value(parameter_file_path, "num_pars_trees"),
        num_rand_trees=get_parameter_value(parameter_file_path, "num_rand_trees"),
        best_treesearch_llh=get_best_iqtree_llh(treesearch_log_file_path),
        best_evaluation_llh=get_best_iqtree_llh(eval_log_file_path),
        # TODO: iqtree compute times
        treesearch_total_time=-1,
        # IqtreeTreesearchTree stuff
        best_tree_newick=read_file_contents(best_tree_file_path)[0],
        # TODO: iqtree seeds
        treesearch_seeds=-1,
        treesearch_trees=read_file_contents(all_treesearch_trees_file_path),
        treesearch_llhs=get_all_iqtree_llhs(treesearch_log_file_path),
        # TODO: iqtree compute times
        treesearch_compute_times=-1,
        # IqtreeEvalTree
        best_eval_tree_newick=read_file_contents(best_eval_tree_file_path)[0],
        # TODO: get run params from iqtree file
        eval_blmins=-1,
        eval_blmaxs=-1,
        eval_trees=read_file_contents(all_eval_trees_file_path),
        eval_llhs=get_all_iqtree_llhs(eval_log_file_path),
        # TODO: eval compute times
        eval_compute_times=-1,
    )
