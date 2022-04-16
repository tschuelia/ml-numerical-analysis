import dataclasses

from .custom_types import *

from .utils import (
    NewickTree,
)

@dataclasses.dataclass
class Program:
    # Runs
    lh_epsilon: Optional[float]

    num_pars_trees: Optional[int]
    num_rand_trees: Optional[int]
    best_treesearch_llh: float
    best_evaluation_llh: Optional[float]
    treesearch_total_time: float

    # Tree search
    best_treesearch_tree: NewickTree
    treeseach_seeds: TreeIndexed[int]
    treesearch_trees: TreeIndexed[NewickTree]
    treesearch_llhs: TreeIndexed[float]
    treesearch_compute_times: TreeIndexed[float]

    # Eval
    best_eval_tree: Optional[NewickTree]

    eval_trees: Optional[TreeIndexed[NewickTree]]
    eval_llhs: Optional[TreeIndexed[float]]
    eval_compute_times: Optional[TreeIndexed[float]]

    def __post_init__(self):
        # check number of trees
        if not self.num_pars_trees and not self.num_rand_trees:
            raise ValueError("Set either num_pars_trees or num_rand_trees.")

        # some sanity checks
        assert len(self.treeseach_seeds) == len(self.treesearch_trees) == len(self.treesearch_llhs) == len(self.treesearch_compute_times)

        if self.eval_trees:
            assert len(self.eval_trees) == len(self.eval_llhs) == len(self.eval_compute_times)

    def get_number_of_trees(self) -> int:
        num_pars = self.num_pars_trees if self.num_pars_trees else 0
        num_rand = self.num_rand_trees if self.num_rand_trees else 0
        return num_pars + num_rand

    def tree_for_index_is_best(self, i: TreeIndex) -> bool:
        return self.treesearch_trees[i].newick_str == self.best_treesearch_tree.newick_str

    def get_number_of_eval_trees(self) -> int:
        return len(self.eval_trees)

    def eval_tree_for_index_is_best(self, i: TreeIndex) -> bool:
        return self.eval_trees[i].newick_str == self.best_eval_tree.newick_str
