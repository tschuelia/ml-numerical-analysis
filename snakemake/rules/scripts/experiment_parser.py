from typing import Dict, List, Tuple

from custom_types import *
from iqtree_parser import get_iqtree_results
from utils import (
    get_all_raxml_seeds,
    get_best_raxml_llh,
    get_best_iqtree_llh,
    get_cleaned_rf_dist,
    get_parameter_value,
    get_raxml_abs_rf_distance,
    get_all_raxml_llhs,
    get_raxml_num_unique_topos,
    get_raxml_rel_rf_distance,
    get_raxml_treesearch_elapsed_time,
)


class Run:
    def __init__(
        self,
        num_raxml_pars_trees: int,
        num_raxml_rand_trees: int,
        blmin: float,
        blmax: float,
        raxml_best_llh: float,
        iqtree_best_llh: float,
        raxml_best_tree: Newick,
        raxml_seeds: TreeIndexed[int],
        raxml_all_trees: TreeIndexed[Newick],
        raxml_llhs_all_trees: TreeIndexed[float],
        raxml_treesearch_elapsed_time: float,
        iqtree_all_trees: TreeIndexed[Newick],
        iqtree_results: TreeIndexed[IqTreeMetrics],
        average_absolute_rf_distance: float,
        average_relative_rf_distance: float,
        num_unique_topos: int,
        rfdist_all_trees: TreeTreeIndexed,
    ):
        self.num_raxml_pars_trees = num_raxml_pars_trees
        self.num_raxml_rand_trees = num_raxml_rand_trees
        self.blmin = blmin
        self.blmax = blmax
        self.raxml_best_llh = raxml_best_llh
        self.iqtree_best_llh = iqtree_best_llh
        self.raxml_best_tree = raxml_best_tree
        self.raxml_seeds = raxml_seeds
        self.raxml_all_trees = raxml_all_trees
        self.raxml_llhs_all_trees = raxml_llhs_all_trees
        self.raxml_treesearch_elapsed_time = raxml_treesearch_elapsed_time
        self.iqtree_all_trees = iqtree_all_trees
        self.iqtree_results = iqtree_results
        self.average_absolute_rf_distance = average_absolute_rf_distance
        self.average_relative_rf_distance = average_relative_rf_distance
        self.num_unique_topos = num_unique_topos
        self.rfdist_all_trees = rfdist_all_trees

    def get_num_raxml_pars_trees(self) -> int:
        return self.num_raxml_pars_trees

    def get_num_raxml_rand_trees(self) -> int:
        return self.num_raxml_rand_trees

    def get_blmin(self) -> float:
        return self.blmin

    def get_blmax(self) -> float:
        return self.blmax

    def get_best_raxml_llh(self) -> float:
        return self.raxml_best_llh

    def get_best_iqtree_llh(self) -> float:
        return self.iqtree_best_llh

    def get_average_absolute_rf_distance(self) -> float:
        return self.average_absolute_rf_distance

    def get_average_relative_rf_distance(self) -> float:
        return self.average_relative_rf_distance

    def get_num_unique_topos(self) -> int:
        return self.num_unique_topos

    def get_num_of_trees(self) -> int:
        return len(self.raxml_all_trees)

    def get_raxml_tree_for_tree_index(self, i: TreeIndex) -> Newick:
        return self.raxml_all_trees[i]

    def get_iqtree_tree_for_tree_index(self, i: TreeIndex) -> Newick:
        return self.iqtree_all_trees[i]

    def get_raxml_llh_for_tree_index(self, i: TreeIndex) -> float:
        return self.raxml_llhs_all_trees[i]

    def get_raxml_treesearch_elapsed_time(self) -> float:
        return self.raxml_treesearch_elapsed_time

    def get_iqtree_llh_for_tree_index(self, i: TreeIndex) -> float:
        results_for_tree_index = self.iqtree_results[i]
        return results_for_tree_index["logL"]

    def tree_for_index_is_best(self, i: TreeIndex) -> bool:
        return self.get_raxml_tree_for_tree_index(i) == self.raxml_best_tree

    def get_raxml_seed_for_tree_index(self, i: TreeIndex) -> bool:
        return self.raxml_seeds[i]

    def get_iqtree_deltaL_for_tree_index(self, i: TreeIndex) -> float:
        results_for_tree_index = self.iqtree_results[i]
        return results_for_tree_index["deltaL"]

    def get_iqtree_test_results_for_tree_index(self, i: TreeIndex) -> Dict:
        results_for_tree_index = self.iqtree_results[i]
        return results_for_tree_index["tests"]

    def get_plain_rfdistance_for_trees(self, tree_indices: Tuple[TreeIndex, TreeIndex]):
        rf_dist = self.rfdist_all_trees[tree_indices]
        return rf_dist[0]

    def get_normalized_rfdistance_for_trees(
        self, tree_indices: Tuple[TreeIndex, TreeIndex]
    ):
        rf_dist = self.rfdist_all_trees[tree_indices]
        return rf_dist[1]


class Experiment:
    def __init__(
        self,
        runs: RunIndexed[Run],
        best_trees: RunIndexed[Newick],
        rfdist_best_trees: RunRunIndexed,
    ):
        self.runs = runs
        self.best_trees = best_trees
        self.rfdist_best_trees = rfdist_best_trees

    def get_best_tree_for_run_index(self, i: RunIndex):
        return self.best_trees[i]

    def get_plain_rfdistance_for_trees(self, run_indices: Tuple[RunIndex, RunIndex]):
        rf_dist = self.rfdist_best_trees[run_indices]
        return rf_dist[0]

    def get_normalized_rfdistance_for_trees(
        self, run_indices: Tuple[RunIndex, RunIndex]
    ):
        rf_dist = self.rfdist_best_trees[run_indices]
        return rf_dist[1]


def read_raxml_all_trees(
    all_raxml_trees_file_path: FilePath,
) -> TreeIndexed[Newick]:
    with open(all_raxml_trees_file_path) as f:
        all_trees = f.readlines()

    return [t.strip() for t in all_trees]


def read_raxml_best_tree(best_raxml_tree_file_path: FilePath) -> Newick:
    with open(best_raxml_tree_file_path) as f:
        # read only first line as second line is an empty line
        return f.readline().strip()


def read_raxml_llhs_all_trees(raxml_treesearch_log: FilePath) -> TreeIndexed[float]:
    return get_all_raxml_llhs(raxml_treesearch_log)


def read_iqtree_all_trees(
    all_iqtree_trees_file_path: FilePath,
) -> TreeIndexed[Newick]:
    with open(all_iqtree_trees_file_path) as f:
        all_trees = f.readlines()

    return [line.strip().split("]")[1] for line in all_trees]


def read_iqtree_results(
    iqtree_results_file_path: FilePath,
) -> TreeIndexed[IqTreeMetrics]:
    return get_iqtree_results(iqtree_results_file_path)


def read_rfdistances_all_trees(
    rfdistances_file_path: FilePath,
) -> TreeTreeIndexed:
    with open(rfdistances_file_path) as f:
        rfdistances = f.readlines()

    res = {}

    for line in rfdistances:
        idx1, idx2, plain, norm = get_cleaned_rf_dist(line)
        res[(idx1, idx2)] = res[(idx2, idx1)] = (plain, norm)

    return res


def read_best_trees_collected(
    best_trees_collected_file_path: FilePath,
) -> RunIndexed[Newick]:
    with open(best_trees_collected_file_path) as f:
        best_trees = f.readlines()

    return [t.strip() for t in best_trees]


def read_rfdistances_best_trees(
    rfdistances_best_file_path: FilePath,
) -> RunRunIndexed:
    with open(rfdistances_best_file_path) as f:
        rfdistances = f.readlines()

    res = {}

    for line in rfdistances:
        runIdx1, runIdx2, plain, norm = get_cleaned_rf_dist(line)
        res[(runIdx1, runIdx2)] = res[(runIdx2, runIdx1)] = (plain, norm)

    return res


def create_Run(
    parameter_file_path: FilePath,
    best_raxml_tree_file_path: FilePath,
    all_raxml_trees_file_path: FilePath,
    raxml_treesearch_log_file_path: FilePath,
    all_iqtree_trees_file_path: FilePath,
    iqtree_results_file_path: FilePath,
    iqtree_test_log_file_path: FilePath,
    rfdistances_file_path: FilePath,
    raxml_rfdistance_logfile_path: FilePath,
) -> Run:

    return Run(
        num_raxml_pars_trees=get_parameter_value(parameter_file_path, "num_pars_trees"),
        num_raxml_rand_trees=get_parameter_value(parameter_file_path, "num_rand_trees"),
        blmin=get_parameter_value(parameter_file_path, "blmin"),
        blmax=get_parameter_value(parameter_file_path, "blmax"),
        raxml_best_llh=get_best_raxml_llh(raxml_treesearch_log_file_path),
        iqtree_best_llh=get_best_iqtree_llh(iqtree_test_log_file_path),
        raxml_best_tree=read_raxml_best_tree(best_raxml_tree_file_path),
        raxml_seeds=get_all_raxml_seeds(raxml_treesearch_log_file_path),
        raxml_all_trees=read_raxml_all_trees(all_raxml_trees_file_path),
        raxml_llhs_all_trees=read_raxml_llhs_all_trees(raxml_treesearch_log_file_path),
        raxml_treesearch_elapsed_time=get_raxml_treesearch_elapsed_time(
            raxml_treesearch_log_file_path
        ),
        iqtree_all_trees=read_iqtree_all_trees(all_iqtree_trees_file_path),
        iqtree_results=read_iqtree_results(iqtree_results_file_path),
        average_absolute_rf_distance=get_raxml_abs_rf_distance(
            raxml_rfdistance_logfile_path
        ),
        average_relative_rf_distance=get_raxml_rel_rf_distance(
            raxml_rfdistance_logfile_path
        ),
        num_unique_topos=get_raxml_num_unique_topos(raxml_rfdistance_logfile_path),
        rfdist_all_trees=read_rfdistances_all_trees(rfdistances_file_path),
    )


def create_Experiment(
    runs: List[Run], best_trees_path: FilePath, rfdist_best_trees_path: FilePath
):
    return Experiment(
        runs,
        read_best_trees_collected(best_trees_path),
        read_rfdistances_best_trees(rfdist_best_trees_path),
    )
