import dataclasses

from typing import Dict, List, Tuple

from custom_types import *
from iqtree_statstest_parser import get_iqtree_results
from utils import (
    get_cleaned_rf_dist,
    get_parameter_value,
    read_file_contents,
    read_rfdistances
)

@dataclasses.dataclass
class Run:
    blmin: float
    blmax: float

def create_Run(parameter_file_path: FilePath):
    return Run(
        blmin = get_parameter_value(parameter_file_path, "blmin"),
        blmax = get_parameter_value(parameter_file_path, "blmax"),
    )


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
    #runs                        : List[Run],
    raxml_best_trees_path       : FilePath,
    raxml_best_eval_trees_path  : FilePath,
    rfdist_raxml_best_trees_path      : FilePath,
    rfdist_raxml_best_eval_trees_path : FilePath,
):
    return Experiment(
        #runs                          = runs,
        raxml_best_trees              = read_file_contents(raxml_best_trees_path),
        raxml_best_eval_trees         = read_file_contents(raxml_best_eval_trees_path),
        rfdist_raxml_best_trees       = read_rfdistances(rfdist_raxml_best_trees_path),
        rfdist_raxml_best_eval_trees  = read_rfdistances(rfdist_raxml_best_eval_trees_path),
    )
# fmt: on