import dataclasses

from snakelib.custom_types import *
from snakelib.utils import (
    NewickTree,
    get_parameter_value,
    read_file_contents,
    parse_newick_string
)

from snakelib.program_parser import Program

from fasttree_utils import (
    get_all_fasttree_llhs,
    get_best_fasttree_llh,
    get_fasttree_treesearch_entire_run,
    get_fasttree_runtimes,
    get_fasttree_run_param_values_from_file,
)



def create_fasttree(
        parameter_file_path: FilePath,
        treesearch_log_file_path: FilePath,
        best_tree_file_path: FilePath,
        all_treesearch_trees_file_path: FilePath,
):
    treesearch_file_content = read_file_contents(all_treesearch_trees_file_path)
    treesearch_trees = [parse_newick_string(newick_str) for newick_str in treesearch_file_content]
    # fmt: off
    fasttree = Program(
        blmin                   = get_parameter_value(parameter_file_path, "blmin"),
        blmax                   = None,
        lh_epsilon              = get_parameter_value(parameter_file_path, "lh_eps"),
        model_epsilon     = None,
        branch_length_smoothing = None,
        spr_lh_epsilon          = None,
        bfgs_factor             = None,

        num_pars_trees          = int(get_parameter_value(parameter_file_path, "num_pars_trees")),
        num_rand_trees          = 0, # for fasttree we cannot chose the starting tree type
        best_treesearch_llh     = get_best_fasttree_llh(treesearch_log_file_path),
        best_evaluation_llh     = None,
        treesearch_total_time   = get_fasttree_treesearch_entire_run(treesearch_log_file_path),

        # Tree search
        best_treesearch_tree        = parse_newick_string(read_file_contents(best_tree_file_path)[0]),
        treeseach_seeds             = get_fasttree_run_param_values_from_file(treesearch_log_file_path, "seed"),
        treesearch_trees            = treesearch_trees,
        treesearch_llhs             = get_all_fasttree_llhs(treesearch_log_file_path),
        treesearch_compute_times    = get_fasttree_runtimes(treesearch_log_file_path),
        starting_tree_types         = None,

        # Eval
        best_eval_tree              = None,
        eval_blmins                 = None,
        eval_blmaxs                 = None,
        eval_lh_epsilons            = None,
        eval_model_param_epsilons   = None,
        eval_raxml_brlen_smoothings = None,
        eval_bfgs_factors           = None,

        eval_trees          = None,
        eval_llhs           = None,
        eval_compute_times  = None,

    )
    # fmt: on

    return fasttree
