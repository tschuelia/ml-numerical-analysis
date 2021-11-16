import dataclasses

from snakelib.custom_types import *
from snakelib.utils import (
    NewickTree,
    get_parameter_value,
    read_file_contents,
    parse_newick_string
)

from snakelib.program_parser import Program

from raxml_utils import (
    get_all_raxml_llhs,
    get_all_raxml_seeds,
    get_best_raxml_llh,
    get_raxml_run_param_values_from_file,
    get_raxml_elapsed_time,
    get_raxml_treesearch_elapsed_time_entire_run,
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
        blmin_eval: float = None,
        blmax_eval: float = None,
        lh_eps_eval: float = None,
        model_param_epsilon_eval: float = None,
        raxml_brlen_smoothings_eval: float = None,
        bfgs_fac_eval: float = None,
):
    treesearch_trees_file_content = read_file_contents(all_treesearch_trees_file_path)
    treesearch_trees = [parse_newick_string(newick_str) for newick_str in treesearch_trees_file_content]

    eval_trees_file_content = read_file_contents(all_eval_trees_file_path)
    eval_trees = [parse_newick_string(newick_str) for newick_str in eval_trees_file_content]

    eval_blmins = [blmin_eval] * len(eval_trees) if blmin_eval else get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "blmin")
    eval_blmaxs = [blmax_eval] * len(eval_trees) if blmax_eval else get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "blmax")
    eval_lh_epsilons = [lh_eps_eval] * len(eval_trees) if lh_eps_eval else get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "lh-epsilon")
    eval_model_param_epsilons = [model_param_epsilon_eval] * len(eval_trees) if model_param_epsilon_eval else get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "param-eps")
    eval_raxml_brlen_smoothings = [raxml_brlen_smoothings_eval] * len(eval_trees) if raxml_brlen_smoothings_eval else get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "brlen-smoothings")
    eval_bfgs_factors = [bfgs_fac_eval] * len(eval_trees) if bfgs_fac_eval else get_raxml_run_param_values_from_file(eval_log_file_path, raxml_command, "bfgs-factor")

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
        eval_blmins                 = eval_blmins,
        eval_blmaxs                 = eval_blmaxs,
        eval_lh_epsilons            = eval_lh_epsilons,
        eval_model_param_epsilons   = eval_model_param_epsilons,
        eval_raxml_brlen_smoothings = eval_raxml_brlen_smoothings,
        eval_bfgs_factors           = eval_bfgs_factors,

        eval_trees          = eval_trees,
        eval_llhs           = get_all_raxml_llhs(eval_log_file_path),
        eval_compute_times  = get_raxml_elapsed_time(eval_log_file_path),
    )
    return raxml
# fmt: on
