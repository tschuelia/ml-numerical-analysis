rule collect_all_raxml_trees_per_combination:
    """
    Rule to collect all parsimony and random tree search trees for one parameter combination in a single file.
    """
    input:
        pars_trees      = expand(f"{full_file_path_raxml_pars}.raxml.bestTree", seed=pars_seeds, allow_missing=True),
        rand_trees      = expand(f"{full_file_path_raxml_rand}.raxml.bestTree", seed=rand_seeds, allow_missing=True),
    output:
        all_raxml_trees = f"{full_file_path_raxml}.allTreesCollected",
    shell:
        "cat {input.pars_trees} {input.rand_trees} > {output.all_raxml_trees}"


rule collect_all_raxml_logs_per_combination:
    """
    Rule to collect all raxml log files for one parameter combination in a single file.
    """
    input:
        pars_trees_logs = expand(f"{full_file_path_raxml_pars}.raxml.treesearch.log", seed=pars_seeds, allow_missing=True),
        rand_trees_logs = expand(f"{full_file_path_raxml_rand}.raxml.treesearch.log", seed=rand_seeds, allow_missing=True),
    output:
        all_raxml_logs  = f"{full_file_path_raxml}.allTreesearchLogs",
    shell:
        "cat {input.pars_trees_logs} {input.rand_trees_logs} > {output.all_raxml_logs}"


rule collect_all_raxml_eval_trees_per_combination:
    """
    Rule to collect all eval trees for one parameter combination in a single file.
    """
    input:
        eval_trees_pars = expand(f"{full_file_path_raxml_eval_pars}.raxml.bestTree", seed=pars_seeds, allow_missing=True),
        eval_trees_rand = expand(f"{full_file_path_raxml_eval_rand}.raxml.bestTree", seed=rand_seeds, allow_missing=True)
    output:
        all_eval_trees  = f"{full_file_path_raxml}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees_pars} {input.eval_trees_rand} > {output.all_eval_trees}"


rule collect_all_raxml_eval_logs_per_combination:
    """
    Rule to collect all eval logs for one parameter combination in a single file.
    """
    input:
        eval_logs_pars = expand(f"{full_file_path_raxml_eval_pars}.raxml.eval.log", seed=pars_seeds, allow_missing=True),
        eval_logs_rand = expand(f"{full_file_path_raxml_eval_rand}.raxml.eval.log", seed=rand_seeds, allow_missing=True),
    output:
        all_eval_logs = f"{full_file_path_raxml}.allEvalLogs"
    shell:
        "cat {input.eval_logs_pars} {input.eval_logs_rand} > {output.all_eval_logs}"


rule collect_all_raxml_eval_trees:
    input:
        trees = expand(f"{full_file_path_raxml}.allEvalTreesCollected", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
    output:
        all_trees=f"{base_dir_raxml}allEvalTreesCollected",
    script:
        "scripts/cat.py"


rule collect_best_overall_eval_tree:
    """
    Rule to find the best overall eval tree. 
    The best eval tree is the one with the highest llh score with lowest runtime.
    """
    input:
        trees = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        logs = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun.json", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
    output:
        best_overall_tree = f"{base_dir_raxml}bestOverallEvalTree",
    script:
        "scripts/save_best_overall_tree.py"