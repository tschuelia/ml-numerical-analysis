rule collect_all_raxml_trees:
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

rule collect_all_raxml_logs:
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

rule collect_all_raxml_eval_trees:
    """
    Rule to collect all eval trees for one parameter combination in a single file.
    """
    input:
        eval_trees = expand(f"{full_file_path_raxml_eval}.raxml.bestTree", 
                                blmin_eval=blmin_eval_opts, blmax_eval=blmax_eval_opts, lh_eps_eval=lh_eps_eval_opts, raxml_param_epsilon_eval=raxml_param_epsilon_eval_opts, raxml_brlen_smoothings_eval=raxml_brlen_smoothings_eval_opts,
                                allow_missing=True)
    output:
        all_eval_trees  = f"{full_file_path_raxml}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"

rule collect_all_raxml_eval_logs:
    """
    Rule to collect all eval logs for one parameter combination in a single file.
    """
    input:
        eval_logs = expand(f"{full_file_path_raxml_eval}.raxml.eval.log", 
                                blmin_eval=blmin_eval_opts, blmax_eval=blmax_eval_opts, lh_eps_eval=lh_eps_eval_opts, raxml_param_epsilon_eval=raxml_param_epsilon_eval_opts, raxml_brlen_smoothings_eval=raxml_brlen_smoothings_eval_opts,
                                allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_raxml}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"

rule collect_best_trees:
    """
    Rule to collect the best tree search tree for each parameter combination and write them all into a single file.
    """
    input:  
        expand(f"{full_file_path_raxml}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
    output:
        best_trees_all_runs = f"{base_dir_raxml}bestTreesCollected",
    shell:
        "cat {input} > {output} "
    
rule collect_best_eval_trees:
    """
    Rule to collect the best eval tree for each parameter combination and write them all into a single file.
    """
    input:  
        expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
    output:
        best_trees_all_runs = f"{base_dir_raxml}bestEvalTreesCollected",
    shell:
        "cat {input} > {output} "

rule collect_best_overall_eval_tree:
    """
    Rule to find the best overall eval tree. 
    The best eval tree is the one with the highest llh score with lowest runtime.
    """
    input:
        trees = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
        logs = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun.json", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
    output:
        best_overall_tree = f"{base_dir_raxml}bestOverallEvalTree",
    script:
        "scripts/save_best_overall_tree.py"