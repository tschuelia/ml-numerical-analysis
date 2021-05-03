rule collect_all_iqtree_trees:
    """
    Rule to collect all parsimony tree search trees for one parameter combination in a single file.
    """
    input:
        pars_trees = expand(f"{full_file_path_iqtree_pars}.treefile", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_trees = f"{full_file_path_iqtree}.allTreesCollected"
    shell:
        "cat {input.pars_trees} > {output.all_iqtree_trees}"
    
rule collect_all_iqtree_logs:
    """
    Rule to collect all iqtree log files for one parameter combination in a single file.
    """
    input:
        pars_trees_logs = expand(f"{full_file_path_iqtree_pars}.iqtree.treesearch.log", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_logs = f"{full_file_path_iqtree}.allTreesearchLogs"
    shell:
        "cat {input.pars_trees_logs} > {output.all_iqtree_logs}"


rule collect_all_iqtree_eval_trees:
    """
    Rule to collect all eval trees for one parameter combination in a single file.
    """
    input:
        eval_trees = expand(f"{full_file_path_iqtree_eval}.treefile", blmin_eval=blmin_eval_opts, blmax_eval=blmax_eval_opts, lh_eps_eval=lh_eps_eval_opts, allow_missing=True)
    output:
        all_eval_trees = f"{full_file_path_iqtree}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"

rule collect_all_iqtree_eval_logs:
    """
    Rule to collect all eval logs for one parameter combination in a single file.
    """
    input:
        eval_logs = expand(f"{full_file_path_iqtree_eval}.iqtree.eval.log", blmin_eval=blmin_eval_opts, blmax_eval=blmax_eval_opts, lh_eps_eval=lh_eps_eval_opts, allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_iqtree}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"

rule collect_best_iqtree_eval_trees:
    """
    Rule to collect the best eval tree for each parameter combination and write them all into a single file.
    """
    input:
        expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts),
    output:
        best_trees_all_runs = f"{base_dir_iqtree}bestEvalTreesCollected",
    script:
        "scripts/cat_trees.py"

rule collect_best_overall_iqtree_eval_tree:
    """
    Rule to find the best overall eval tree. 
    The best eval tree is the one with the highest llh score with lowest runtime.
    """
    input:
        trees = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts),
        logs = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun.json", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts),
    output:
        best_overall_tree = f"{base_dir_iqtree}bestOverallEvalTree",
    script:
        "scripts/save_best_overall_tree.py"