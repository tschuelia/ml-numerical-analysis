rule collect_all_iqtree_trees_per_combination:
    """
    Rule to collect all parsimony tree search trees for one parameter combination in a single file.
    """
    input:
        pars_trees = expand(f"{full_file_path_iqtree_pars}.treefile", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_trees = f"{full_file_path_iqtree}.allTreesCollected"
    shell:
        "cat {input.pars_trees} > {output.all_iqtree_trees}"


rule collect_all_iqtree_logs_per_combination:
    """
    Rule to collect all iqtree log files for one parameter combination in a single file.
    """
    input:
        pars_trees_logs = expand(f"{full_file_path_iqtree_pars}.iqtree.treesearch.log", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_logs = f"{full_file_path_iqtree}.allTreesearchLogs"
    shell:
        "cat {input.pars_trees_logs} > {output.all_iqtree_logs}"


rule collect_all_iqtree_eval_trees_per_combination:
    """
    Rule to collect all eval trees for one parameter combination in a single file.
    """
    input:
        eval_trees = expand(f"{full_file_path_iqtree_eval}.treefile", seed=pars_seeds, allow_missing=True)
    output:
        all_eval_trees = f"{full_file_path_iqtree}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"


rule collect_all_iqtree_eval_logs_per_combination:
    """
    Rule to collect all eval logs for one parameter combination in a single file.
    """
    input:
        eval_logs = expand(f"{full_file_path_iqtree_eval}.iqtree.eval.log", seed=pars_seeds, allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_iqtree}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"


rule collect_all_iqtree_eval_trees:
    input:
        trees = expand(f"{full_file_path_iqtree}.allEvalTreesCollected", lh_eps=lh_eps_opts),
    output:
        all_trees=f"{base_dir_iqtree}allEvalTreesCollected",
    shell:
        "cat {input.trees} > {output.all_trees}"


rule collect_best_overall_iqtree_eval_tree:
    """
    Rule to find the best overall eval tree. 
    The best eval tree is the one with the highest llh score with lowest runtime.
    """
    input:
        trees = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun", lh_eps=lh_eps_opts),
        logs = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun.json", lh_eps=lh_eps_opts),
    output:
        best_overall_tree = f"{base_dir_iqtree}bestOverallEvalTree",
    script:
        "scripts/save_best_overall_tree.py"