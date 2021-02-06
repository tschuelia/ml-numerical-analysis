rule collect_all_iqtree_trees:
    input:
        pars_trees = expand(f"{full_file_path_iqtree_pars}.treefile", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_trees = f"{full_file_path_iqtree}.allTreesCollected"
    shell:
        "cat {input.pars_trees} > {output.all_iqtree_trees}"
    
rule collect_all_iqtree_logs:
    input:
        pars_trees_logs = expand(f"{full_file_path_iqtree_pars}.iqtree.treesearch.log", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_logs = f"{full_file_path_iqtree}.allTreesearchLogs"
    shell:
        "cat {input.pars_trees_logs} > {output.all_iqtree_logs}"


rule collect_all_iqtree_eval_trees:
    input:
        eval_trees = expand(f"{full_file_path_iqtree_eval}.treefile", blmin_eval=blmin_opts, blmax_eval=blmax_opts, lh_eps_eval=lh_eps_opts, allow_missing=True)
    output:
        all_eval_trees = f"{full_file_path_iqtree}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"

rule collect_all_iqtree_eval_logs:
    input:
        eval_logs = expand(f"{full_file_path_iqtree_eval}.iqtree.eval.log", blmin_eval=blmin_opts, blmax_eval=blmax_opts, lh_eps_eval=lh_eps_opts, allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_iqtree}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"