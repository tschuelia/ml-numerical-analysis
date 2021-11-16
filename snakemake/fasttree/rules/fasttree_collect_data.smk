rule collect_all_fasttree_trees_per_combination:
    input:
        trees = expand(f"{full_file_path_fasttree_tree}.treefile", seed=seeds, allow_missing=True),
    output:
        all_fasttree_trees=f"{full_file_path_fasttree}.allTreesCollected"
    shell:
        "cat {input.trees} > {output.all_fasttree_trees}"

rule collect_all_fasttree_logs_per_combination:
    input:
        trees_logs = expand(f"{full_file_path_fasttree_tree}.fasttree.treesearch.log", seed=seeds,
            allow_missing=True),
    output:
        all_fasttree_logs=f"{full_file_path_fasttree}.allTreesearchLogs"
    shell:
        "cat {input.trees_logs} > {output.all_fasttree_logs}"

rule collect_best_fasttree_trees:
    input:
        trees = expand(f"{full_file_path_fasttree}.bestTreeOfRun", blmin=blmin_opts, lh_eps=lh_eps_opts),
    output:
        best_trees_all_runs = f"{base_dir_fasttree}bestTreesCollected"
    shell:
        "cat {input} > {output} "

rule collect_best_overall_fasttree_tree:
    input:
        trees = expand(f"{full_file_path_fasttree}.bestTreeOfRun", blmin=blmin_opts, lh_eps=lh_eps_opts),
        logs = expand(f"{full_file_path_fasttree}.bestTreeOfRun.json", blmin=blmin_opts, lh_eps=lh_eps_opts),
    output:
        best_overall_tree = f"{base_dir_fasttree}bestOverallTree",
    script:
        "scripts/save_best_overall_tree.py"

rule collect_all_fasttree_trees:
    input:
        trees = expand(f"{full_file_path_fasttree_tree}.treefile",seed=seeds, blmin=blmin_opts, lh_eps=lh_eps_opts)
    output:
        all_trees = f"{base_dir_fasttree}allTreesCollected"
    shell:
        "cat {input.trees} > {output.all_trees}"