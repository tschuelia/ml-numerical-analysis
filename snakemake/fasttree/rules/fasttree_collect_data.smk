rule collect_all_fasttree_trees:
    input:
        trees = expand(f"{full_file_path_fasttree_tree}.treefile", seed=seeds, allow_missing=True),
    output:
        all_fasttree_trees=f"{full_file_path_fasttree}.allTreesCollected"
    shell:
        "cat {input.trees} > {output.all_fasttree_trees}"

rule collect_all_fasttree_logs:
    input:
        trees_logs = expand(f"{full_file_path_fasttree_tree}.fasttree.treesearch.log", seed=seeds,
            allow_missing=True),
    output:
        all_fasttree_logs=f"{full_file_path_fasttree}.allTreesearchLogs"
    shell:
        "cat {input.trees_logs} > {output.all_fasttree_logs}"

