rule save_best_fasttree_tree_to_file:
    input:
        trees = rules.collect_all_fasttree_trees.output.all_fasttree_trees,
        logs = rules.collect_all_fasttree_logs.output.all_fasttree_logs,
    output:
        best_tree = f"{full_file_path_fasttree}.bestTreeOfRun",
        best_log = f"{full_file_path_fasttree}.bestTreeOfRun.json",
    script:
        "scripts/save_best_tree.py"