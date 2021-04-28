rule save_best_iqtree_tree_to_file:
    input:
        trees = rules.collect_all_iqtree_trees.output.all_iqtree_trees,
        logs = rules.collect_all_iqtree_logs.output.all_iqtree_logs,
    output:
        best_tree = f"{full_file_path_iqtree}.bestTreeOfRun",
        best_log = f"{full_file_path_iqtree}.bestTreeOfRun.json",
    script:
        "scripts/save_best_tree.py"

rule save_best_iqtree_eval_tree_to_file:
    input:
        trees = rules.collect_all_iqtree_eval_trees.output.all_eval_trees,
        logs = rules.collect_all_iqtree_eval_logs.output.all_eval_logs,
    output:
        best_tree = f"{full_file_path_iqtree}.bestEvalTreeOfRun",
        best_log = f"{full_file_path_iqtree}.bestEvalTreeOfRun.json",
    script:
        "scripts/save_best_tree.py"