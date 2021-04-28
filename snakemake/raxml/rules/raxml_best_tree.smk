rule save_best_raxml_tree_to_file:
    input:
        trees = rules.collect_all_raxml_trees.output.all_raxml_trees,
        logs = rules.collect_all_raxml_logs.output.all_raxml_logs,
    output:
        best_tree = f"{full_file_path_raxml}.bestTreeOfRun",
        best_log = f"{full_file_path_raxml}.bestTreeOfRun.json",
    script:
        "scripts/save_best_tree.py"

rule save_best_raxml_eval_tree_and_log_to_file:
    input:
        trees = rules.collect_all_raxml_eval_trees.output.all_eval_trees,
        logs = rules.collect_all_raxml_eval_logs.output.all_eval_logs,
    output:
        best_tree = f"{full_file_path_raxml}.bestEvalTreeOfRun",
        best_log = f"{full_file_path_raxml}.bestEvalTreeOfRun.json",
    script:
        "scripts/save_best_tree.py"