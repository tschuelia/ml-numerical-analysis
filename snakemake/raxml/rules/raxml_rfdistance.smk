rule rf_distance_treesearch_trees:
    input:
        trees = rules.collect_all_raxml_trees.output.all_raxml_trees
    output:
        rfDist      = f"{full_file_path_raxml}.raxml.rfDistances",
        rfDist_log  = f"{full_file_path_raxml}.raxml.rfDistances.log",
    params:
        prefix = full_file_path_raxml
    log:
        f"{full_file_path_raxml}.raxml.rfDistances.snakelog",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input.trees} "
        "--prefix {params.prefix} "
        ">> {output.rfDist_log} "

rule rf_distance_best_treesearch_trees:
    input:
        rules.collect_best_trees.output.best_trees_all_runs
    output:
        rfDist      = f"{base_dir_raxml}bestTrees.raxml.rfDistances",
        rfDist_log  = f"{base_dir_raxml}bestTrees.raxml.rfDistances.log",
    log:
        f"{base_dir_raxml}bestTrees.raxml.rfDistances.snakelog",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {base_dir_raxml}/bestTrees "
        ">> {output.rfDist_log} "

rule topo_rf_distance_all_settings_with_best_of_eval:
    input:
        rules.collect_best_eval_trees.output.best_trees_all_runs
    output:
        rfDist      = f"{base_dir_raxml}bestEvalTrees.raxml.rfDistances",
        rfDist_log  = f"{base_dir_raxml}bestEvalTrees.raxml.rfDistances.log",
    log:
        f"{base_dir_raxml}bestEvalTrees.raxml.rfDistances.snakelog",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {base_dir_raxml}/bestEvalTrees "
        ">> {output.rfDist_log} "