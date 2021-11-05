rule filter_unique_tree_topologies:
    input:
        all_trees = rules.collect_best_eval_trees.output.best_trees_all_runs,
    output:
        filtered_trees  = f"{base_dir_raxml}filteredEvalTrees",
        clusters        = f"{base_dir_raxml}filteredEvalTreesClusters",
    params:
        raxml_command = raxml_command
    script:
        "scripts/filter_tree_topologies.py"


rule iqtree_statstest_on_raxml_eval_trees:
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.filter_unique_tree_topologies.output.filtered_trees,
        best_tree       = rules.collect_best_overall_eval_tree.output.best_overall_tree,

    output:
        summary     = f"{base_dir_raxml}significance.iqtree",
        iqtree_log  = f"{base_dir_raxml}significance.iqtree.log",

    params:
        model           = config["parameters"]["model"]["iqtree"],
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix          = f"{base_dir_raxml}significance",
    log:
        f"{base_dir_raxml}significance.iqtree.snakelog",
    shell:
        "{iqtree_command} "
        "-s {input.msa} "
        "-m {params.model} "
        "-pre {params.prefix} "
        "-z {input.filtered_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "-seed 0 "
        "> {output.iqtree_log} "