rule iqtree_statstest_on_fasttree_trees:
    input:
        msa         = config["data"]["input"],
        all_trees   = rules.collect_best_fasttree_trees.output.best_trees_all_runs,
        best_tree   = rules.collect_best_overall_fasttree_tree.output.best_overall_tree,

    output:
        summary     = f"{base_dir_fasttree}significance.iqtree",
        iqtree_log  = f"{base_dir_fasttree}significance.iqtree.log",

    params:
        model           = config["parameters"]["model"]["iqtree"],
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix          = f"{base_dir_fasttree}significance",
    log:
        f"{base_dir_fasttree}significance.iqtree.snakelog",
    shell:
        "{iqtree_command} "
        "-s {input.msa} "
        "-m {params.model} "
        "-pre {params.prefix} "
        "-z {input.all_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "> {output.iqtree_log} "