rule fasttree_filter_unique_tree_topologies:
     input:
         all_trees = rules.collect_all_fasttree_trees.output.all_trees,
     output:
         filtered_trees  = f"{base_dir_fasttree}filteredTrees",
         clusters        = f"{base_dir_fasttree}filteredTreesClusters",
     params:
         raxml_command = raxml_command
     script:
         "scripts/filter_tree_topologies.py"

rule iqtree_statstest_on_fasttree_trees:
    input:
        msa                 = config["data"]["input"],
        filtered_trees      = rules.fasttree_filter_unique_tree_topologies.output.filtered_trees,
        best_tree           = rules.collect_best_overall_fasttree_tree.output.best_overall_tree,

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
        "-z {input.filtered_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "-seed 0 "
        "> {output.iqtree_log} "