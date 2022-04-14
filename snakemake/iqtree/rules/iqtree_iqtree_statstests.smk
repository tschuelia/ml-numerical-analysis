model = config["parameters"]["model"]["iqtree"]
partitioned = "/" in model

"""
IQ-Tree significance tests for only the eval trees
"""
rule iqtree_filter_unique_tree_topologies:
    input:
        all_trees = rules.collect_all_iqtree_eval_trees.output.all_trees,
    output:
        raxml_rfdist_log = f"{base_dir_iqtree}filteredEvalTrees.raxml.log",
        filtered_trees  = f"{base_dir_iqtree}filteredEvalTrees",
        clusters        = f"{base_dir_iqtree}filteredEvalTreesClusters",
    params:
        raxml_command = raxml_command
    script:
        "scripts/filter_tree_topologies.py"


rule iqtree_statstest_on_iqtree_eval_trees:
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.iqtree_filter_unique_tree_topologies.output.filtered_trees,
        best_tree       = rules.collect_best_overall_iqtree_eval_tree.output.best_overall_tree,

    output:
        summary     = f"{base_dir_iqtree}significance.iqtree",
        iqtree_log  = f"{base_dir_iqtree}significance.iqtree.log",

    params:
        model           = model,
        model_str       = "-p" if partitioned else "-m",
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix          = f"{base_dir_iqtree}significance",
    log:
        f"{base_dir_iqtree}significance.iqtree.snakelog",
    shell:
        "{iqtree_command} "
        "-s {input.msa} "
        "{params.model_str} {params.model} "
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


"""
Since the IQ-Tree test results change with the order of the input trees we compare
each tree to the best tree in a pairwise comparison instead 
"""
max_workers = 2

rule pairwise_iqtree_statstest_on_iqtree_eval_trees:
    """
    In order to save time, we use the filtered trees from above
    Trees with the same topology should yield the same test results so we can map this back 
    later on
    """
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.iqtree_filter_unique_tree_topologies.output.filtered_trees,
        best_tree       = rules.collect_best_overall_iqtree_eval_tree.output.best_overall_tree,
    output:
        iqtree_results = f"{base_dir_iqtree}pairwiseSignificanceEval.json",
    params:
        iqtree_command  = iqtree_command,
        model           = model,
        partitioned     = partitioned,
        max_workers     = max_workers,
    script:
        "scripts/pairwise_iqtree.py"
