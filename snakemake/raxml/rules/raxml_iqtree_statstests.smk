"""
IQ-Tree significance tests for only the eval trees
"""
rule raxml_filter_unique_tree_topologies:
    input:
        all_trees = rules.collect_all_raxml_eval_trees.output.all_trees,
    output:
        raxml_rfdist_log    = f"{base_dir_raxml}filteredEvalTrees.raxml.log",
        filtered_trees      = f"{base_dir_raxml}filteredEvalTrees",
        clusters            = f"{base_dir_raxml}filteredEvalTreesClusters",
    params:
        raxml_command = raxml_command
    script:
        "scripts/filter_tree_topologies.py"


rule iqtree_statstest_on_raxml_eval_trees:
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.raxml_filter_unique_tree_topologies.output.filtered_trees,
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


"""
IQ-Tree significance tests the eval trees and the re-evaluated starting trees
"""
rule raxml_filter_unique_tree_topologies_eval_and_starting:
    input:
        all_trees = rules.collect_raxml_starting_eval_and_eval_trees.output.all_trees,
    output:
        raxml_rfdist_log    = f"{base_dir_raxml}filteredEvalAndStartingTrees.raxml.log",
        filtered_trees      = f"{base_dir_raxml}filteredEvalAndStartingTrees",
        clusters            = f"{base_dir_raxml}filteredEvalAndStartingTreesClusters",
    params:
        raxml_command = raxml_command
    script:
        "scripts/filter_tree_topologies.py"


rule iqtree_statstest_on_raxml_eval_and_starting_trees:
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.raxml_filter_unique_tree_topologies_eval_and_starting.output.filtered_trees,
        best_tree       = rules.collect_best_overall_eval_tree.output.best_overall_tree,

    output:
        summary     = f"{base_dir_raxml}significanceEvalAndStarting.iqtree",
        iqtree_log  = f"{base_dir_raxml}significanceEvalAndStarting.iqtree.log",

    params:
        model           = config["parameters"]["model"]["iqtree"],
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix          = f"{base_dir_raxml}significanceEvalAndStarting",
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


"""
Since the IQ-Tree test results change with the order of the input trees we compare
each tree to the best tree in a pairwise comparison instead 
"""
max_workers = 2

rule pairwise_iqtree_statstest_on_raxml_eval_trees:
    """
    In order to save time, we use the filtered trees from above
    Trees with the same topology should yield the same test results so we can map this back 
    later on
    """
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.raxml_filter_unique_tree_topologies.output.filtered_trees,
        best_tree       = rules.collect_best_overall_eval_tree.output.best_overall_tree,
    output:
        iqtree_results = f"{base_dir_raxml}pairwiseSignificanceEval.json",
    params:
        iqtree_command  = iqtree_command,
        model           = config["parameters"]["model"]["iqtree"],
        max_workers     = max_workers,
    script:
        "scripts/pairwise_iqtree.py"


rule pairwise_iqtree_statstest_on_raxml_eval_and_starting_trees:
    """
    In order to save time, we use the filtered trees from above
    Trees with the same topology should yield the same test results so we can map this back 
    later on
    """
    input:
        msa             = config["data"]["input"],
        filtered_trees  = rules.raxml_filter_unique_tree_topologies_eval_and_starting.output.filtered_trees,
        best_tree       = rules.collect_best_overall_eval_tree.output.best_overall_tree,
    output:
        iqtree_results = f"{base_dir_raxml}pairwiseSignificanceEvalAndStarting.json",
    params:
        iqtree_command  = iqtree_command,
        model           = config["parameters"]["model"]["iqtree"],
        max_workers     = max_workers,
    script:
        "scripts/pairwise_iqtree.py"