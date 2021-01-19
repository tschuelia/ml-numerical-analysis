
rule topo_rf_distance_per_setting:
    """
    This rule computes the pairwise topological RF distances between found ML trees for one parameter combination.
    """
    input:
        mlTrees = rules.raxml_tree.output.raxml_ml_trees
    output:
        rfDistances = f"{full_file_path_raxml}.raxml.rfDistances",
    params:
        prefix = full_file_path_raxml
    log:
        f"{full_file_path_raxml}.raxml.rfDistances.log",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input.mlTrees} "
        "--prefix {params.prefix} "
        "> {log} "

rule _collect_best_trees:
    """ 
    This rule collects all bestTrees for all parameter combinations in one file.
    """
    input:  
        expand(f"{full_file_path_raxml}.raxml.bestTree", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
    output:
        best_trees_all_runs = f"{outdir}/bestTreesCollected",
    shell:
        "cat {input} > {output} "

rule topo_rf_distance_all_settings:
    """
    This rule computes the pairwise RF distances between all bestTrees.
    """
    input:
        rules._collect_best_trees.output.best_trees_all_runs
    output:
        f"{outdir}/bestTrees.raxml.rfDistances",
    resources:
        mem_mb=5000, # 5GB (in MB)
        runtime=180 # 3 hours (in minutes)
    log:
        f"{outdir}/bestTrees.raxml.rfDistances.log",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {outdir}/bestTrees "
        "> {log} "

