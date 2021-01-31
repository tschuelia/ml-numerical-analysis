rule topo_rf_distance_per_setting:
    """
    This rule computes the pairwise topological RF distances between all found ML trees for one parameter combination.
    """
    input:
        rules.collect_all_raxml_trees.output.all_raxml_trees
    output:
        rfDistances = f"{full_file_path_raxml}.raxml.rfDistances",
    params:
        prefix = full_file_path_raxml
    log:
        f"{full_file_path_raxml}.raxml.rfDistances.log",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {params.prefix} "
        "> {log} "


rule collect_best_trees:
    input:  
        expand(f"{full_file_path_raxml}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
    output:
        best_trees_all_runs = f"{outdir}/bestTreesCollected",
    shell:
        "cat {input} > {output} "

rule topo_rf_distance_all_settings_with_best_of_run:
    input:
        rules.collect_best_trees.output.best_trees_all_runs
    output:
        f"{outdir}/bestTrees.raxml.rfDistances",
    log:
        f"{outdir}/bestTrees.raxml.rfDistances.log",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {outdir}/bestTrees "
        "> {log} "

rule collect_best_eval_trees:
    input:  
        expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
    output:
        best_trees_all_runs = f"{outdir}/bestEvalTreesCollected",
    shell:
        "cat {input} > {output} "

rule topo_rf_distance_all_settings_with_best_of_eval:
    input:
        rules.collect_best_eval_trees.output.best_trees_all_runs
    output:
        f"{outdir}/bestEvalTrees.raxml.rfDistances",
    log:
        f"{outdir}/bestEvalTrees.raxml.rfDistances.log",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {outdir}/bestEvalTrees "
        "> {log} "