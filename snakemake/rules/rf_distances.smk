
rule topo_distance_per_setting:
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

rule _collect_trees:
    """ 
    This rule collects all mlTrees for all parameter combinations in one file.
    """
    input:  
        expand(f"{full_file_path_raxml}.raxml.mlTrees", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
    output:
        all_trees_all_runs = f"{outdir}/mltrees",
    shell:
        "cat {input} > {output} "

rule topo_rf_distance_all_settings:
    """
    This rule computes the pairwise RF distances between all obtained ML trees for all parameter combinations.
    """
    input:
        rules._collect_trees.output.all_trees_all_runs
    output:
        f"{outdir}/mltrees.raxml.rfDistances",
    log:
        f"{outdir}/mltrees.raxml.rfDistances.log",
    shell:
        "{raxml_command} "
        "--rfdist "
        "--tree {input} "
        "--prefix {outdir}/mltrees "
        "> {log} "

