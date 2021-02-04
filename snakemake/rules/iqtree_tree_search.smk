rule iqtree_pars_tree:
    input:
        msa = config["data"]["input"],
    output:
        treefile    = f"{full_file_path_iqtree_pars}.treefile",
        iqtree_done = f"{full_file_path_iqtree_pars}.done",
        iqtree_log = f"{full_file_path_iqtree_pars}.iqtree.treesearch.log"
    params:
        model       = config["parameters"]["model"]["iqtree"],
        threads     = config["parameters"]["iqtree"]["threads"],
        prefix      = full_file_path_iqtree_pars,
    log:
        f"{full_file_path_iqtree_pars}.snakelog"
    shell:
        # prevent snakemake from rerunning finished iqtree runs
        "grep -q 'BEST SCORE FOUND' {params.prefix}.iqtree.treesearch.log || "
        "{iqtree_command} "
        "-m {params.model} "
        "-s {input.msa} "
        "-ninit 1 "
        "-blmin {wildcards.blmin} "
        "-blmax {wildcards.blmax} "
        "-seed {wildcards.seed} "
        "-pre {params.prefix} "
        "-nt {params.threads} "
        "> {output.iqtree_log} "
        "&& touch {params.prefix}.done"

rule collect_all_iqtree_trees:
    input:
        pars_trees = expand(f"{full_file_path_iqtree_pars}.treefile", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_trees = f"{full_file_path_iqtree}.allTreesCollected"
    shell:
        "cat {input.pars_trees} > {output.all_iqtree_trees}"
    
rule collect_all_iqtree_logs:
    input:
        pars_trees_logs = expand(f"{full_file_path_iqtree_pars}.iqtree.treesearch.log", seed=pars_seeds, allow_missing=True),
    output:
        all_iqtree_logs = f"{full_file_path_iqtree}.allTreesearchLogs"
    shell:
        "cat {input.pars_trees_logs} > {output.all_iqtree_logs}"

rule save_best_iqtree_tree_to_file:
    input:
        trees = rules.collect_all_iqtree_trees.output.all_iqtree_trees,
        logs = rules.collect_all_iqtree_logs.output.all_iqtree_logs,
    output:
        best_tree = f"{full_file_path_iqtree}.bestTreeOfRun"
    params:
        program = "iqtree"
    script:
        "scripts/save_best_tree.py"

rule re_eval_best_iqtree_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_iqtree}.bestTreeOfRun",
    
    output:
        treefile    = f"{full_file_path_iqtree_eval}.treefile",
        iqtree_done = f"{full_file_path_iqtree_eval}.done",
        eval_log    = f"{full_file_path_iqtree_eval}.iqtree.eval.log"
    params:
        model           = config["parameters"]["model"]["iqtree"],
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix          = full_file_path_iqtree_eval,
    log:
        f"{full_file_path_iqtree_eval}.snakelog"
    shell:
        # prevent snakemake from rerunning finished iqtree runs
        "grep -q 'BEST SCORE FOUND' {params.prefix}.iqtree.eval.log || "
        "{iqtree_command} "
        "-m {params.model} "
        "-s {input.msa} "
        "-te {input.best_tree_of_run} "
        "-blmin {wildcards.blmin_eval} "
        "-blmax {wildcards.blmax_eval} "
        "-pre {params.prefix} "
        "-nt {params.threads} "
        "> {output.eval_log} "
        "&& touch {params.prefix}.done"

rule collect_all_iqtree_eval_trees:
    input:
        eval_trees = expand(f"{full_file_path_iqtree_eval}.treefile", blmin_eval=blmin_opts, blmax_eval=blmax_opts, allow_missing=True)
    output:
        all_eval_trees = f"{full_file_path_iqtree}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"

rule collect_all_iqtree_eval_logs:
    input:
        eval_logs = expand(f"{full_file_path_iqtree_eval}.iqtree.eval.log", blmin_eval=blmin_opts, blmax_eval=blmax_opts, allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_iqtree}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"

rule save_best_iqtree_eval_tree_to_file:
    input:
        trees = rules.collect_all_iqtree_eval_trees.output.all_eval_trees,
        logs = rules.collect_all_iqtree_eval_logs.output.all_eval_logs,
    output:
        best_tree = f"{full_file_path_iqtree}.bestEvalTreeOfRun"
    params:
        program = "iqtree"
    script:
        "scripts/save_best_tree.py"