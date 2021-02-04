rule iqtree_pars_tree:
    input:
        msa = config["data"]["input"],
    output:
        iqtree_done = touch(f"{full_file_path_iqtree_pars}.done"),
    params:
        model       = config["parameters"]["model"]["iqtree"],
        threads     = config["parameters"]["iqtree"]["threads"],
        prefix_tmp  = full_file_path_iqtree_pars_tmp,
        tmp_dir     = full_dir_iqtree_pars_tmp,
        search_log  = f"{full_file_path_iqtree_pars_tmp}.iqtree.treesearch.log"
    log:
        f"{full_file_path_iqtree_pars}.snakelog"
    shell:
        "mkdir -p {params.tmp_dir}; "
        "{iqtree_command} "
        "-m {params.model} "
        "-s {input.msa} "
        "-ninit 1 "
        "-blmin {wildcards.blmin} "
        "-blmax {wildcards.blmax} "
        "-seed {wildcards.seed} "
        "-pre {params.prefix_tmp} "
        "-nt {params.threads} "
        ">> {params.search_log} "

rule reveal_hidden_treesearch_log:
    input: 
        rules.iqtree_pars_tree.output.iqtree_done # this is just a dummy and is unused
    output:
        iqtree_log = f"{full_file_path_iqtree_pars}.iqtree.treesearch.log",
        treefile   = f"{full_file_path_iqtree_pars}.treefile"
    params:
        prefix_tmp = full_dir_iqtree_pars_tmp,
        prefix = full_dir_iqtree_pars
    shell:
        "cp {params.prefix_tmp}/* {params.prefix}"


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
        iqtree_done = touch(f"{full_file_path_iqtree_eval}.done"),
    params:
        model           = config["parameters"]["model"]["iqtree"],
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix_tmp      = full_file_path_iqtree_eval_tmp,
        tmp_dir         = full_dir_iqtree_eval_tmp,
        eval_log        = f"{full_file_path_iqtree_eval_tmp}.iqtree.eval.log"
    log:
        f"{full_file_path_iqtree_eval}.snakelog"
    shell:
        "mkdir -p {params.tmp_dir}; "
        "{iqtree_command} "
        "-m {params.model} "
        "-s {input.msa} "
        "-te {input.best_tree_of_run} "
        "-blmin {wildcards.blmin_eval} "
        "-blmax {wildcards.blmax_eval} "
        "-pre {params.prefix_tmp} "
        "-nt {params.threads} "
        ">> {params.eval_log} "

rule reveal_hidden_eval_log:
    input: 
        rules.re_eval_best_iqtree_tree.output.iqtree_done # this is just a dummy and is unused
    output:
        iqtree_log = f"{full_file_path_iqtree_eval}.iqtree.eval.log",
        treefile   = f"{full_file_path_iqtree_eval}.treefile",
    params:
        prefix_tmp = full_dir_iqtree_eval_tmp,
        prefix = full_dir_iqtree_eval
    shell:
        "cp {params.prefix_tmp}/* {params.prefix}"

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