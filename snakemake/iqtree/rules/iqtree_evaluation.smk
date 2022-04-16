model = config["parameters"]["model"]["iqtree"]
partitioned = "/" in model

rule re_eval_best_iqtree_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_iqtree_pars}.treefile",
    output:
        iqtree_done = touch(f"{full_file_path_iqtree_eval}.done"),
    params:
        model           = model,
        model_str       = "-p" if partitioned else "-m",
        threads         = config["parameters"]["iqtree"]["threads"],
        prefix_tmp      = full_file_path_iqtree_eval_tmp,
        tmp_dir         = full_dir_iqtree_eval_tmp,
        eval_log        = f"{full_file_path_iqtree_eval_tmp}.iqtree.eval.log",
    log:
        f"{full_file_path_iqtree_eval}.snakelog"
    shell:
        "mkdir -p {params.tmp_dir}; "
        "{iqtree_command} "
        "{params.model_str} {params.model} "
        "-s {input.msa} "
        "-te {input.best_tree_of_run} "
        "-pre {params.prefix_tmp} "
        "-nt {params.threads} "
        ">> {params.eval_log} "

rule reveal_hidden_eval_files:
    input: 
        rules.re_eval_best_iqtree_tree.output.iqtree_done # this is just a dummy and is unused
    output:
        iqtree_log = f"{full_file_path_iqtree_eval}.iqtree.eval.log",
        treefile   = f"{full_file_path_iqtree_eval}.treefile",
    params:
        prefix_tmp = full_dir_iqtree_eval_tmp,
        prefix = full_dir_iqtree_eval
    shell:
        "mv {params.prefix_tmp}/* {params.prefix}"