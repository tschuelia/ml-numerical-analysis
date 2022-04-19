model = config["parameters"]["model"]["iqtree"]
partitioned = "/" in model

rule iqtree_pars_tree:
    input:
        msa = config["data"]["input"],
    output:
        iqtree_done = touch(f"{full_file_path_iqtree_pars}.done"),
    params:
        model       = config["parameters"]["model"]["iqtree"],
        model_str   = "-p" if partitioned else "-m",
        threads     = config["parameters"]["iqtree"]["threads"],
        prefix_tmp  = full_file_path_iqtree_pars_tmp,
        tmp_dir     = full_dir_iqtree_pars_tmp,
        search_log  = f"{full_file_path_iqtree_pars_tmp}.iqtree.treesearch.log"
    log:
        f"{full_file_path_iqtree_pars}.snakelog"
    shell:
        "mkdir -p {params.tmp_dir}; "
        "{iqtree_command} "
        "{params.model_str} {params.model} "
        "-s {input.msa} "
        "-ninit 1 "
        "-eps {wildcards.lh_eps} "
        "-seed {wildcards.seed} "
        "-pre {params.prefix_tmp} "
        "-nt {params.threads} "
        ">> {params.search_log} "

rule reveal_hidden_treesearch_files:
    input: 
        rules.iqtree_pars_tree.output.iqtree_done # this is just a dummy and is unused
    output:
        iqtree_log = f"{full_file_path_iqtree_pars}.iqtree.treesearch.log",
        treefile   = f"{full_file_path_iqtree_pars}.treefile"
    params:
        prefix_tmp = full_file_path_iqtree_pars_tmp, #full_dir_iqtree_pars_tmp,
        prefix = full_dir_iqtree_pars
    shell:
        "cp {params.prefix_tmp}* {params.prefix}"

