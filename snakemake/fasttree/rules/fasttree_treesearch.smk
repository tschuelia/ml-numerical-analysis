rule fasttree_tree:
    input:
        msa = config["data"]["input"],
    output:
        fasttree_done = touch(f"{full_file_path_fasttree_tree}.done"),
    params:
        model       = config["parameters"]["model"]["fasttree"],
        prefix_tmp  = full_file_path_fasttree_tree_tmp,
        tmp_dir     = full_dir_fasttree_tree_tmp,
        search_log  = f"{full_file_path_fasttree_tree_tmp}.fasttree.treesearch.log"
    log:
        f"{full_file_path_fasttree_tree}.snakelog"
    shell:
        "mkdir -p {params.tmp_dir}; "
        "{fasttree_command} "
        "-{params.model} "
        "-gamma "
        "-nt "
        "-seed {wildcards.seed} "
        "-blmin {wildcards.blmin} "
        "-lheps {wildcards.lh_eps} "
        "< {input.msa} "
        "> {params.prefix_tmp}.treefile "
        "2> {params.search_log} "


rule reveal_hidden_fasttree_treesearch_files:
    input:
        rules.fasttree_tree.output.fasttree_done # this is just a dummy and is unused
    output:
        iqtree_log = f"{full_file_path_fasttree_tree}.fasttree.treesearch.log",
        treefile   = f"{full_file_path_fasttree_tree}.treefile"
    params:
        prefix_tmp = full_dir_fasttree_tree_tmp,
        prefix = full_dir_fasttree_tree
    shell:
        "cp {params.prefix_tmp}/* {params.prefix}"

