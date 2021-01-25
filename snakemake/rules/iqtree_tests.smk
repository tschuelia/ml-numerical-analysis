rule perform_iqtree_tests:
    input: 
        msa         = config["data"]["input"],
        all_trees    = f"{full_file_path_raxml}.allTreesCollected",
        best_tree   = f"{full_file_path_raxml}.bestTreeOfRun",
    output:
        summary     = f"{full_file_path_iqtree}.iqtree",
        best_iqtree = f"{full_file_path_iqtree}.treefile",
        all_trees   = f"{full_file_path_iqtree}.trees",
    params:
        model       = config["parameters"]["model"]["iqtree"],
        threads     = config["parameters"]["iqtree"]["threads"],
        prefix      = full_file_path_iqtree,
    log:
        f"{full_file_path_iqtree}.iqtree_tests.log"
    shell:
        "{iqtree_command} "
        "-blmin {wildcards.blmin} "
        "-s {input.msa} "
        "-m {params.model} "
        "-pre {params.prefix} "
        "-z {input.all_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "> {log} "



