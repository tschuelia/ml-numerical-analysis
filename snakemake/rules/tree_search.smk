rule raxml_tree:
    input:
        msa             = config["data"]["input"],
    output:
        raxml_best_tree       = f"{full_file_path_raxml}.raxml.bestTree",
        raxml_best_model      = f"{full_file_path_raxml}.raxml.bestModel",
        raxml_ml_trees        = f"{full_file_path_raxml}.raxml.mlTrees",
    params:
        # general params
        model           = config["parameters"]["model"]["raxml-ng"],
        # raxml-ng specific params
        threads         = config["parameters"]["raxml-ng"]["threads"],
        seed            = config["parameters"]["raxml-ng"]["seed"],
        num_pars_trees  = config["parameters"]["raxml-ng"]["num_pars_trees"],
        num_rand_trees  = config["parameters"]["raxml-ng"]["num_rand_trees"],
        prefix          = full_file_path_raxml,
    resources:
        mem_mb=10000 # 10GB (in MB)
        runtime=1440 # 1 day (in minutes)
    log:
        raxml_treesearch_log    = f"{full_file_path_raxml}.raxml.treesearch.log",
    shell:
        "{raxml_command} " 
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {wildcards.blmin} "
        "--blmax {wildcards.blmax} "
        "--threads {params.threads} "
        "--seed {params.seed} "
        "--tree pars{{{params.num_pars_trees}}},rand{{{params.num_rand_trees}}} "
        "> {log} "