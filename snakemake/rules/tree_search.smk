rule raxml_tree:
    input:
        msa             = config["data"]["input"],
    output:
        best_tree       = f"{full_file_path_raxml}.raxml.bestTree",
        best_model      = f"{full_file_path_raxml}.raxml.bestModel",
        ml_trees        = f"{full_file_path_raxml}.raxml.mlTrees",
    params:
        # general params
        model           = config["parameters"]["model"]["raxml-ng"],
        # raxml-ng specific params
        threads         = config["parameters"]["raxml-ng"]["threads"],
        seed            = config["parameters"]["raxml-ng"]["seed"],
        num_pars_trees  = config["parameters"]["raxml-ng"]["num_pars_trees"],
        num_rand_trees  = config["parameters"]["raxml-ng"]["num_rand_trees"],
        prefix          = full_file_path_raxml,
    log:
        f"{full_file_path_raxml}.raxml.treesearch.log",
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